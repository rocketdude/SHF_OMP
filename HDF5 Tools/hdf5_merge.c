 /*@@
   @file      hdf5_merge.c
   @date      Thu 10 Jan 2002
   @author    Thomas Radke
   @desc
              This utility program takes a list of Cactus HDF5 datafiles,
              merges them at the group hierarchy level and dumps the resulting
              tree to a new HDF5 file.
              Compile using: gcc -O3 -lhdf5 -o hdf5_merge hdf5_merge.c
   @enddesc
   @version   $Id: hdf5_merge.c,v 1.2 2008/01/23 13:21:42 tradke Exp $
 @@*/

//#include "cctk.h"

// some macros to fix compatibility issues as long
// as 1.8.0 is in beta phase
#define H5_USE_16_API 1

#include <hdf5.h>
#include <search.h>

#if (H5_VERS_MAJOR == 1 && (H5_VERS_MINOR == 8) && (H5_VERS_RELEASE == 0))
#warning "Hacking HDF5 1.8.0 compatiblity with 1.6.x; fix once 1.8.0 stable"
#endif

#include <time.h> 
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define DEBUG

/* the rcs ID and its dummy function to use it */
//static const char *rcsid = "$Header: /cactusdevcvs/CactusExternal/HDF5/src/util/hdf5_merge.c,v 1.2 2008/01/23 13:21:42 tradke Exp $";
//CCTK_FILEVERSION(CactusExternal_HDF5_util_hdf5_merge_c)


/*****************************************************************************/
/*                           macro definitions                               */
/*****************************************************************************/
/* macro to do an HDF5 call, check its return code, and print a warning
   in case of an error */
#define CHECK_ERROR(hdf5_call)                                                \
          do                                                                  \
          {                                                                   \
            int _error_code = hdf5_call;                                      \
                                                                              \
                                                                              \
            if (_error_code < 0)                                              \
            {                                                                 \
              fprintf (stderr, "WARNING: line %d: HDF5 call '%s' returned "   \
                               "error code %d\n",                             \
                                __LINE__, #hdf5_call, _error_code);           \
              nerrors++;                                                      \
            }                                                                 \
          } while (0)


/*****************************************************************************/
/*                           global variables                                */
/*****************************************************************************/
/* NOTE: although it isn't good programming practice
         we make these variables global for convenience
         since they are accessed from recursively or
         indirectly called routines which only get passed
         a single user-supplied argument */
static char *current_file  = NULL; /* name of the current input file */
static char *pathname = NULL;      /* pathname of the current object */
static unsigned int nerrors = 0;   /* global error counter */
static int do_create, do_copy;     /* for two-pass runs, what to do in this pass */
static int verbose = 0;            /* output dataset names as we go */
static int create_groups = 0;      /* create one group per iteration */
static int two_passes = 0;         /* first create all groups and datasets, then fill them with data */
static int change_nioprocs = 0;    /* set Cactus' nioprocs */
static time_t program_startup;     /* used to determine if a dataset has been created by us */

struct treenode
{
  int iteration;
  const char *objectname;
  const char *source_file;
};
static void *iterations_seen = NULL;

/*****************************************************************************/
/*                           local function prototypes                       */
/*****************************************************************************/
static herr_t CopyObject (hid_t copy_from, const char *objectname, void *arg);
static herr_t CopyAttribute (hid_t src, const char *attr_name, void *arg);
static int CheckIfExists(hid_t from, const char *objectname);
static void usage(char *argv[]);
static struct treenode *find_iteration(int iteration, const char *objectname);
static struct treenode *add_iteration(int iteration, const char *objectname, const char *filename);
static int compare(const void *a, const void *b);
static int files_from_same_set(const char *a, const char *b);

static int compare(const void *a, const void *b)
{
  struct treenode *n, *m;
  const char *s, *t;

  n = (struct treenode *)a;
  m = (struct treenode *)b;

  if(n->iteration != m->iteration)
    return n->iteration - m->iteration;

  // compare beginning of strings stopping at ' ' or at the end of the string
  s = n->objectname, t = m->objectname;
  while(*s == *t && *s != '\0' && *s != ' ')
    s++,t++;
  return *s - *t;
}

static struct treenode *find_iteration(int iteration, const char *objectname)
{
  struct treenode n = {iteration, objectname, NULL};
  struct treenode **m;

  /* tfind and tsearch return a pointer to a structure whose first element is a
   * pointer to the payload */
  m = tfind(&n, &iterations_seen, compare);

  return m ? *m : NULL;
}

static struct treenode *add_iteration(int iteration, const char *objectname, const char *fn)
{
  struct treenode *n;

  n = find_iteration(iteration, objectname);
  if(NULL == n)
  {
    n = malloc(sizeof(*n));
    assert(n != NULL);

    n->iteration = iteration;
    n->source_file = strdup(fn);
    n->objectname = strdup(objectname);

    tsearch(n, &iterations_seen, compare);
  }

  return n;
}

/* checks if files *a and *b are from the same set of output files */
static int files_from_same_set(const char *a, const char *b)
{
  const char *ext_a, *ext_b;

  ext_a = strstr(a, ".file_");
  ext_b = strstr(b, ".file_");

  if(ext_a && ext_b)
  {
    /* both are chunked, the part before and after file_ must agree */
    if(strncmp(a, b, ext_a-a+1) != 0)
      return 0;
    /* skip beyond .file_[0-9]*. */
    ext_a = strchr(ext_a+1, '.');
    ext_b = strchr(ext_b+1, '.');
    if(!ext_a || !ext_b)
      return 0;
    /* compare remainder */
    return strcmp(ext_a, ext_b) == 0;
  }
  else if(!ext_a && !ext_b)
  {
    /* both are unchunked so they must have the very same name */
    return strcmp(a, b) == 0;
  }
  else
  {
    /* one file is chunked the other unchunked, so they can never be from the same set */
    return 0;
  }
}


 /*@@
   @routine    main
   @date       Sat 24 Feb 2001
   @author     Thomas Radke
   @desc
               Main routine of the HDF5 file merger
   @enddesc

   @calls      CopyObject

   @var        argc
   @vdesc      number of command line arguments
   @vtype      int
   @vio        in
   @endvar
   @var        argv
   @vdesc      command line arguments
   @vtype      char *[]
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 for success, negative return values indicate an error
   @endreturndesc
@@*/
int main (int argc, char *argv[])
{
  int i, pass;
  hid_t *infiles, outfile;

  /* used to track which objects are created by me */
  program_startup = time(NULL);

  /* parse options */
  for (i = 1 ; i < argc ; i++)
  {
    int j;

    /* check if this option is valid */
    if(strlen(argv[i]) != 2 || *argv[i] != '-') /* not an option after all */
      continue;

    /* parse option */
    switch(argv[i][1])
    {
      case 'v':
        verbose = 1;
        break;
      case 'g':
        create_groups = 1;
        break;
      case 't':
        two_passes = 1;
        break;
      case 'u':
         change_nioprocs = 1;
        break;
      default:
        fprintf(stderr, "unknown option '%s'.", argv[i]);
        usage(argv);
        return (1);
        break;
    }

    /* remove option from list of arguments */
    for(j = i+1 ; j < argc ; j++)
      argv[j-1] = argv[j];
    argc -= 1;
    i -= 1; /* re-parse the current option */
  }

  /* give some help if called with incorrect number of parameters */
  if (argc < 3)
  {
    usage(argv);
    return (0);
  }

  H5E_BEGIN_TRY
  {
    /* open the input files */
    infiles = (hid_t *) malloc ((argc - 2) * sizeof (hid_t));
    for (i = 0; i < argc - 2; i++)
    {
      infiles[i] = H5Fopen (argv[i + 1], H5F_ACC_RDONLY, H5P_DEFAULT);
      if (infiles[i] < 0)
      {
        fprintf (stderr, "ERROR: Cannot open HDF5 input file '%s' !\n\n",
                 argv[i + 1]);
        return (1);
      }
      H5Fclose (infiles[i]);
    }

    /* try to open an existing outfile file in append mode,
       if this fails create it as a new file */
    outfile = H5Fopen (argv[argc-1], H5F_ACC_RDWR, H5P_DEFAULT);
    if (outfile < 0)
    {
      outfile = H5Fcreate (argv[argc-1], H5F_ACC_TRUNC, H5P_DEFAULT,
                           H5P_DEFAULT);
    }
    if (outfile < 0)
    {
      fprintf (stderr, "ERROR: Cannot open HDF5 output file '%s' !\n\n",
               argv[argc - 1]);
      return (1);
    }
  } H5E_END_TRY

  /*printf ("\n  -------------------------\n"
            "  Cactus 4 HDF5 File Merger\n"
            "  -------------------------\n");*/

  /* do the copying by iterating over all objects */
  do_create = 1;
  do_copy = !two_passes;
  for (pass = 0 ; pass < (two_passes ? 2 : 1) ; pass++)
  {
    for (i = 0; i < argc - 2; i++)
    {
      if(two_passes)
      {
        printf ("\n  Merging objects from input file '%s' into output file '%s' (pass %d)\n",
                argv[i + 1], argv[argc-1], pass+1);
      }
      else
      {
        /*printf ("\n  Merging objects from input file '%s' into output file '%s'\n",
                argv[i + 1], argv[argc-1]);*/
      }
      pathname = (char *)"/";
      current_file = argv[i + 1];
      CHECK_ERROR (infiles[i] = H5Fopen (current_file, H5F_ACC_RDONLY, H5P_DEFAULT));
      CHECK_ERROR (H5Giterate (infiles[i], "/", NULL, CopyObject, &outfile));
      /* finally, close the open file */
      CHECK_ERROR (H5Fclose (infiles[i]));
    }

    do_create = 0;
    do_copy = 1;
  }

  CHECK_ERROR (H5Fclose (outfile));
  free (infiles);

  /* report status */
  if (nerrors == 0)
  {
    /*printf ("\n\n   *** All input files successfully merged. ***\n\n");*/
  }
  else
  {
    fprintf (stderr, "\n\n   *** WARNING: %d errors occured during "
                     "file merging. ***\n\n", nerrors);
  }

  return (0);
}


/*****************************************************************************/
/*                           local routines                                  */
/*****************************************************************************/

 /*@@
   @routine    usage
   @date       Wed Oct 28 12:03:18 EDT 2009
   @author     Roland Haas
   @desc
               Print usage information of the HDF5 file merger
   @enddesc

   @calls      fprintf

   @var        argv
   @vdesc      command line arguments
   @vtype      char *[]
   @vio        in
   @endvar

   @returntype void
@@*/
void usage(char *argv[])
{
    fprintf (stderr, "Usage: %s [-g] [-t] [-v] <infile1> [<infile2> ...] <outfile>\n",argv[0]);
    fprintf (stderr, "       -g : create groups for each iteration\n");
    fprintf (stderr, "       -u : set niprocs=1 (unchunk)\n");
    fprintf (stderr, "       -t : copy datasets in two passes\n");
    fprintf (stderr, "       -v : output each dataset name as it is copied\n");
    fprintf (stderr, "       Cactus' hdf5_merge uses -v by default\n");
    fprintf (stderr, "   eg, %s -g -t alp.time*.h5 alp.h5"
                     "\n\n", argv[0]);
}
 /*@@
   @routine    CheckIfExists
   @date       Fri Oct 30 16:35:31 EDT 2009
   @author     Roland Haas
   @desc
               Check if an object of a given name exists
   @enddesc

   @calls      H5Gget_objinfo

   @var        from
   @vdesc      identifier for the group the object belongs to
   @vtype      hid_t
   @vio        in
   @endvar
   @var        objectname
   @vdesc      name of the object
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 - object does not exist
               1 - object exists
   @endreturndesc
@@*/
static int CheckIfExists(hid_t from, const char *objectname)
{
  int retval;
  H5G_stat_t objectinfo;

  /* check whether an object by that name already exists */
  H5E_BEGIN_TRY
  {
    retval = H5Gget_objinfo (from, objectname, 0, &objectinfo) >= 0;
  } H5E_END_TRY

  return retval;
}

 /*@@
   @routine    CopyObject
   @date       Sat 24 Feb 2001
   @author     Thomas Radke
   @desc
               Iterator recursively called by H5Giterate() for every object
               in the input file
               It copies the current object to the output file if it didn't
               already exist there.
   @enddesc

   @calls      CopyAttribute

   @var        from
   @vdesc      identifier for the group the current object belongs to
   @vtype      hid_t
   @vio        in
   @endvar
   @var        objectname
   @vdesc      name of the current object
   @vtype      const char *
   @vio        in
   @endvar
   @var        _to
   @vdesc      user-supplied argument indicating the output object identifier
   @vtype      hid_t
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 - continue the iteration for following group objects
               1 - short-curcuit, no further iteration of this group
   @endreturndesc
@@*/
static herr_t CopyObject (hid_t from,
                          const char *objectname,
                          void *_to)
{
  hid_t to, datatype, dataspace;
  H5G_stat_t objectinfo;
  char *current_pathname;
  size_t objectsize;
  void *data;

  /* get the output object identifier */
  to = *(hid_t *) _to;

  /* check the type of the current object */
  CHECK_ERROR (H5Gget_objinfo (from, objectname, 0, &objectinfo));
  if (objectinfo.type == H5G_GROUP)
  {
    if (verbose)
      printf ("   iterating through group '%s%s'\n", pathname, objectname);

    /* build the full pathname for the current to object to process */
    current_pathname = pathname;
    pathname = (char *) malloc (strlen (current_pathname) +
                                strlen (objectname) + 2);
    sprintf (pathname, "%s%s/", current_pathname, objectname);

    /* iterate over all objects in the (first) input file */
    CHECK_ERROR (from = H5Gopen (from, objectname));
    if(do_create && !CheckIfExists(to, objectname))
    {
      CHECK_ERROR (to = H5Gcreate (to, objectname, 0));
      CHECK_ERROR (H5Aiterate (from, NULL, CopyAttribute, &to));
    }
    else
    {
      CHECK_ERROR (to = H5Gopen (to, objectname));
    }
    CHECK_ERROR (H5Giterate (from, ".", NULL, CopyObject, &to));
    CHECK_ERROR (H5Gclose (to));
    CHECK_ERROR (H5Gclose (from));

    /* reset the pathname */
    free (pathname);
    pathname = current_pathname;
  }
  else if (objectinfo.type == H5G_DATASET)
  {

    /* create a group for the timestep if it does not yet exist */
    char timestep[128]; /* uses to hold the string it123456789 */
    struct treenode *n;
    char *p;
    int iter;
    static hid_t tsgroup = -1;
    static int cached_iter = -1;
    
    p=strstr(objectname, "it=");
    if(NULL != p)
    {
      /* mark where we first saw this iteration and only accept further datasets
       * for it from there */
      iter = atoi(p+3); // skip 'it=', atoi stops at non-digit characters
      n = add_iteration(iter, objectname, current_file);
    }
    else
    {
      iter = -1;
      n = NULL;
    }

    /* copy dataset only if it is from the correct file or if it does not carry
     * an iteration tag (which likely makes it metadata) */
    if (NULL != n && !files_from_same_set(n->source_file, current_file))
    {
      if (verbose)
      {
        printf ("   object '%s%s' will not be copied (iteration already exists from earlier file)\n",
                 pathname, objectname);
      }
      return (0);
    }

    CHECK_ERROR (from = H5Dopen (from, objectname));
    CHECK_ERROR (datatype = H5Dget_type (from));
    CHECK_ERROR (dataspace = H5Dget_space (from));

    /* put each iteration in a separate group if so requested */
    if (iter != -1 && create_groups)
    {
      /* NB: this can and has to be called during both passes to get the correct
       * to handle */
      if (iter != cached_iter)
      {
        if (tsgroup != -1)
        {
          CHECK_ERROR (H5Gclose (tsgroup));
          tsgroup = -1;
        }

        assert(snprintf(timestep, sizeof(timestep)/sizeof(timestep[0]), "it%09d", iter) < (int)(sizeof(timestep)/sizeof(timestep[0])));
        if (!do_create || CheckIfExists(to, timestep))
          CHECK_ERROR (tsgroup = H5Gopen (to, timestep));
        else
          CHECK_ERROR (tsgroup = H5Gcreate (to, timestep, 0));
        
        cached_iter=iter;
      }
      to = tsgroup;
    }

    /* check if object already exists (and was not created by us) */
    if (do_create && CheckIfExists(to, objectname))
    {
      if (verbose)
      {
        printf ("   object '%s%s' will not be copied (already exists)\n",
                 pathname, objectname);
      }
      return (0);
    } 
    else if (do_copy && !do_create)
    {
      /* we only get here during pass two of two-pass runs */
      CHECK_ERROR (H5Gget_objinfo (to, objectname, 0, &objectinfo));
      if (objectinfo.mtime < program_startup)
      {
          /* object was not created by us, so we don't modify it (see first
           * pass testabove) */
          return (0);
      }
    }

    if (verbose)
    {
      if (do_copy)
        printf ("   copying dataset '%s%s'\n", pathname, objectname);
      else
        printf ("   creating dataset '%s%s'\n", pathname, objectname);
    }

    /* first pass: create datasets */
    if (do_create) 
    {
      CHECK_ERROR (to = H5Dcreate (to, objectname, datatype, dataspace,
                  H5P_DEFAULT));
      CHECK_ERROR (H5Aiterate (from, NULL, CopyAttribute, &to));
    }
    else
    {
      CHECK_ERROR (to = H5Dopen (to, objectname));
    }

    /* second pass: fill with data */
    if (do_copy)
    {
      objectsize = H5Sget_select_npoints (dataspace) * H5Tget_size (datatype);
      if (objectsize > 0)
      {
        data = malloc (objectsize);
        if (data == NULL)
        {
          fprintf (stderr, "failled to allocate %d bytes of memory, giving up\n",
                   (int) objectsize);
          exit (-1);
        }
        CHECK_ERROR (H5Dread (from, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                              data));
        CHECK_ERROR (H5Dwrite (to, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,data));
        free (data);
      }
    }
    CHECK_ERROR (H5Dclose (to));
    CHECK_ERROR (H5Dclose (from));
    CHECK_ERROR (H5Sclose (dataspace));
    CHECK_ERROR (H5Tclose (datatype));
  }
  else
  {
    fprintf (stderr, "WARNING: Found object '%s%s' which is not neither a "
                     "group nor a dataset ! Object will not be copied.\n",
                     pathname, objectname);
    nerrors++;
  }

  return (0);
}


 /*@@
   @routine    CopyAttribute
   @date       Sat 24 Feb 2001
   @author     Thomas Radke
   @desc
               Iterator recursively called by H5Aiterate() for every attribute
               of an object (dataset or group)
   @enddesc

   @var        from
   @vdesc      identifier for the group or dataset to read the attribute from
   @vtype      hid_t
   @vio        in
   @endvar
   @var        attrname
   @vdesc      name of the current attribute
   @vtype      const char *
   @vio        in
   @endvar
   @var        _to
   @vdesc      user-supplied argument indicating the group or dataset
               to copy the attribute to
   @vtype      hid_t
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 - continue the iteration for following attributes
   @endreturndesc
@@*/
static herr_t CopyAttribute (hid_t from,
                             const char *attrname,
                             void *_to)
{
  hid_t attr, datatype, dataspace, to;
  size_t attrsize;
  void *value;


  /* get the target group/dataset identifier */
  to = *(hid_t *) _to;

  /* open the attribute given by its name, get type, dataspace, and value
     and just copy it */
  CHECK_ERROR (attr = H5Aopen_name (from, attrname));
  CHECK_ERROR (datatype = H5Aget_type (attr));
  CHECK_ERROR (dataspace = H5Aget_space (attr));
  attrsize = H5Tget_size (datatype);
  if (H5Sis_simple (dataspace) > 0)
  {
    attrsize *= H5Sget_simple_extent_npoints (dataspace);
  }
  if (attrsize > 0)
  {
    value = malloc (attrsize);
    CHECK_ERROR (H5Aread (attr, datatype, value));
    CHECK_ERROR (H5Aclose (attr));
    if (change_nioprocs &&
        strcmp(attrname, "nioprocs") == 0 &&
        strcmp(pathname, "/Parameters and Global Attributes/") == 0)
    {
        *((int*)value) = 1;
    }
    CHECK_ERROR (attr = H5Acreate (to, attrname, datatype, dataspace,
                                   H5P_DEFAULT));
    CHECK_ERROR (H5Awrite (attr, datatype, value));
    free (value);
  }
  CHECK_ERROR (H5Aclose (attr));
  CHECK_ERROR (H5Sclose (dataspace));
  CHECK_ERROR (H5Tclose (datatype));

  return (0);
}
