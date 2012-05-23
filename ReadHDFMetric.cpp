/*
 *  ReadHDFMetric.cpp
 *  
 *
 *  Created by Darius Bunandar on 5/11/12.
 *  Copyright 2012 The University of Texas at Austin. All rights reserved.
 *
 *	This function is to read metric data files in HDF5
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <cstring>
#include <regex.h>

#include <hdf5.h>
using namespace	std;

/*****************************************************************************
 *************************     Macro Definitions   ***************************
 *****************************************************************************/

// value for an unset parameter
#define PARAMETER_UNSET	-424242.0

// fuzzy factor for comparing dataset timesteps with user-specified value
#define FUZZY_FACTOR	1e-6

// the datatype for the 'start' argument in H5Sselect_hyperslab()
#if (H5_VERS_MAJOR == 1 && \
(H5_VERS_MINOR < 6 || (H5_VERS_MINOR == 6 && H5_VERS_RELEASE < 4)))
#define h5size_t hssize_t
#else
#define h5size_t hsize_t
#endif

// macro to check the return code of calls to the HDF5 library
#define CHECK_HDF5(fn_call) {                                                 \
		  const int _retval = fn_call;                                        \
		  if (_retval < 0) {                                                  \
			cerr << "HDF5 call " << #fn_call                                  \
				 << " in file " << __FILE__ << " line " << __LINE__           \
				 << " returned error code " << _retval << endl;               \
		  }                                                                   \
		 }

// the name of the group containing file metadata
#define METADATA_GROUP "Parameters and Global Attributes"

/*****************************************************************************
 *************************       Global Data         *************************
 *****************************************************************************/

// the slab coordinate as selected by the user
static double slab_coord[3] = {PARAMETER_UNSET, PARAMETER_UNSET, PARAMETER_UNSET};

// flag for outputting data in full 3D
static bool out3D = false;

// the specific timestep selected by the user
static double timestep = PARAMETER_UNSET;

// the regular expression to match against dataset names
static const char* regex = NULL;
static regex_t preg;

// the total number of datasets sliced
static unsigned int slices_extracted = 0;

// whether to print omitted datasets
static bool verbose = false; 

// rank of the output slices
static hsize_t outrank = 0;

// output file id
static hid_t outfile = -1;

// output data buffer
double **buffer;

/*****************************************************************************
 *************************     External directive  ***************************
 *****************************************************************************/
extern "C" {
	int readhdfmetric_(int argc, char *const argv[], double **dataout);
}

/*****************************************************************************
 *************************     Function Prototypes   *************************
 *****************************************************************************/
static herr_t ProcessDataset (hid_t group, const char *name, void *_file);

/*****************************************************************************
 *************************     Function Starts Here  *************************
 *****************************************************************************/

int readhdfmetric_(int argc, char *const argv[], double **dataout)
{
	
	buffer = dataout;
	
	int i;
	bool help = false;
	int num_slab_options = 0;
	
	// evaluate command line parameters
	for (i = 1; i < argc; i++) {
		if (strcmp (argv[i], "--help") == 0) {
			help = true; break;
		} else if (strcmp (argv[i], "--verbose") == 0) {
			verbose = true;
		} else if (strcmp (argv[i], "--out3d-cube") == 0) {
			outrank = 3;
			out3D = true;
		} else if (strcmp (argv[i], "--timestep") == 0 and i+1 < argc) {
			timestep = atof (argv[++i]);
		} else if (strcmp (argv[i], "--match") == 0 and i+1 < argc) {
			regex = argv[++i];
			if (regcomp (&preg, regex, REG_EXTENDED)) {
				cerr << "Error: invalid regular expression '" << regex << "' given"
				<< endl << endl;
				return (-1);
			}
		} else if (strcmp (argv[i], "--out1d-xline-yz") == 0 and i+2 < argc) {
			outrank = 1;
			slab_coord[1] = atof (argv[++i]);
			slab_coord[2] = atof (argv[++i]); num_slab_options++;
		} else if (strcmp (argv[i], "--out1d-yline-xz") == 0 and i+2 < argc) {
			outrank = 1;
			slab_coord[0] = atof (argv[++i]);
			slab_coord[2] = atof (argv[++i]); num_slab_options++;
		} else if (strcmp (argv[i], "--out1d-zline-xy") == 0 and i+2 < argc) {
			outrank = 1;
			slab_coord[0] = atof (argv[++i]);
			slab_coord[1] = atof (argv[++i]); num_slab_options++;
		} else if (strcmp (argv[i], "--out2d-yzplane-x") == 0 and i+1 < argc) {
			outrank = 2;
			slab_coord[0] = atof (argv[++i]); num_slab_options++;
		} else if (strcmp (argv[i], "--out2d-xzplane-y") == 0 and i+1 < argc) {
			outrank = 2;
			slab_coord[1] = atof (argv[++i]); num_slab_options++;
		} else if (strcmp (argv[i], "--out2d-xyplane-z") == 0 and i+1 < argc) {
			outrank = 2;
			slab_coord[2] = atof (argv[++i]); num_slab_options++;
		} else {
			break;
		}
	}
	
	/* give some help if called with incorrect number of parameters */
	if (help or i >= argc-1 or num_slab_options != (out3D ? 0 : 1)) {
		const string indent (strlen (argv[0]) + 1, ' ');
		cerr << endl << "   ------------------"
		<< endl << "   Carpet HDF5 Slicer"
		<< endl << "   ------------------" << endl
		<< endl
		<< "Usage: " << endl
		<< argv[0] << " [--help]" << endl
		<< indent << "[--match <regex string>]" << endl
		<< indent << "[--timestep <cctk_time value>]" << endl
		<< indent << "[--verbose]" << endl
		<< indent << "<--out1d-line value value> | <--out2d-plane value> | <out3d-cube>" << endl
		<< indent << "<hdf5_infiles> <hdf5_outfile>" << endl << endl
		<< "  where" << endl
		<< "    [--help]                         prints this help" << endl
		<< "    [--match <regex string>]         selects HDF5 datasets by their names" << endl
		<< "                                     matching a regex string using POSIX" << endl
		<< "                                     Extended Regular Expression syntax" << endl
		<< "    [--timestep <cctk_time value>]   selects all HDF5 datasets which" << endl
		<< "                                     (fuzzily) match the specified time" << endl
		<< "    [--verbose]                      lists skipped HDF5 datasets on stderr" << endl
		<< endl
		<< "  and either <--out1d-line value value> or <--out2d-plane value> or <--out3d-cube>" << endl
		<< "  must be specified as in the following:" << endl
		<< "    --out1d-xline-yz  <origin_y> <origin_z>" << endl
		<< "    --out1d-yline-xz  <origin_x> <origin_z>" << endl
		<< "    --out1d-zline-xy  <origin_x> <origin_y>" << endl
		<< endl
		<< "    --out2d-yzplane-x  <origin_x>" << endl
		<< "    --out2d-xzplane-y  <origin_y>" << endl
		<< "    --out2d-xyplane-z  <origin_z>" << endl
		<< endl
		<< "    --out3d-cube" << endl
#if 0
		<< "    --out2d-yzplane-xi <origin_xi>" << endl
		<< "    --out2d-xzplane-yi <origin_yi>" << endl
		<< "    --out2d-xyplane-zi <origin_zi>" << endl
#endif
		<< endl
		<< "  eg, to extract a 2D xy-plane at z = 0:" << endl
		<< "    " << argv[0] << " --out2d-xyplane-z 0 alp.file_*.h5 alp.z=0.h5" << endl
		<< endl
		<< "  or the same plane but only for datasets of refinement level 0:" << endl
		<< "    " << argv[0] << " --match 'ADMBASE::alp it=[0-9]+ tl=0 rl=0' --out2d-xyplane-z 0 alp.file_*.h5 alp.z=0.h5" << endl;

		return (0);
	}
	
	// don't use output file, here we output to variable
	// so the last variable of argv just becomes a dummy
	
	// browse though input file(s)
	vector<hid_t> filelist;
	for (; i < argc-1; i++) {
		hid_t file;
		
		H5E_BEGIN_TRY {
			file = H5Fopen (argv[i], H5F_ACC_RDONLY, H5P_DEFAULT);
		} H5E_END_TRY;
		
		if (file < 0) {
			cerr << "Could not open input file '" << argv[i] << "'" << endl << endl;
			continue;
		}
		
		cout << "  iterating through input file '" << argv[i] << "'..." << endl;
		CHECK_HDF5 (H5Giterate (file, "/", NULL, ProcessDataset, &file));
		
		filelist.push_back (file);
	}
	
	if (slices_extracted == 0) {
		cerr << endl << "No matching datasets were found for the selected "
		"slice parameters." << endl << endl;
	}
	
	// close all files
	for (size_t j = 0; j < filelist.size(); j++) {
		if (filelist[j] >= 0) {
			CHECK_HDF5 (H5Fclose (filelist[j]));
		}
	}
	
	cout << endl << "Done." << endl << endl;
	
	return (0);
}


/*@@
 @routine    ProcessDataset
 @date       Fri 17 October 2008
 @author     Thomas Radke
 @desc
 The worker routine which is called by H5Giterate().
 It checks whether the current HDF5 object is a dataset matching
 the user's slab criteria.
 @enddesc
 
 @var        group
 @vdesc      HDF5 object to start the iteration
 @vtype      hid_t
 @vio        in
 @endvar
 @var        datasetname
 @vdesc      name of the object at the current iteration
 @vtype      const char *
 @vio        in
 @endvar
 @var        _file
 @vdesc      pointer to the descriptor of the currently opened file
 @vtype      void *
 @vio        in
 @endvar
 @@*/
static herr_t ProcessDataset (hid_t group, const char *datasetname, void *_file)
{
	// we are interested in datasets only - skip anything else
	H5G_stat_t object_info;
	CHECK_HDF5 (H5Gget_objinfo (group, datasetname, 0, &object_info));
	if (object_info.type != H5G_DATASET) {
		return (0);
	}
	
	hid_t dataset, datatype, attr;
	CHECK_HDF5 (dataset   = H5Dopen (group, datasetname));
	
	// check the dataset's datatype - make sure it is either integer or real
	CHECK_HDF5 (datatype  = H5Dget_type (dataset));
	H5T_class_t typeclass;
	CHECK_HDF5 (typeclass = H5Tget_class (datatype));
	
	// read the timestep and variable name
	CHECK_HDF5 (attr = H5Aopen_name (dataset, "time"));
	double time;
	CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_DOUBLE, &time));
	CHECK_HDF5 (H5Aclose (attr));
	CHECK_HDF5 (attr = H5Aopen_name (dataset, "name"));
	hid_t stringdatatype;
	CHECK_HDF5 (stringdatatype = H5Aget_type (attr));
	string varname(H5Tget_size (stringdatatype) + 1, 0);
	CHECK_HDF5 (H5Aread (attr, stringdatatype, &varname[0]));
	CHECK_HDF5 (H5Aclose (attr));
	
	// read the dimensions
	hid_t dataspace = -1;
	CHECK_HDF5 (dataspace = H5Dget_space (dataset));
	int ndims;
	ndims = H5Sget_simple_extent_ndims (dataspace);
	hsize_t dims[3];
	int iorigin[3];
	double origin[3], delta[3];
	int itimestep = 0, level = 0;
	
	bool is_okay = false;
	if (typeclass != H5T_FLOAT and typeclass != H5T_INTEGER) {
		if (verbose) {
			cerr << "skipping dataset '" << datasetname << "':" << endl
			<< "  is not of integer or floating-point datatype" << endl;
		}
	} else if (ndims != 3) {
		if (verbose) {
			cerr << "skipping dataset '" << datasetname << "':" << endl
			<< "  dataset has " << ndims << " dimensions" << endl;
		}
	} else if (regex && regexec (&preg, datasetname, 0, NULL, 0)) {
		if (verbose) {
			cerr << "skipping dataset '" << datasetname << "':" << endl
			<< "  name doesn't match regex" << endl;
		}
	} else if (timestep != PARAMETER_UNSET and
			   fabs (timestep - time) > FUZZY_FACTOR) {
		if (verbose) {
			cerr << "skipping dataset '" << datasetname << "':" << endl
			<< "  timestep (" << time << ") doesn't match" << endl;
		}
	} else {
		CHECK_HDF5 (attr = H5Aopen_name (dataset, "iorigin"));
		CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_INT, iorigin));
		CHECK_HDF5 (H5Aclose (attr));
		CHECK_HDF5 (attr = H5Aopen_name (dataset, "origin"));
		CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_DOUBLE, origin));
		CHECK_HDF5 (H5Aclose (attr));
		CHECK_HDF5 (attr = H5Aopen_name (dataset, "delta"));
		CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_DOUBLE, delta));
		CHECK_HDF5 (H5Aclose (attr));
		CHECK_HDF5 (attr = H5Aopen_name (dataset, "timestep"));
		CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_INT, &itimestep));
		CHECK_HDF5 (H5Aclose (attr));
		CHECK_HDF5 (attr = H5Aopen_name (dataset, "level"));
		CHECK_HDF5 (H5Aread (attr, H5T_NATIVE_INT, &level));
		CHECK_HDF5 (H5Aclose (attr));
		CHECK_HDF5 (H5Sget_simple_extent_dims (dataspace, dims, NULL));
		
		int i;
		is_okay = out3D;
		if (not out3D) {
			for (i = 0; i < 3; i++) {
				if (slab_coord[i] != PARAMETER_UNSET) {
					is_okay = slab_coord[i] >= origin[i] and
                    slab_coord[i] <= origin[i] + (dims[2-i]-1)*delta[i];
					if (not is_okay) break;
				}
			}
		}
		if (not is_okay) {
			if (verbose) {
				cerr << "skipping dataset '" << datasetname << "':" << endl
				<< "  slab " << slab_coord[i] << " is out of dataset range ["
				<< origin[i] << ", "
				<< origin[i] + (dims[2-i]-1)*delta[i] << "]"
				<< endl;
			}
		}
	}
	
	if (not is_okay) {
		CHECK_HDF5 (H5Tclose (stringdatatype));
		CHECK_HDF5 (H5Tclose (datatype));
		CHECK_HDF5 (H5Dclose (dataset));
		return (0);
	}
	
	slices_extracted++;
	
	h5size_t slabstart[3] = {0, 0, 0};
	hsize_t slabcount[3] = {dims[0], dims[1], dims[2]};
	hsize_t outslabcount[3];
	double slice_origin[3], slice_delta[3];
	int slice_iorigin[3];
	int j = 0;
	for (int i = 0; i < 3; i++) {
		if (slab_coord[i] != PARAMETER_UNSET) {
			slabstart[2-i] = (h5size_t) ((slab_coord[i] - origin[i]) / delta[i] + 0.5);
			slabcount[2-i] = 1;
		} else {
			outslabcount[outrank-j-1] = dims[2-i];
			slice_origin[j] = origin[i];
			slice_delta[j] = delta[i];
			slice_iorigin[j] = iorigin[i];
			j++;
		}
	}
	assert(j == outrank);
	
	hid_t slabspace;
	CHECK_HDF5 (slabspace = H5Screate_simple (3, slabcount, NULL));
	CHECK_HDF5 (H5Sselect_hyperslab (dataspace, H5S_SELECT_SET,
									 slabstart, NULL, slabcount, NULL));
	const hssize_t npoints = H5Sget_select_npoints (dataspace);
	// make sure the vector allocates at least one element
	char *data = new char[(npoints + 1) * H5Tget_size(datatype)];
	CHECK_HDF5 (H5Dread (dataset, datatype, slabspace, dataspace, H5P_DEFAULT,
						 data));
	CHECK_HDF5 (H5Dclose (dataset));
	CHECK_HDF5 (H5Sclose (slabspace));
	CHECK_HDF5 (H5Sclose (dataspace));
	
	//	*buffer = new double[outslabcount[0]];
	//	for(int ii=0;ii< outslabcount[0]; ii++) buffer[ii] = new double [outslabcount[1]];
	
	int kk = 0; 
	for(int ii=0; ii<outslabcount[0]; ii++){
		for(int jj=0; jj<outslabcount[1]; jj++){
			buffer[jj][ii]=data[kk];
			kk++;
		}
	}
	
	delete[] data;
	
	return 0;
}
