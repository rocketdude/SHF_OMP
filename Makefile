# MAKEFILE
# Define variables

PROG = SHES

MODS =
#Place objects related to modules first
OBJS = \
       IO.o \
       main.o \
       EvolveData.o \
       Functions.o \
       GetInitialData.o \
       GetResults.o \
       OneDimensionalInterpolation.o \
       ReinitializeData.o \
       SpectralTransform.o \
       TricubicInterpolation.o

F90 = gfortran-mp-4.7
CC =

#F90 DEPENDENCIES
SWITCH = \
   -O3 -fopenmp
LIBS = \
   -I/usr/local/include -L/usr/local/lib -framework veclib \
   -fexternal-blas -fblas-matmul-limit=2 \
   -ISHTOOLS/modules \
   -LSHTOOLS/lib \
   -lSHTOOLS2.7 \
   -I/opt/local/include -L/opt/local/lib \
   -lfftw3 -lm -m64 \
   -L/usr/local/hdf5/lib -I/usr/local/hdf5/include \
   -lhdf5_fortran -lhdf5
FLAGS =

#C++ DEPENDENCIES
CCSWITCH =
CCLIBS =
CCFLAGS =

all: $(PROG)

$(PROG) : $(OBJS)
	$(F90) -o $(PROG) $(SWITCH) $(OBJS) $(LIBS)

%.mod: %.o $.f90
	$(F90) -c $(FLAGS) $(SWITCH) $(LIBS) $<
%.o: %.f90
	$(F90) -c $(FLAGS) $(SWITCH) $(LIBS) $<
%.o: %.cpp
	$(CC) -c $(CCFLAGS) $(CCSWITCH) $(CCLIBS) $<

# Cleaning everything
clean:
	rm -f $(PROG) $(OBJS) $(MODS) *.kmo       

# END MAKEFILE
