# MAKEFILE
# Define variables

PROG = SHF

MODS = DynMetricArray.mod
#Place objects related to modules first
OBJS = DynMetricArray.o IO.o main.o EvaluatedSdr.o EvaluateF.o EvaluateG.o EvaluateJacobian.o EvolveData.o FindU.o Functions.o GetFilter.o GetInitialData.o GetMatrices.o GetMetric.o GetSchwarzschildMetric.o Invert3Metric.o progress.o ReadHDF5Metric.o ReadNInterpolateMetric.o ReinitializeData.o SphHarmonicY.o Spline.o TricubicInterpolation.o

F90 = mpif90
CC = mpicc

#F90 DEPENDENCIES
SWITCH = -O3 -xW -mcmodel=large -openmp
LIBS = -Wl,-rpath,$$TACC_MKL_LIB \
   -L$$TACC_MKL_LIB -lmkl -lguide \
   -I$$TACC_HDF5_INC -I$$TACC_HDF5_LIB -L$$TACC_HDF5_LIB -lhdf5_fortran -lhdf5 -lz -cxxlib
FLAGS = -fpp

#C++ DEPENDENCIES
CCSWITCH = -O3 -xW -mcmodel=large -openmp
CCLIBS = -I$$TACC_MKL_INC -Wl,-rpath,$$TACC_MKL_LIB \
   -L$$TACC_MKL_LIB -lmkl -lguide \
   -I$$TACC_HDF5_INC -Wl,-rpath,$$TACC_HDF5_LIB -L$$TACC_HDF5_LIB -lhdf5 -lz
CCFLAGS =

all: $(PROG)

$(PROG) : $(OBJS)
	$(F90) -o $(PROG) $(SWITCH) $(OBJS) $(LIBS)

%.mod: %.o $.f90
	$(F90) -c $(FLAGS) $(SWITCH) $<
%.o: %.f90
	$(F90) -c $(FLAGS) $(SWITCH) $<
%.o: %.cpp
	$(CC) -c $(CCFLAGS) $(CCSWITCH) $(CCLIBS) $<

# Cleaning everything
clean:
	rm -f $(PROG) $(OBJS) $(MODS) *.kmo       

# END MAKEFILE
