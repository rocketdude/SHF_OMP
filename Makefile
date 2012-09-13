# MAKEFILE
# Define variables

PROG = SHF

MODS = DynMetricArray.mod
#Place objects related to modules first
OBJS = DynMetricArray.o \
       IO.o \
       main.o \
       EvaluatedSdr.o \
       EvaluateJacobian.o \
       EvaluateMatrixofMetric.o \
       EvolveData.o \
       FindU.o \
       Functions.o \
       GetAllMetricComponents.f90 \
       GetFilter.o \
       GetInitialData.o \
       GetMatrices.o \
       GetMetricAtCurrentTime.o \
       GetMetricComponent.o \
       GetSchwarzschildMetric.o \
       Invert3Metric.o \
       OneDimensionalInterpolation.o \
       ReadHDF5MetricData.o \
       ReinitializeData.o \
       SphHarmonicY.o \
       TemporalMetricInterpolation.o \
       TricubicInterpolation.o

F90 = mpif90
CC = mpicc

#F90 DEPENDENCIES
SWITCH = -O3 -xW -mcmodel=large -openmp
LIBS = -Wl,-rpath,$$TACC_MKL_LIB \
   -L$$TACC_MKL_LIB -lmkl -lguide \
   -I$$TACC_HDF5_INC -I$$TACC_HDF5_LIB -L$$TACC_HDF5_LIB -lhdf5_fortran -lhdf5 -lz
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
	$(F90) -c $(FLAGS) $(SWITCH) $(LIBS) $<
%.o: %.f90
	$(F90) -c $(FLAGS) $(SWITCH) $(LIBS) $<
%.o: %.cpp
	$(CC) -c $(CCFLAGS) $(CCSWITCH) $(CCLIBS) $<

# Cleaning everything
clean:
	rm -f $(PROG) $(OBJS) $(MODS) *.kmo       

# END MAKEFILE
