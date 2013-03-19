# MAKEFILE
# Define variables

PROG = SHF

MODS = dynmetricarray.mod
#Place objects related to modules first
OBJS = \
       DynMetricArray.o \
       IO.o \
       main.o \
       CalculateArea.o \
       EvaluateJacobian.o \
       EvaluateMatrixofMetric.o \
       EvolveData.o \
       FindU.o \
       Functions.o \
       GetAllMetricComponents.o \
       GetInitialData.o \
       GetMetricAtCurrentTime.o \
       GetMetricComponent.o \
       GetKerrSchildMetric.o \
       Invert3Metric.o \
       OneDimensionalInterpolation.o \
       ReadHDF5MetricData.o \
       ReinitializeData.o \
       SpectralTransform.o \
       TricubicInterpolation.o

F90 = mpif90
CC = mpicc

#F90 DEPENDENCIES
SWITCH = -O3 -mkl -xhost -mcmodel=large -openmp
INC = \
   -I$$MKLROOT/include \
   -I$$TACC_HDF5_INC \
   -I$$TACC_FFTW3_INC 
LIBS = \
   -free -ISHTOOLS/modules \
   -LSHTOOLS/lib \
   -lSHTOOLS2.7 \
   -mkl=parallel \
   -Wl,-rpath,$$TACC_HDF5_LIB -L$$TACC_HDF5_LIB \
   -lhdf5_fortran -lhdf5 -lz \
   -L$$TACC_FFTW3_LIB -lfftw3 -lm
FLAGS = -fpp

#C++ DEPENDENCIES
CCSWITCH =
CCINC =
CCLIBS =
CCFLAGS =

all: $(PROG)

$(PROG) : $(OBJS)
	$(F90) -o $(PROG) $(SWITCH) $(OBJS) $(LIBS)

%.mod: %.o $.f90
	$(F90) -c $(FLAGS) $(SWITCH) $(INC) $< $(LIBS)
%.o: %.f90
	$(F90) -c $(FLAGS) $(SWITCH) $(INC) $< $(LIBS)
%.o: %.cpp
	$(CC) -c $(CCFLAGS) $(CCSWITCH) $(CCINC) $< $(CCLIBS)

# Cleaning everything
clean:
	rm -f $(PROG) $(OBJS) $(MODS) *.kmo       

# END MAKEFILE
