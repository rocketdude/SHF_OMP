# MAKEFILE
# Define variables

PROG = SHF

OBJS = IO.o main.o CalculateSchwarzschildMetric.o CalculateInitialData.o EvaluatedSdr.o EvaluateF.o EvaluateG.o EvolveData.o FindU.o Functions.o InitializeFilter.o InitializeMatrices.o ReinitializeData.o ReadHDFMetric.o SphHarmonicY.o Spline.o

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

main.o: IO.o CalculateSchwarzschildMetric.o CalculateInitialData.o EvaluatedSdr.o EvaluateF.o EvaluateG.o EvolveData.o FindU.o Functions.o InitializeFilter.o InitializeMatrices.o ReadHDFMetric.o ReinitializeData.o SphHarmonicY.o Spline.o

%.o: %.f90
	$(F90) -c $(FLAGS) $(SWITCH) $<
%.o: %.cpp
	$(CC) -c $(CCFLAGS) $(CCSWITCH) $(CCLIBS) $<

# Cleaning everything
clean:
	rm -f $(PROG) $(OBJS) *.kmo

# END MAKEFILE