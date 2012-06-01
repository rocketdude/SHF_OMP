# MAKEFILE
# Define variables

PROG = SHF

MODS = DynMetricArray.mod
#Place objects related to modules first
OBJS = DynMetricArray.o IO.o main.o EvaluatedSdr.o EvaluateJacobian.o EvaluateMatrixofMetric.o EvolveData.o FindU.o Functions.o GetFilter.o GetInitialData.o GetMatrices.o GetMetric.o GetMetricComponent.o GetSchwarzschildMetric.o Invert3Metric.o progress.o ReadHDF5MetricData.o ReinitializeData.o SphHarmonicY.o Spline.o TricubicInterpolation.o

F90 = h5fc
CC = h5c++

#F90 DEPENDENCIES
SWITCH = -I/usr/local/include -L/usr/local/lib -lhdf5 -framework veclib -fexternal-blas -fblas-matmul-limit=2 -O2 -fopenmp 
LIBS = 
FLAGS =

#C++ DEPENDENCIES
CCSWITCH =
CCLIBS = 
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
