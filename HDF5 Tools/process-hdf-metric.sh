#!/bin/bash

##This script is used to process the HDF5 metric.
##We separate the data into a single file containing 1 metric component, 1 rl, and 1 iteration.
##This shell script uses:
##1. hdf5 1.8 library for the h5copy utility
##2. hdf5_slicer.cc compiled as 'slicer' (mpicc -o slicer hdf5_slicer.cc -<Include HDF5LIBS>)
##3. Modified hdf5_merge.c compiled as 'merge' (mpicc -o merge hdf5_merge.c -<Include HDF5LIBS>)

##On TACC--Ranger or Lonestar--use these after loading the hdf5 library:
#module swap pgi intel
#module load hdf5/1.8.8
#mpicc -o slicer hdf5_slicer.cc -I$TACC_HDF5_INC -Wl,-rpath,$TACC_HDF5_LIB -L$TACC_HDF5_LIB -lhdf5 -lz
#mpicc -O3 -o merge hdf5_merge.c -I$TACC_HDF5_INC -Wl,-rpath,$TACC_HDF5_LIB -L$TACC_HDF5_LIB -lhdf5 -lz

##User should define these variables
RL=9 #Relevant level
CMAX=12 #Total number of processors outputting the files
ITSTART=0 #The lowest iteration number you want
ITMAX=400 #Maximum iteration number you want
LSNAME=KRANC2BSSNCHIMATTER #Dataset name for lapse and shift functions
SPNAME=ADMBASE #Dataset name for the spatial metric

##Initialize the variables
IT=$ITSTART
C=0

echo "Maximum number of processors is $CMAX"
echo "Processing up to iteration $ITMAX"
##
##Process the alpha and beta[xyz], both have dataset names $LSNAME
##
for METRIC in alpha beta1 beta2 beta3; do
  while [ "$IT" -le "$ITMAX" ]; do
  if [ $((IT % 4)) -eq 0 ]; then
    echo "Processing $METRIC at it = $IT"
    while [ "$C" -lt "$CMAX" ]; do         
      #Extract dataset from the file and save it into a new file
      h5copy -i "$METRIC.h5" -o "$METRIC.c=$C.it=$IT.rl=$RL.h5" \
      -s "$LSNAME::$METRIC it=$IT tl=0 rl=$RL c=$C" \
      -d "$LSNAME::$METRIC it=$IT tl=0 rl=$RL c=$C"
      #Increase the value of C
      let "C+=1"
    done
  #Combine the extracted dataset into a single file using the hdf5_slicer
  ./slicer --match "$LSNAME::$METRIC it=$IT tl=0 rl=$RL" \
  --out3d-cube $METRIC.c=*.it=$IT.rl=$RL.h5 $METRIC.it=$IT.rl=$RL.h5
  #Merge together the number of processor components into a single group
  #based on iteration number
  #./merge -g -u "$METRIC.it=$IT.rl=$RL.h5" "$METRIC.it=$IT.rl=$RL.merged.h5"
  #Reinitialize the value of C
  let "C=0"
  fi
  #Increase the iteration number
  let "IT+=1"
  done
  #Reinitialize the value of IT
  let "IT=$ITSTART"
done

##
##Process the spatial metric g[xyz][xyz], all have dataset names $SPNAME
##
for METRIC in gxx gyy gzz gxy gxz gyz; do
  while [ "$IT" -le "$ITMAX" ]; do
  if [ $((IT % 4)) -eq 0 ]; then  
    echo "Processing $METRIC at it = $IT"
      while [ "$C" -lt "$CMAX" ]; do         
      #Extract dataset from the file and save it into a new file
      h5copy -i "$METRIC.h5" -o "$METRIC.c=$C.it=$IT.rl=$RL.h5" \
      -s "$SPNAME::$METRIC it=$IT tl=0 rl=$RL c=$C" \
      -d "$SPNAME::$METRIC it=$IT tl=0 rl=$RL c=$C"
      #Increase the value of C
      let "C+=1"
    done
  #Combine the extracted dataset into a single file using the hdf5_slicer
  ./slicer --match "$SPNAME::$METRIC it=$IT tl=0 rl=$RL" \
  --out3d-cube $METRIC.c=*.it=$IT.rl=$RL.h5 $METRIC.it=$IT.rl=$RL.h5
  #Merge together the number of processor components into a single group
  #based on iteration number
  #./merge -g -u "$METRIC.it=$IT.rl=$RL.h5" "$METRIC.it=$IT.rl=$RL.merged.h5"
  #Reinitialize the value of C
  let "C=0"
  fi
  #Increase the iteration number
  let "IT+=1"
  done
  #Reinitialize the value of IT
  let "IT=$ITSTART"
done
#Remove the chunked data
rm *.c=*
exit 0
