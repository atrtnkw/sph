#!/bin/sh
#PJM -L rscgrp=regular-flat
#PJM -L elapse=00:30:00
#PJM -g jh190021
#PJM -N reduceElement
#PJM -j
#PJM -L node=48
#PJM --mpi proc=3072
##PJM -L node=1
##PJM --mpi proc=64
#PJM --omp thread=1

idir=../r004m/run.b0.90_0.60_h010-100/time000.00-003.00
odir=./
ibgn=48
iend=48
tsnp=0.0625
flag=-1

####
echo "$idir"  > input.list
echo "$odir" >> input.list
echo "$ibgn" >> input.list
echo "$iend" >> input.list
echo "$tsnp" >> input.list
echo "$flag" >> input.list

if ! test -e $odir
then
    mkdir -p $odir
fi

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
