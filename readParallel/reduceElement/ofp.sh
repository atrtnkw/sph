#!/bin/sh
#PJM -L rscgrp=regular-flat
##PJM -L rscgrp=debug-flat
#PJM -L elapse=00:30:00
#PJM -g xg18i004
#PJM -N alignParticleRecord
#PJM -j
#PJM -L node=48
#PJM --mpi proc=3072
##PJM -L node=1
##PJM --mpi proc=64
#PJM --omp thread=1

idir=../../../ddet2
#odir=../r004m/run.b1.00_h030-060_0.60_h010-100/rad/init
odir=./
ibgn=0
iend=40
tsnp=0.0625

####
echo "$idir"  > input.list
echo "$odir" >> input.list
echo "$ibgn" >> input.list
echo "$iend" >> input.list
echo "$tsnp" >> input.list

if ! test -e $odir
then
    mkdir -p $odir
fi

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
