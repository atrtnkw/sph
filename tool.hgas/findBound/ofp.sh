#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=01:00:00
#PJM -g xg18i004
#PJM -N findBound
#PJM -j
#PJM -L node=12
#PJM --mpi proc=768
#PJM --omp thread=1

echo "../../mod0001"  > input.list
echo "./hoge" >> input.list
echo "197 200" >> input.list
echo "1" >> input.list

# idir
# odir
# tbgn tend
# printmode(0:all, 1:summary only)

odir=`awk '{if(NR==2) print $0}' input.list`
if ! test -e $odir
then
    mkdir -p $odir
fi

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
