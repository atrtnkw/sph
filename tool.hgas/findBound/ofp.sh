#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=01:00:00
#PJM -g xg18i004
#PJM -N findBound
#PJM -j
#PJM -L node=48
#PJM --mpi proc=3072
#PJM --omp thread=1

echo "../r004m/run.b1.00_h030-060_0.60_h010-100/time010.00-050.00/"  > input.list
echo "./hoge" >> input.list
echo "10 10" >> input.list
echo "1" >> input.list
echo "0" >> input.list
echo "0" >> input.list

# idir
# odir
# tbgn tend
# searchmode(0:bound, 1:unbound)
# printmode(0:all, 1:summary only)
# excludedID(-1/0/1)

odir=`awk '{if(NR==2) print $0}' input.list`
if ! test -e $odir
then
    mkdir -p $odir
fi

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
