#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=01:00:00
#PJM -g xg18i004
#PJM -N sphericalize
#PJM -j
#PJM -L node=48
#PJM --mpi proc=3072
#PJM --omp thread=1

#echo "../r004m/run.b1.00_h030-060_0.60/findunbound"  > input.list
#echo "./"    >> input.list
#echo "1e13 1e8"  >> input.list
#echo "50 50" >> input.list

echo "../r004m/run.b1.00_h030-060_0.60/findbound"  > input.list
echo "./"    >> input.list
echo "1e10 1e5"  >> input.list
echo "50 50" >> input.list

# idir
# odir
# tbgn tend

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
