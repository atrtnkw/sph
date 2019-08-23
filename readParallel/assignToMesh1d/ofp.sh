#!/bin/sh
##PJM -L rscgrp=regular-flat
#PJM -L rscgrp=debug-flat
#PJM -L elapse=00:30:00
#PJM -g xg18i004
#PJM -N assignToMesh
#PJM -j
#PJM -L node=96
#PJM --mpi proc=3072
#PJM --omp thread=1

echo "../r004m/run.b1.00_h030-060_0.60_h010-100/time010.00-050.00"  > input.list
#echo "../r004m/run.b1.00_h030-060_0.60_h010-100/snr" >> input.list
echo "./" >> input.list
echo "20 20" >> input.list

# idir
# odir
# tbgn tend

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
