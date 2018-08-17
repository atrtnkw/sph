#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=00:30:00
#PJM -g xg18i004
#PJM -N project
#PJM -j
#PJM -L node=48
#PJM --mpi proc=3072
#PJM --omp thread=1

echo "../r004m/run.b1.00_h050-050_0.90/mod0001"  > input.list
echo "../r004m/run.b1.00_h050-050_0.90/anim02e11/" >> input.list
echo "42 50" >> input.list
echo "0 0 1 0" >> input.list
echo "0 0" >> input.list
echo "2e11" >> input.list

# idir
# odir
# tbgn tend
# normal vector and distance of the projection surface
# parallel displacement on the projection surface
# width of the projection surface

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
