#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=24:00:00
#PJM -g xg18i004
#PJM -N project
#PJM -j
#PJM -L node=12
#PJM --mpi proc=768
#PJM --omp thread=1

echo "../../tempdir"  > input.list
#echo "../r001m/run.b1.00_h100-100_0.60/anim10e09/mod0001" >> input.list
echo "../r001m/run.b1.00_h100-100_0.60/anim20e09/mod0001" >> input.list
#echo "16 100" >> input.list
echo "18 100" >> input.list
echo "0 0 1 0" >> input.list
echo "0 0" >> input.list
#echo "10e9" >> input.list
echo "20e9" >> input.list

# idir
# odir
# tbgn tend
# normal vector and distance of the projection surface
# parallel displacement on the projection surface
# width of the projection surface

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
