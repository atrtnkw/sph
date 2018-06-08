#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=01:00:00
#PJM -g xg18i004
#PJM -N project
#PJM -j
#PJM -L node=48
#PJM --mpi proc=3072
#PJM --omp thread=1

#echo "../../../ddet2/mod0064"  > input.list
#echo "../r004m/run.b1.00_h030-060_0.60/anim02e09/mod0064" >> input.list
#echo "120 128" >> input.list
#echo "0 0 1 0" >> input.list
#echo "-6e8 -6e8" >> input.list
#echo "2e9" >> input.list

echo "../../../ddet2/mod0001"  > input.list
#echo "../r004m/run.b1.00_h030-060_0.60/anim20e10/mod0001" >> input.list
echo "../r004m/run.b1.00_h030-060_0.60/anim20e10.yz/mod0001" >> input.list
echo "50 50" >> input.list
#echo "0 0 1 0" >> input.list
echo "1 0 0 0" >> input.list
echo "0 0" >> input.list
echo "2e11" >> input.list

# idir
# odir
# tbgn tend
# normal vector and distance of the projection surface
# parallel displacement on the projection surface
# width of the projection surface

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
