#!/bin/sh
##PJM -L rscgrp=regular-cache
#PJM -L rscgrp=debug-cache
#PJM -L elapse=00:30:00
#PJM -g xg18i004
#PJM -N remnant1
#PJM -j
#PJM -L node=2
#PJM --mpi proc=128
#PJM --omp thread=1

echo "../../../wdmrg2"  > input.list
echo "../r008k/run.b0.90-0.60_norot/xyplane/" >> input.list
echo "0 500 100" >> input.list
echo "0 0 1 0" >> input.list # for xyplane
#echo "0 1 0 0" >> input.list # for xzplane
echo "0 0" >> input.list
echo "4e9" >> input.list

# idir
# odir
# tbgn tend
# normal vector and distance of the projection surface
# parallel displacement on the projection surface
# width of the projection surface

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
