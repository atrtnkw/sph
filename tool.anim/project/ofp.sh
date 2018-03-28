#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=00:15:00
##PJM -g xg17i022
#PJM -g xg17i056
#PJM -N project
#PJM -j
#PJM -L node=48
#PJM --mpi proc=3072
#PJM --omp thread=1

echo "../r016m/run.s1.00_h100-100/"  > input.list
echo "../r016m/run.s1.00_h100-100/anim10e09/" >> input.list
echo "86 128" >> input.list
echo "0 1 0 0" >> input.list
echo "0 0" >> input.list
echo "10e9" >> input.list

# idir
# odir
# tbgn tend
# normal vector and distance of the projection surface
# parallel displacement on the projection surface
# width of the projection surface

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
