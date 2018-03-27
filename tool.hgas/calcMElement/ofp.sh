#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=01:00:00
##PJM -g xg17i022
#PJM -g xg17i056
#PJM -N project
#PJM -j
#PJM -L node=3
#PJM --mpi proc=144
#PJM --omp thread=1

echo "../r001m/run.s1.00_h100-100/"  > input.list
echo "../r001m/run.s1.00_h100-100/anly/" >> input.list
echo "16 16" >> input.list

# idir
# odir
# tbgn tend
# normal vector and distance of the projection surface
# parallel displacement on the projection surface
# width of the projection surface

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
