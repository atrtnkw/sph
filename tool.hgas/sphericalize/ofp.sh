#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=01:00:00
#PJM -g xg18i004
#PJM -N sphericalize
#PJM -j
#PJM -L node=3
#PJM --mpi proc=144
#PJM --omp thread=1

echo "../r004m/run.s1.00_h015-060/"  > input.list
echo "../r004m/run.s1.00_h015-060/init1d/" >> input.list
echo "128 128" >> input.list

# idir
# odir
# tbgn tend
# normal vector and distance of the projection surface
# parallel displacement on the projection surface
# width of the projection surface

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
