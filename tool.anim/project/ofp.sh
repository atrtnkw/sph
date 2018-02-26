#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=24:00:00
##PJM -g xg17i022
#PJM -g xg17i056
#PJM -N project
#PJM -j
#PJM -L node=12
#PJM --mpi proc=768
#PJM --omp thread=1

echo "/work/1/xg17i022/x10049/git-sph/imbh/data/r016k/run.onwd1.30_bh3e1_b05.00_rt4/"  > input.list
echo "/work/1/xg17i022/x10049/git-sph/imbh/data/project" >> input.list
echo "304 336" >> input.list
echo "0 0 1 0" >> input.list
echo "0 0" >> input.list
echo "2e9" >> input.list

# idir
# odir
# tbgn tend
# normal vector and distance of the projection surface
# parallel displacement on the projection surface
# width of the projection surface

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
