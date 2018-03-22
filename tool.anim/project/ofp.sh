#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=24:00:00
##PJM -g xg17i022
#PJM -g xg17i056
#PJM -N project
#PJM -j
#PJM -L node=24
#PJM --mpi proc=1536
#PJM --omp thread=1

#echo "/work/1/xg17i022/x10049/git-sph/imbh/data/r064m/run.hewd0.45_bh3e2_b06.00/"  > input.list
echo "/work/1/xg17i022/x10049/git-sph/imbh/"  > input.list
#echo "/work/1/xg17i022/x10049/git-sph/imbh/data/r064m/run.hewd0.45_bh3e2_b06.00/anim10e09/" >> input.list
echo "/work/1/xg17i022/x10049/git-sph/imbh/anim10e09/" >> input.list
echo "92 100" >> input.list
echo "0 0 1 0" >> input.list
echo "0 0" >> input.list
echo "10e9" >> input.list

# idir
# odir
# tbgn tend
# normal vector and distance of the projection surface
# parallel displacement on the projection surface
# width of the projection surface

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
