#!/bin/sh
#PJM -L rscgrp=regular-flat
#PJM -L elapse=00:30:00
#PJM -g xg17i022
#PJM -N xyplane
#PJM -j
#PJM -L node=1
#PJM --mpi proc=4
#PJM --omp thread=16
#PJM -o stdouterr.log
#PJM -e stdouterr.log

echo "../r001m/run.hewd0.45_bh3e2_b07.00/t00"  > input.list
echo "hoge/hoge" >> input.list
echo "1 768" >> input.list
echo "100 100" >> input.list
echo "-4e9 -4e9 0." >> input.list
echo "8e9 256" >> input.list

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list

# idir
# otype
# file_flag
# tbgn tend
# xmin[0] xmin[1] xmin[2]
# width nx
