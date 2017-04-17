#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=24:00:00
#PJM -g xg17i022
#PJM -N bounce
#PJM -j
#PJM -L node=1
#PJM --mpi proc=64
#PJM --omp thread=1

itype=../r256k/run.cowd1.0_bh3e2_b04.00/t00/sph_t0015
fflag=1
nfile=768
xmin0=-2e9
ymin0=-2e9
xmax=1e11
width=4e9
nnxxx=1024

echo "$itype"  > input.list
echo "$fflag $nfile" >> input.list
echo "$xmin0 $ymin0 $xmax" >> input.list
echo "$width $nnxxx" >> input.list

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list

