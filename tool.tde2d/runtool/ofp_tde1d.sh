#!/bin/sh
##PJM -L rscgrp=regular-cache
##PJM -L elapse=24:00:00
#PJM -L rscgrp=debug-cache
#PJM -L elapse=00:30:00
#PJM -g xg18i004
#PJM -N fimbh2
#PJM -j
#PJM -L node=1
#PJM --mpi proc=64
#PJM --omp thread=1

idir=data/r064m/run.hewd0.45_bh3e2_b07.00/initTde2d
ifile=column109.dat
#
fexe=flash4

module load phdf5/1.8.17

if test -e "$idir"/"$ifile"
then
    cp "$idir"/"$ifile" init.dat
    mpiexec.hydra -n ${PJM_MPI_PROC} ./"$fexe"
else
    echo "$idir/$ifile is not found"
fi
