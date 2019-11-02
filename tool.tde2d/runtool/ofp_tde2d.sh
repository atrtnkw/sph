#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=24:00:00
##PJM -L rscgrp=debug-cache
##PJM -L elapse=00:30:00
#PJM -g xg18i004
#PJM -N tde2d
#PJM -j
#PJM -L node=1
#PJM --mpi proc=64
#PJM --omp thread=1

idir=data/r064m/run.hewd0.45_bh3e2_b07.00/initTde2d
ifile=initTde2d.dat

###############################

if test -e "$idir"/"$ifile"
then
    module load phdf5/1.8.17
    cp $idir/$ifile init.dat
    mpiexec.hydra -n ${PJM_MPI_PROC} ./flash4
else
    echo "Error: File $idir/$ifile is not found"
fi
