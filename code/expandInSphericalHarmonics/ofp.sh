#!/bin/sh
#PJM -L rscgrp=regular-flat
#PJM -L elapse=02:00:00
#PJM -g xg18i004
#PJM -N harmonics
#PJM -j
#PJM -L node=2
#PJM --mpi proc=128
#PJM --omp thread=1

lfile=input.list1

echo "../../../wdmrg/"  > $lfile
echo "../r008k/run.b0.90-0.60/harmonics/" >> $lfile
echo "0 10000 1" >> $lfile

# idir
# odir
# tbgn tend dtsp

mpiexec.hydra -n ${PJM_MPI_PROC} ./run $lfile
