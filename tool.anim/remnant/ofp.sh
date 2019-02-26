#!/bin/sh
#PJM -L rscgrp=regular-flat
#PJM -L elapse=02:00:00
#PJM -g xg18i004
#PJM -N remnant1
#PJM -j
#PJM -L node=2
#PJM --mpi proc=128
#PJM --omp thread=1

lfile=input.list1

echo "../../../wdmrg/t00000-01000"  > $lfile
echo "./small.entr" >> $lfile
echo "0 1000 1" >> $lfile

echo "0 0 1 0" >> $lfile # for xyplane
echo "0 0" >> $lfile
echo "6e9" >> $lfile

# idir
# odir
# tbgn tend
# normal vector and distance of the projection surface
# parallel displacement on the projection surface
# width of the projection surface

mpiexec.hydra -n ${PJM_MPI_PROC} ./run $lfile
