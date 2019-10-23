#!/bin/sh
##PJM -L rscgrp=regular-flat
#PJM -L rscgrp=debug-flat
#PJM -L elapse=00:30:00
#PJM -g xg18i004
#PJM -N extractHeDet
#PJM -j
#PJM -L node=48
#PJM --mpi proc=3072
#PJM --omp thread=1

odir=../r004m/run.b1.00_h030-060_0.60/surfaceDetonation

echo "../r004m/run.b1.00_h030-060_0.60/time000.00-001.00"  > input.list
echo "$odir"  >> input.list
echo "0 16"     >> input.list
echo "0"        >> input.list
echo "4.6e8"    >> input.list
echo "90. 270." >> input.list

# idir
# odir
# tbgn tend
# targetWD
# rcore
# phi0 phi1

if ! test -e $odir
then
    mkdir $odir
fi
cp input.list $odir
mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
