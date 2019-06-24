#!/bin/sh
##PJM -L rscgrp=regular-flat
#PJM -L rscgrp=debug-flat
#PJM -L elapse=00:30:00
#PJM -g xg18i004
#PJM -N extractCoreDet
#PJM -j
#PJM -L node=48
#PJM --mpi proc=3072
#PJM --omp thread=1

odir=../r004m/run.b1.00_h050-050_0.90_h060-100/coreDetonation

echo "../r004m/run.b1.00_h050-050_0.90_h060-100/time000.00-005.00"  > input.list
echo "$odir"       >> input.list
echo "64 80 64"  >> input.list
echo "1"           >> input.list
echo "5e8 30."     >> input.list
echo "1.0 0.9"     >> input.list
echo "1.62234e+09" >> input.list
echo "3.94296e+08" >> input.list

# idir
# odir
# tbgn tend nsnap
# targetWD
# Radius Phi
# M0 M1 [Msun]
# BinaryRadius   [cm]
# BinaryVelocity [cms^-1}

if ! test -e $odir
then
    mkdir $odir
fi
cp input.list $odir
mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
