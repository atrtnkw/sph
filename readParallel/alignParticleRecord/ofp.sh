#!/bin/sh
##PJM -L rscgrp=regular-flat
#PJM -L rscgrp=debug-flat
#PJM -L elapse=00:05:00
#PJM -g xg18i004
#PJM -N alignParticleRecord
#PJM -j
#PJM -L node=48
#PJM --mpi proc=3072
##PJM -L node=1
##PJM --mpi proc=64
#PJM --omp thread=1

idir=../../../ddet2
odir=hoge
id1st=0
idint=32

####
if ! test -e $idir/id.list
then
    echo "Error: Not found $idir/id.list" >&2
    exit
else
    echo "Found $idir/id.list" >&2
fi
nptcl=`head -n1 $idir/id.list | awk '{print $1}'`

echo "$idir"          > input.list
echo "$odir"         >> input.list
echo "$nptcl"        >> input.list
echo "$id1st $idint" >> input.list

if ! test -e $odir
then
    mkdir -p $odir
fi

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
