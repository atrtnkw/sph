#!/bin/sh
#PJM -L rscgrp=debug-flat
#PJM -L elapse=00:05:00
#PJM -g xg18i004
#PJM -N extractParticle
#PJM -j
#PJM -L node=48
#PJM --mpi proc=3072
##PJM -L node=1
##PJM --mpi proc=64
#PJM --omp thread=1

idir=../r004m/run.b1.00_h030-060_0.60_h010-100/rad/time000.00-005.00
odir=hoge
time=0
id1st=0
idint=320
mode=1

echo "$idir"          > input.list
echo "$odir"         >> input.list
echo "$time"         >> input.list
echo "$id1st $idint" >> input.list
echo "$mode"         >> input.list

if ! test -e $odir
then
    mkdir -p $odir
fi

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list

ptime=`printf "%04d" $time`
cat "$odir"/ext*.dat | sort -n > "$odir"/ext_t"$ptime".dat
