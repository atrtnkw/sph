#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=24:00:00
#PJM -g xg17i022
#PJM -N bounce
#PJM -j
#PJM -L node=1
#PJM --mpi proc=64
#PJM --omp thread=1

## The number of processes must be 2^n.
## The number of "nnxxx" must be 2^n.

itbgn=14
itend=16
iform=../r001m/run.cowd1.0_bh3e2_b04.50/t00/sph_t
odir=../r001m/run.cowd1.0_bh3e2_b04.50/fitting
#
fflag=1
nfile=768
xmin0=-2e9
ymin0=-2e9
xmax=1e11
width=4e9
nnxxx=512

if ! test -e $odir
then
    mkdir $odir
fi

for itime in $(seq -f "%04g" $itbgn 1 $itend)
do
    itype="$iform""$itime"
    otype="$odir"/t"$itime"
    
    echo "$itype"               > input.list
    echo "$fflag $nfile"       >> input.list
    echo "$xmin0 $ymin0 $xmax" >> input.list
    echo "$width $nnxxx"       >> input.list
    echo "$otype"              >> input.list
    
    mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list

done

