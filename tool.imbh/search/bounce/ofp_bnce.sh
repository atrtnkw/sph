#!/bin/sh
##PJM -L rscgrp=regular-cache
#PJM -L rscgrp=regular-flat
#PJM -L elapse=24:00:00
#PJM -g xg17i022
#PJM -N bounce
#PJM -j
#PJM -L node=16
#PJM --mpi proc=128
#PJM --omp thread=8

## The number of processes must be 2^n.
## The number of "nnxxx" must be 2^n.

itbgn=122
itend=132
idtsp=1
idir=../r004m/run.cowd1.05_bh3e2_b04.50_nw4
odir=../r004m/run.cowd1.05_bh3e2_b04.50_nw4/fitting.mach02.+-
#
fflag=1
nfile=1536
xmin0=-3e9
ymin0=-3e9
xmax=1e11
width=6e9
nnxxx=512
#

if ! test -e $odir
then
    mkdir $odir
fi

cp $0 $odir/ofp_bnce.sh

for itime in $(seq -f "%04g" $itbgn $idtsp $itend)
do
    for dnum in $(seq -f "%02g" 0 1 99)
    do
        pnfile=`printf "%06d" $nfile`
        if ! test -e "$idir"/t"$dnum"/sph_t"$itime"_p"$pnfile"_i000000.dat
        then
            echo "Not found"$idir"/t"$dnum"/sph_t"$itime"_p"$pnfile".dat"
            continue
        fi
        
        itype="$idir"/t"$dnum"/sph_t"$itime"
        otype="$odir"/t"$itime"
        
        echo "$itype"               > input.list
        echo "$fflag $nfile"       >> input.list
        echo "$xmin0 $ymin0 $xmax" >> input.list
        echo "$width $nnxxx"       >> input.list
        echo "$otype"              >> input.list
        
        mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
        
        break
    done
done

