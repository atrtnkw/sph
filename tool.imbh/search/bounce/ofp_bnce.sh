#!/bin/sh
##PJM -L rscgrp=regular-cache
#PJM -L rscgrp=regular-flat
##PJM -L elapse=24:00:00
#PJM -L elapse=00:30:00
#PJM -g xg17i022
#PJM -N bounce
#PJM -j
##PJM -L node=32
#PJM -L node=64
#PJM --mpi proc=128
#PJM --omp thread=16
#PJM -o stdouterr.log
#PJM -e stdouterr.log

## The number of processes must be 2^n.
## The number of "nnxxx" must be 2^n.

itbgn=103
itend=103
idtsp=1
#idir=../r032m/run.hewd0.45_bh3e2_b07.00
idir=../../../imbh.3
odir=../r064m/run.hewd0.45_bh3e2_b07.00/fitting.chkmach04_r008m
nfile=3072
#nfile=1536
#nfile=768
#
fflag=1
xmin0=-5e9
ymin0=-5e9
xmax=1e11
width=10e9
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

mv stdouterr.log $odir
