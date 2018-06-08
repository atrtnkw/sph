#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=24:00:00
#PJM -g xg18i004
#PJM -N mk3d
#PJM -j
#PJM -L node=1
#PJM --mpi proc=1
#PJM --omp thread=1

NPARALLEL=1
NPROCESS=$NPARALLEL
ifile="../r004m/run.b1.00_h050-050_0.90_h056-100/nohot/s1.00_h050-050.data"
tfile="../r004m/run.b1.00_h050-050_0.90_h056-100/inhot/s1.00_h050-050"
iflag=0
size=1e8
#size=5e7
#size=3e7
spotx=0.0
spoty=0.0
#spotz=4.2e8 #s1.00_h100-100
spotz=3.6e8 #s1.00_h050-100
#spotz=3.6e8 #s1.00_h025-100

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

nptcl=`wc -l $ifile | awk '{print $1}'`
echo $nptcl       > input.list
echo $ifile      >> input.list
echo $tfile.data >> input.list
echo $iflag      >> input.list
echo $size       >> input.list
echo $spotx $spoty $spotz >> input.list
cp input.list $tfile.log

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
