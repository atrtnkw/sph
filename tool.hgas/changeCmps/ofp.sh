#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=20:00:00
#PJM -g xg17i022
#PJM -N chcmp
#PJM -j
#PJM -L node=1
#PJM --mpi proc=1
#PJM --omp thread=1

NPARALLEL=1
NPROCESS=$NPARALLEL
ifile=../r001m/rlx2_s0.60/final.dat
#ofile=../r001m/nons/s0.60_h100-100.data
ofile=hoge.data
rmin=6.67303218e8
fhe=1.0
fc=0.0
fo=0.0
fne=0.0
fmg=0.0

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

nptcl=`wc -l $ifile | awk '{print $1}'`
echo $nptcl  > input.list
echo $ifile >> input.list
echo $ofile >> input.list
echo $rmin  >> input.list
echo $fhe   >> input.list
echo $fc    >> input.list
echo $fo    >> input.list
echo $fne   >> input.list
echo $fmg   >> input.list

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
