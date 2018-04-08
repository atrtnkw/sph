#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=01:00:00
#PJM -g xg18i004
#PJM -N mktd
#PJM -j
#PJM -L node=96
#PJM --mpi proc=6144
#PJM --omp thread=1

itype=../r064m/run.hewd0.45_bh3e2_b06.00/t05/sph_t0100
otype=../r064m/run.hewd0.45_bh3e2_b06.00/init.nuc/hewd0.45_bh3e2_b06.00_t0100
nvecw="1.6 1.0 0.0 1e9"

echo $itype >  input.list
echo $otype >> input.list
echo $nvecw >> input.list

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
