#!/bin/sh
#PJM -L rscgrp=regular-cache
#PJM -L elapse=48:00:00
#PJM -g xg18i004
#PJM -N subCh
#PJM -j
#PJM -L node=1
#PJM --mpi proc=64
#PJM --omp thread=1

module load phdf5/1.8.17

mpiexec.hydra -n ${PJM_MPI_PROC} ./flash4
