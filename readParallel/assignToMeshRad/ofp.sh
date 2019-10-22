#!/bin/sh
##PJM -L rscgrp=regular-flat
#PJM -L rscgrp=debug-flat
#PJM -L elapse=00:10:00
#PJM -g xg18i004
#PJM -N assignToMeshRad
#PJM -j
#PJM -L node=96
#PJM --mpi proc=3072
#PJM --omp thread=1

echo "../r004m/run.b1.00_h030-060_0.60_h010-100/rad/time005.00-050.00"  > input.list
#echo "../r004m/run.b1.00_h030-060_0.60_h010-100/rad/anly/heelem_t0010" >> input.list
echo "coelem_t0010" >> input.list
echo "../r004m/run.b1.00_h030-060_0.60_h010-100/rad/element/sphelem_t0010" >> input.list
echo "10" >> input.list
echo "1"  >> input.list

# idir
# otype
# etype
# itime
# flag(0:CompanionHeCO,1:CompanionCO)

mpiexec.hydra -n ${PJM_MPI_PROC} ./run input.list
