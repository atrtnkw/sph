#PBS -N test_heos
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q short-a

NPARALLEL=1
NPROCESS=$NPARALLEL
fexe=run
tin=1e9

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./"$fexe" $tin
