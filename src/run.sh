#PBS -N test
#PBS -l mppwidth=24
#PBS -j oe
#PBS -q short-b

#NPARALLEL=24
#NPROCESS=$NPARALLEL
NPARALLEL=24
NPROCESS=1
tempodir=snap
odir=pex.saitoh.c03.3
ifile=init/saitoh_init/pex.64.init.c03
fexe=run

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

mkdir $tempodir
cp $0 $tempodir

nptcl=`echo "$npernode / 24 * $NPROCESS * $OMP_NUM_THREADS" | bc`
aprun -n $NPROCESS -d $OMP_NUM_THREADS ./"$fexe" $ifile

mv $tempodir $odir
