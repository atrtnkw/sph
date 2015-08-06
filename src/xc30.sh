#PBS -N test
#PBS -l mppwidth=24
#PBS -j oe
#PBS -q short-b

NPARALLEL=24
NPROCESS=$NPARALLEL

odir=pex
ifile=../init.ideal/pex.init
tempodir=snap
fexe=run

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

mkdir $tempodir
cp $0 $tempodir

nptcl=`echo "$npernode / 24 * $NPROCESS * $OMP_NUM_THREADS" | bc`
aprun -n $NPROCESS -d $OMP_NUM_THREADS ./"$fexe" $ifile

mv $tempodir $odir
