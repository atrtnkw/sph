#PBS -N test
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q short-b

NPARALLEL=1
NPROCESS=$NPARALLEL

odir=shock_1d
ifile=../init.ideal/shock_1d.init
tempodir=snap
fexe=run

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

mkdir $tempodir
cp $0 $tempodir

nptcl=`echo "$npernode / 24 * $NPROCESS * $OMP_NUM_THREADS" | bc`
aprun -n $NPROCESS -d $OMP_NUM_THREADS ./"$fexe" $ifile

mv $tempodir $odir
