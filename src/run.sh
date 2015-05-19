#PBS -N test
#PBS -l mppwidth=48
#PBS -j oe
#PBS -q short-b

NPARALLEL=48
NPROCESS=48
tempodir=snap
odir=kh
ifile=init/kh.init
fexe=run

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

mkdir $tempodir
cp $0 $tempodir

nptcl=`echo "$npernode / 24 * $NPROCESS * $OMP_NUM_THREADS" | bc`
aprun -n $NPROCESS -d $OMP_NUM_THREADS ./"$fexe" $ifile

mv $tempodir $odir
