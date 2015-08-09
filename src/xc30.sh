#PBS -N test
#PBS -l mppwidth=192
#PBS -j oe
#PBS -q short-b
##PBS -q debug

NPARALLEL=192
NPROCESS=$NPARALLEL

odir=r001k
ifile=../init.heos/r001k/b1.10-1.00.init
#odir=pex
#ifile=../init.ideal/glass/pex_ns064.init
tempodir=snap
fexe=run

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

mkdir $tempodir
cp $0 $tempodir

nptcl=`echo "$npernode / 24 * $NPROCESS * $OMP_NUM_THREADS" | bc`
aprun -n $NPROCESS -d $OMP_NUM_THREADS ./"$fexe" $ifile

mv $tempodir $odir
