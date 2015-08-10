#PBS -N test
#PBS -l mppwidth=144
#PBS -j oe
#PBS -q bulk-b
##PBS -q short-b
##PBS -q debug

NPARALLEL=144
NPROCESS=$NPARALLEL
odir=test
#ifile=../../init.heos/r064k/nons/s1.00_t000.init
ifile=../init.heos/r064k/nons/s1.00_t000.init
#ifile=../init.heos/r008k/nons/s1.10_t000.init
#
tempodir=snap
fexe=run

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

mkdir $tempodir
cp $0 $tempodir

nptcl=`echo "$npernode / 24 * $NPROCESS * $OMP_NUM_THREADS" | bc`
aprun -n $NPROCESS -d $OMP_NUM_THREADS ./"$fexe" $ifile

mv $tempodir $odir
