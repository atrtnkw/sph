#PBS -N test
#PBS -l mppwidth=24
#PBS -j oe
#PBS -q bulk-b
##PBS -q short-b
##PBS -q debug

NPARALLEL=24
NPROCESS=$NPARALLEL
odir=try03
#ifile=../result.heos/r001k/b1.10-1.00/s01.bsep/t01/t0400
ifile=../init.heos/r001k/rlxb_b1.10-1.00/b1.10-1.00.init
#
tempodir=snap
fexe=run

tag=`echo $ifile | awk '{print substr($0, length($0)-3, 4)}'`
if test $tag = init
then
    flag=0
else
    flag=1
fi

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

mkdir $tempodir
cp $0 $tempodir

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./"$fexe" $flag $ifile

mv $tempodir $odir
