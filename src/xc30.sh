#PBS -N test
#PBS -l mppwidth=144
#PBS -j oe
#PBS -q bulk-b
##PBS -q short-b
##PBS -q debug

NPARALLEL=144
NPROCESS=$NPARALLEL
odir=test
ifile=../init.heos/r008k/rlxb_b1.10-1.00/b1.10-1.00.init
#ifile=b1.10-1.00/test.a/t0090
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
