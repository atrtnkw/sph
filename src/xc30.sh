#PBS -N rb0081110
#PBS -l mppwidth=144
#PBS -j oe
##PBS -q large-b
#PBS -q bulk-b
##PBS -q short-b
##PBS -q debug

NPARALLEL=144
NPROCESS=$NPARALLEL
odir=rlxb_b1.10-1.00
ifile=r008k/nonb/b1.10-1.00.init
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
