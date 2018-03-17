#PBS -N mk3d
#PBS -l mppwidth=1
#PBS -j oe
##PBS -q short-b+
#PBS -q bulk-b+

NPARALLEL=1
NPROCESS=$NPARALLEL
ifile=../r008m/run.s1.00_h003-020/nohot/s1.00_h003-020.data
tfile=../r008m/run.s1.00_h003-020/init/s1.00_h003-020
iflag=0
size=1e8
#size=5e7
#size=3e7
spotx=0.0
spoty=0.0
#spotz=4.2e8 #s1.00_h100-100
#spotz=3.6e8 #s1.00_h050-100
spotz=3.6e8 #s1.00_h025-100

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

nptcl=`wc -l $ifile | awk '{print $1}'`
echo $nptcl       > input.list
echo $ifile      >> input.list
echo $tfile.data >> input.list
echo $iflag      >> input.list
echo $size       >> input.list
echo $spotx $spoty $spotz >> input.list
cp input.list $tfile.log

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list
