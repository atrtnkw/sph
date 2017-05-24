#PBS -N mk3d
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q short-b+

NPARALLEL=1
NPROCESS=$NPARALLEL
ifile=../r002m/rlx2_s1.05/s1.05_h030-060.data
tfile=../r002m/run.s1.05_h030-060/init/s1.05_h030-060
iflag=0
size=1e8
#size=5e7
#size=3e7
spotx=0.0
spoty=0.0
#spotz=4.2e8 #s1.00_h100-100
#spotz=3.9e8 #s1.00_h150-100
#spotz=3.7e8 #s1.00_h200-100
spotz=3.6e8  #s1.05_h050-100

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
