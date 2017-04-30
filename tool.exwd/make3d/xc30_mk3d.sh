#PBS -N mk3d
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q short-b+

NPARALLEL=1
NPROCESS=$NPARALLEL
ifile=../r001m/rlx2_s1.00/s1.00_h060-060.data
ofile=../r001m/init/s1.00_h060-060.data
iflag=0
size=1e8
#size=5e7
#size=3e7
spotx=0.0
spoty=0.0
spotz=4.2e8

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

nptcl=`wc -l $ifile | awk '{print $1}'`
echo $nptcl  > input.list
echo $ifile >> input.list
echo $ofile >> input.list
echo $iflag >> input.list
echo $size  >> input.list
echo $spotx $spoty $spotz >> input.list

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list
