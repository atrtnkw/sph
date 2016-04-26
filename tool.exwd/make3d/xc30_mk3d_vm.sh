#PBS -N mk3d
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q short-a

NPARALLEL=1
NPROCESS=$NPARALLEL
ifile=../r064k/run.b1.1-1.0.mrg/unfy/sph_t0144.dat
ofile=../r064k/init/b1.1-1.0.ig.t0144.data
iflag=1
hsize=5e8
spotx=+4.818394e+07
spoty=-1.103376e+08
spotz=+1.723009e+07

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

nptcl=`wc -l $ifile | awk '{print $1}'`
echo $nptcl  > input.list
echo $ifile >> input.list
echo $ofile >> input.list
echo $iflag >> input.list

echo 

echo $hsize >> input.list
echo $spotx >> input.list
echo $spoty >> input.list
echo $spotz >> input.list

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list
