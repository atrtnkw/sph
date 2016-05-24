#PBS -N mk3d
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q short-a

NPARALLEL=1
NPROCESS=$NPARALLEL
#ifile=../../../exwd/t002/sph_t0153.dat
#ofile=../../../exwd/init/b1.1-0.9_he0100-050_t0153.data
#iflag=0
#size=1e8
#spotx=-6.735527e+08
#spoty=+8.732118e+07
#spotz=-1.899386e+07
ifile=../../../exwd.2/hedet/t001/sph_t0117.dat
ofile=../../../exwd.2/hedet/b1.1-0.9_he0100-050_2nd.data
iflag=1
size=5e8
spotx=-3.5e8
spoty=-8e8
spotz=-1.899386e+07

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
