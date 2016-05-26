#PBS -N mk3d
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q short-a

NPARALLEL=1
NPROCESS=$NPARALLEL
ifile=../../../imbh.2/rlx2_s1.20/final.dat
ofile=../../../imbh.2/rlx2_s1.20/s1.20_onemg_t100.data
iflag=0
rmin=0.0
fhe=0.0
fc=0.0
fo=0.6
fne=0.35
fmg=0.05

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

nptcl=`wc -l $ifile | awk '{print $1}'`
echo $nptcl  > input.list
echo $ifile >> input.list
echo $ofile >> input.list
echo $rmin  >> input.list
echo $fhe   >> input.list
echo $fc    >> input.list
echo $fo    >> input.list
echo $fne   >> input.list
echo $fmg   >> input.list

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list
