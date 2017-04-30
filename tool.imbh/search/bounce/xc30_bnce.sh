#PBS -N bounce
#PBS -l mppwidth=16
#PBS -j oe
##PBS -q large-b+
##PBS -q bulk-b+
#PBS -q short-b+
##PBS -q debug

## The number of processes must be 2^n.
## The number of "nnxxx" must be 2^n.

NPROCESS=16
export OMP_NUM_THREADS=1

itype=../r256k/run.cowd1.0_bh1e2_b04.00/t00/sph_t0056
fflag=1
nfile=72
xmin0=-2e9
ymin0=-2e9
xmax=1e11
width=4e9
#nnxxx=1024
nnxxx=512

cd $PBS_O_WORKDIR

echo "$itype"  > input.list
echo "$fflag $nfile" >> input.list
echo "$xmin0 $ymin0 $xmax" >> input.list
echo "$width $nnxxx" >> input.list

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list
