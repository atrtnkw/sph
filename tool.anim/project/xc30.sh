#PBS -N project
#PBS -l mppwidth=72
#PBS -j oe
##PBS -q bulk-b+
#PBS -q short-b+

NPROCESS=72
OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

echo "/work/tanikwat/git-sph/ddet/data/r512k/run.s1.00_h100-100/t00"  > input.list
echo "/work/tanikwat/git-sph/ddet/data/project" >> input.list
echo "0 4" >> input.list
echo "0 1 0 0" >> input.list
echo "1.2e9" >> input.list

# idir
# odir
# tbgn tend
# normal vector (3), constant (1)
# width

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list
