#PBS -N project
#PBS -l mppwidth=144
#PBS -j oe
##PBS -q bulk-b+
#PBS -q short-b+

NPROCESS=144
OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

echo "/work/tanikwat/git-sph/ddet/data/r008m/run.s1.00_h100-100/unfy"  > input.list
echo "/work/tanikwat/git-sph/ddet/data/r008m/run.s1.00_h100-100/anim_12e08" >> input.list
echo "81 128" >> input.list
echo "0 1 0 0" >> input.list
echo "12e8" >> input.list

# idir
# odir
# tbgn tend
# normal vector (3), constant (1)
# width

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list
