#PBS -N project
#PBS -l mppwidth=144
#PBS -j oe
##PBS -q bulk-b+
#PBS -q short-b+

NPROCESS=144
OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

echo "/work/tanikwat/git-sph/ddet/data/r001m/run.s1.00_h090-090/unfy"  > input.list
echo "/work/tanikwat/git-sph/ddet/data/r001m/run.s1.00_h090-090/anim02e09" >> input.list
echo "80 128" >> input.list
echo "0 1 0 0" >> input.list
echo "0 0" >> input.list
echo "2e9" >> input.list

# idir
# odir
# tbgn tend
# normal vector and distance of the projection surface
# parallel displacement on the projection surface
# width of the projection surface

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list
