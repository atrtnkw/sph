#PBS -N project
#PBS -l mppwidth=144
#PBS -j oe
##PBS -q bulk-b+
#PBS -q short-b+

NPROCESS=144
OMP_NUM_THREADS=1
#adir="/work/tanikwat/git-sph/ddet/data/r004m/run.s1.00_h030-060/anim02e09"
adir="/work/tanikwat/git-sph/ddet/data/r004m/run.s1.00_h015-060/anim02e09"

cd $PBS_O_WORKDIR

#echo "/work/tanikwat/git-sph/ddet/data/r004m/run.s1.00_h030-060/"  > input.list
echo "/work/tanikwat/git-sph/ddet/data/r004m/run.s1.00_h015-060/"  > input.list
echo $adir >> input.list
echo "18 109" >> input.list
echo "0 1 0 0" >> input.list
echo "0 0" >> input.list
echo "2e09" >> input.list

# idir
# odir
# tbgn tend
# normal vector and distance of the projection surface
# parallel displacement on the projection surface
# width of the projection surface

if ! test -e $adir
then
    mkdir -p $adir
fi

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list
