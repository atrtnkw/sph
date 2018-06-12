#PBS -N project
#PBS -l mppwidth=144
#PBS -j oe
##PBS -q bulk-b+
#PBS -q short-b+

NPROCESS=144
OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

echo "/work/tanikwat/git-sph/ddet/data/r004m/run.s1.00_h015-060/"  > input.list
echo "./" >> input.list
echo "128 128" >> input.list

# idir
# odir
# tbgn tend

if ! test -e $adir
then
    mkdir -p $adir
fi

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list
