#PBS -N hotspot
#PBS -l mppwidth=24
#PBS -j oe
#PBS -q bulk-a

NPARALLEL=24
NTHREAD=1

cd $PBS_O_WORKDIR

echo "/work/tanikwat/git-sph/nswd/data/r256k/run.bns-wd1.0/unfy"  > input.list
echo "/work/tanikwat/git-sph/nswd/data/r256k/run.bns-wd1.0/data/hotspot" >> input.list
echo "659 660" >> input.list

NPROCESS=`echo "$NPARALLEL / $NTHREAD" | bc`
export OMP_NUM_THREADS=$NTHREAD
aprun -n $NPROCESS -d $NTHREAD ./run input.list

#export OMP_NUM_THREADS=24
#aprun -n 1 -d 24 ./run input.list

# idir
# otype
# tbgn tend
