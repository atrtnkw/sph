#PBS -N disk
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q bulk-a

cd $PBS_O_WORKDIR

echo "/work/tanikwat/git-sph/nswd/data/r256k/run.bns-wd1.0/unfy"  > input.list
echo "/work/tanikwat/git-sph/nswd/data/r256k/run.bns-wd1.0/data/disk" >> input.list
#echo "/work/tanikwat/git-sph/nswd/data/disk/disk" >> input.list
echo "640 670" >> input.list

aprun -n 1 -d 1 ./run input.list

# idir
# otype
# tbgn tend
