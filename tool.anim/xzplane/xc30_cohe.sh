#PBS -N xzplane
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q bulk-a

cd $PBS_O_WORKDIR

echo "/work/tanikwat/git-sph/nswd/data/r512k/run.bns-wd1.0/unfy"  > input.list
echo "/work/tanikwat/git-sph/nswd/data/r512k/run.bns-wd1.0/anim/xz03e09" >> input.list
echo "1 132" >> input.list
echo "273 311" >> input.list
echo "-3e9 0. -3e9" >> input.list
echo "6e9 256" >> input.list

aprun -n 1 -d 1 ./run input.list

# idir
# otype
# tbgn tend
# xmin[0] xmin[1] xmin[2]
# width nx
