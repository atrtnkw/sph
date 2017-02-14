#PBS -N xzplane
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q bulk-a

cd $PBS_O_WORKDIR

rdir="/work/tanikwat/git-sph/cohe/data/r008k/run.bco1.00-he0.45"
echo "$rdir/unfy"  > input.list
echo "$rdir/anim/xz05e09" >> input.list
echo "0 1" >> input.list
echo "0 800" >> input.list
echo "-5e9 0. -5e9" >> input.list
echo "10e9 256" >> input.list

aprun -n 1 -d 1 ./run input.list

# idir
# otype
# tbgn tend
# xmin[0] xmin[1] xmin[2]
# width nx
