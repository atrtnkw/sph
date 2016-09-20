#PBS -N xyplane
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q bulk-a

cd $PBS_O_WORKDIR

echo "/work/tanikwat/git-sph/nswd/data/r128k/run.bns-wd1.0/unfy"  > input.list
echo "/work/tanikwat/git-sph/nswd/data/r128k/run.bns-wd1.0/anim/xy03e09" >> input.list
echo "1031 1500" >> input.list
echo "-3e9 -3e9 0." >> input.list
echo "6e9 256" >> input.list
#echo "-5e9 -5e9 0." >> input.list
#echo "10e9 256" >> input.list

aprun -n 1 -d 1 ./run input.list

# idir
# otype
# tbgn tend
# xmin[0] xmin[1] xmin[2]
# width nx
