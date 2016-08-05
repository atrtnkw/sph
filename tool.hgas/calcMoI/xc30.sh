#PBS -N calcMoI
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q bulk-a

cd $PBS_O_WORKDIR

echo "/work/tanikwat/git-sph/nswd/data/r016k/run.std/unfy"  > input.list
echo "/work/tanikwat/git-sph/nswd/data/r016k/run.std/data" >> input.list
echo "2e9" >> input.list
echo "0 500" >> input.list

aprun -n 1 -d 1 ./run input.list
