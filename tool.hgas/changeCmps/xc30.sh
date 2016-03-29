#PBS -N mk3d
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q short-a

NPARALLEL=1
NPROCESS=$NPARALLEL
ifile=data/s1.15_rlx2.dat
ofile=../init3d/s1.15_heco002.data
iflag=0
hefrc=0.02

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

nptcl=`wc -l $ifile | awk '{print $1}'`
echo $nptcl  > input.list
echo $ifile >> input.list
echo $ofile >> input.list
echo $iflag >> input.list

nhelium=`echo "$nptcl * $hefrc" | bc -l | awk '{printf("%d\n", $1);}'`
awk '{print $4**2+$5**2+$6**2;}' $ifile | sort -g \
    | tail -n"$nhelium" | head -n1 | awk '{print sqrt($1);}' >> input.list

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list
