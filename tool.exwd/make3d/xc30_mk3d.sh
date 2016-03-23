#PBS -N mk3d
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q short-a

NPARALLEL=1
NPROCESS=$NPARALLEL
#ifile=../r016k/run.b1.1-1.0/unfy/sph_t0035.dat
#ofile=../r016k/init/t0035_b1.1-1.0.data
ifile=../../t0035_b1.1-1.0_ph01/sph_t0064.dat
ofile=../../t0035_b1.1-1.0_ph01/b1.1-1.0_phase01.data
iflag=1
hefrc=0.1

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
