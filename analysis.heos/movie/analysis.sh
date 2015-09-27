if test $# -ne 5
then
    echo "sh $0 <nthread> <m> <n/1k> <ifile> <ofile>"
    exit
fi

nthrd=$1
tmass=$2
nreso=$3
ifile=$4
ofile=$5
nproc=1

nptcl=`echo "$tmass * 10 * $nreso * 1024" | bc`

export OMP_NUM_THREADS=$1

echo "Nproc: $nproc"
echo "Nthrd: $OMP_NUM_THREADS"
echo "Nptcl: $nptcl"

mpirun -np $nproc ./run $nptcl $ifile $ofile
