if test $# -ne 2
then
    echo "$0 <nproc> <ifile>"
    exit
fi

nproc=$1
ifile=$2

odir=snap

if ! test -e $odir
then
    mkdir $odir
fi

export OMP_NUM_THREADS=1
mpirun -np $nproc ./run $ifile
