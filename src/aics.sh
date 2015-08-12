if test $# -ne 2
then
    echo "$0 <nproc> <ifile>"
    exit
fi

export OMP_NUM_THREADS=1

nproc=$1
ifile=$2

tag=`echo $2 | awk '{print substr($0, length($0)-3, 4)}'`
if test $tag = init
then
    flag=0
else
    flag=1
fi

odir=snap
if ! test -e $odir
then
    mkdir $odir
fi

mpirun -np $nproc ./run $flag $ifile
