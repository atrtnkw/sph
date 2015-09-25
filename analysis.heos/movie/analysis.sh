if test $# -ne 3
then
    echo "sh $0 <nproc> <ifile> <ofile>"
    exit
fi

nproc=$1
ifile=$2
ofile=$3

echo "Input file is being made..."
wc -l $ifile | awk '{print $1}' > tmp.dat
cat $ifile >> tmp.dat

echo "The program starts..."
mpirun -np $nproc ./run tmp.dat $ofile
