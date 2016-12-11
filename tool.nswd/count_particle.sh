if test $# -ne 3
then
    echo "sh $0 <idir> <tbgn> <tend>" >&2
    exit
fi

idir=$1
tbgn=$2
tend=$3

echo "Directory $idir" >&2

for time in $(seq $tbgn 1 $tend)
do
    ptim=`printf "%04d" $time`
#    wc -l "$idir"/sph_t"$ptim".dat | awk '{print time, $1}' time=$time
    cat "$idir"/sph_t"$ptim"*.dat | wc -l | awk '{print time, $1}' time=$time
done
