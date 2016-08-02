if test $# -ne 3
then
    echo "sh $0 <idir> <tbgn> <tend>" >&2
    exit
fi

idir=$1
tbgn=$2
tend=$3

echo "Directory $idir" >&2

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    wc -l "$idir"/sph_t"$time".dat | awk '{print time, $1}' time=$time
done
