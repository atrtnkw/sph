if test $# -ne 1
then
    echo "sh $0 <idir>"
    exit
fi

idir=$1
tbgn=0
tend=500

echo "Directory $idir" >&2

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    wc -l "$idir"/sph_t"$time".dat | awk '{print time, $1}' time=$time
done
