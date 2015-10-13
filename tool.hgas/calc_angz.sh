if test $# -ne 3
then
    echo "sh $0 <idir> <tbgn> <tend>"
    exit
fi

idir=$1
tbgn=$2
tend=$3

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    file="$idir"/sph_t"$time".dat
    awk -f calc_angz.awk id=0 time=$time $file
done
