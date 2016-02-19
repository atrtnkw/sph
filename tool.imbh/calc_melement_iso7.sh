if test $# -ne 4
then
    echo "sh $0 <idir> <tbgn> <tend> <dtsp>" >&2
    exit
fi

idir=$1
tbgn=$2
tend=$3
dtsp=$4

echo "Directory $idir" >&2

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    awk -f ~/git-sph/tool.imbh/calc_melement_iso7.awk time=$time dtsp=$dtsp "$idir"/sph_t"$time".dat
done
