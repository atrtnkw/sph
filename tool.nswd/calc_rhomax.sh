if test $# -ne 3
then
    echo "sh $0 <idir> <tbgn> <tend>"
    exit
fi

idir=$1
tbgn=$2
tend=$3

echo "Directory $idir" >&2

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    awk 'BEGIN{rhomax=0.;}{if($16>rhomax){rhomax=$16;line=$0;}}END{print line;}' \
        "$idir"/sph_t"$time".dat
done
