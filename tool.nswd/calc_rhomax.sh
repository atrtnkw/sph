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
    awk 'BEGIN{rhomax=0.;}{if($16>rhomax){rhomax=$16;line=$0;}}END{print line;}' \
        "$idir"/sph_t"$time".dat
done
