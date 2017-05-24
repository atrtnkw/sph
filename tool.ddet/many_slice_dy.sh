if test $# -ne 6
then
    echo "sh $0 <dy> <nprc> <idir> <tbgn> <tend> <dtsp>"
    exit
fi

dy=$1
nprc=`printf "%06d" $2`
idir=$3
tbgn=$4
tend=$5
dtsp=$6

for time in $(seq -f "%04g" $tbgn $dtsp $tend)
do
    for dnum in $(seq -f "%02g" 0 1 99)
    do
        if test -e sph_t"$time".dat.dy
        then
            continue
        fi
        if ! test -e "$idir"/t"$dnum"/sph_t"$time"_p"$nprc"_i000000.dat
        then
            echo "Not found "$idir"/t"$dnum"/sph_t"$time".dat"
            continue
        fi
        echo "Process $time ..." >&2
        awk -f ~/git-sph/tool.hgas/slice.awk dy="$dy" "$idir"/t"$dnum"/sph_t"$time"*.dat \
            > sph_t"$time".dat.dy
    done
done
