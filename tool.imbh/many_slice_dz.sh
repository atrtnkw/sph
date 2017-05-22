if test $# -ne 5
then
    echo "sh $0 <dz> <nprc> <tbgn> <tend> <dtsp>"
    exit
fi

dz=$1
nprc=`printf "%06d" $2`
tbgn=$3
tend=$4
dtsp=$5

for time in $(seq -f "%04g" $tbgn $dtsp $tend)
do
    for dnum in $(seq -f "%02g" 0 1 99)
    do
        if test -e sph_t"$time".dat.dz
        then
            continue
        fi
        if ! test -e t"$dnum"/sph_t"$time"_p"$nprc"_i000000.dat
        then
            echo "Not found t"$dnum"/sph_t"$time".dat"
            continue
        fi
        echo "Process $time ..." >&2
        awk -f ~/git-sph/tool.hgas/slice.awk dz="$dz" t"$dnum"/sph_t"$time"*.dat \
            > sph_t"$time".dat.dz
    done
done
