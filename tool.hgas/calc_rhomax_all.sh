if test $# -ne 3
then
    echo "sh $0 <idir> <tbgn> <tend>" 1>&2
    exit
fi

idir=$1
tbgn=$2
tend=$3

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    file="$idir"/sph_t"$time".dat
    if test -e $file
    then
        ttemp=`expr $time + 0`
        printf "%d " $ttemp
        awk '{print $16;}' $file | sort -g | tail -n1
    fi
done
