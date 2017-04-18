if test $# -ne 4
then
    echo "sh $0 <idir> <tbgn> <tend> <dt>" 1>&2
    exit
fi

tdir=$1
tbgn=$2
tend=$3
dsnp=$4

for idir in $tdir
do
    echo "start directory $idir"
    for time in $(seq $tbgn 1 $tend)
    do
        res=`echo "$time % $dsnp" | bc`
        if test $res -eq 0
        then
            continue
        fi
        ptim=`printf "%04d" $time`
        file0="$idir"/sph_t"$ptim".dat
        file1="$idir"/sph_t"$ptim"*000000.dat
        if test -e $file0
        then
            echo "NOT delete $idir/sph_t$ptim.dat"
        elif test -e $file1
        then
            echo "delete $idir/sph_t$ptim.dat"
            rm -f "$idir"/sph_t"$ptim"*.dat
        fi
    done
done
