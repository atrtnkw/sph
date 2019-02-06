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
        file1="$idir"/t"$ptim"*000000.hexa
        if test -e $file1
        then
            echo "delete $idir/t$ptim.hexa"
            rm -f "$idir"/t"$ptim"*.hexa
        fi
    done
done
