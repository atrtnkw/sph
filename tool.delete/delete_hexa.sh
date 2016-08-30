if test $# -ne 2
then
    echo "sh $0 <idir> <dt>" 1>&2
    exit
fi

tdir=$1
dsnp=$2

for idir in $tdir
do
    echo "start directory $idir"
    for time in $(seq -f "%04g" 0 $dsnp 9999)
    do
        file="$idir"/t"$time"_p000000.hexa
        if test -e $file
        then
            tbgn=$time
            break
        fi
    done
    for time in $(seq -f "%04g" $tbgn $dsnp 9999)
    do
        file="$idir"/t"$time"_p000000.hexa
        if test -e $file
        then
            tend=$time
        fi
        if ! test -e $file
        then
            break
        fi
    done
    echo "tbgn: $tbgn tend: $tend" 1>&2
    for time in $(seq -f "%04g" $tbgn $dsnp $tend)
    do
        if test $time -ne $tbgn -a  $time -ne $tend
        then
            rm -f "$idir"/t"$time"_p*.hexa
        fi
    done
    touch "$idir"/deleted.log
done
