if test $# -ne 1
then
    echo "sh $0 <idir>" 1>&2
    exit
fi

idir=$1

for time in $(seq -f "%04g" 0 1 9999)
do
    file="$idir"/t"$time"_p000000.hexa
    if test -e $file
    then
        tbgn=$time
        break
    fi
done

for time in $(seq -f "%04g" $tbgn 1 9999)
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

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    if test $time -ne $tbgn -a  $time -ne $tend
    then
        rm -f "$idir"/t"$time"_p*.hexa
    fi
done
