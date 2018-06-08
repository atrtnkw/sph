if test $# -ne 1
then
    echo "sh $0 <basedir>" >&2
    exit
fi

bdir=$1

dbgn=0
dend=99
prevt=-1

if test -e tmp.time
then
    rm -f tmp.time
fi

for ndir in $(seq -f "%02g" $dbgn 1 $dend)
do
    idir="$bdir"/t"$ndir"
#    echo "start $idir" >&2
    ifile="$idir"/time.log
    if ! test -e $ifile
    then
        continue
    fi
    awk '{if($1=="#"||$3>prevt) print $0;}' prevt=$prevt $ifile >> tmp.time
    prevt=`tail -n1 tmp.time | awk '{print $3}'`
done

cat tmp.time

rm -f tmp.time
