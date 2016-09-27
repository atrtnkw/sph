if test $# -ne 3
then
    echo "sh $0 <basedir> <txx> <tyy>" >&2
    exit
fi

bdir=$1
dbgn=$2
dend=$3

prevt=-1

if test -e tmp.time
then
    rm -f tmp.time
fi

for ndir in $(seq -f "%02g" $dbgn 1 $dend)
do
    idir="$bdir"/t"$ndir"
    echo "start $idir" >&2
    ifile="$idir"/quad.log
    awk '{if($1>prevt) print $0;}' prevt=$prevt $ifile >> tmp.time
    prevt=`tail -n1 tmp.time | awk '{print $1}'`
done

cat tmp.time

rm -f tmp.time
