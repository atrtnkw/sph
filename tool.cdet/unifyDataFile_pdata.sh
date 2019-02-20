if test $# -ne 3
then
    echo "sh $0 <bdir> <odir> <nfile>" >&2
    exit
fi

bdir=$1
odir=$2
nfile=$3

for n in $(seq -f "%04g" 1 1 $nfile)
do
    if test -e tmp.dat
    then
        rm -f tmp.dat
    fi
    tend=-1
    for ndir in $(seq -f "%02g" 0 1 99 )
    do
        idir="$bdir"/t"$ndir"
        if ! test -e $idir
        then
            continue
        fi
        ifile="$idir"/pdata."$n"
        echo $ifile >&2
        awk '{if($1>tend&&NR%16==0) print $0}' tend=$tend $ifile >> tmp.dat
        tend=`tail -n1 tmp.dat | awk '{print $1}'`
    done
    ofile="$odir"/id000000"$n".in
    mv tmp.dat $ofile
done

rm -f tmp.dat
