if test $# -ne 3
then
    echo "sh $0 <tbgn> <tend> <dtsp>" >&2
    exit
fi

tbgn=$1
tend=$2
dtsp=$3
tdir=`dirname $0`
nfmesh=6

alias python3='~/work/yt-conda/bin/python3'
python3 "$tdir"/small2d.py $tbgn $tend $dtsp

for time in $(seq $tbgn $dtsp $tend)
do
    if test $time -eq $tend
    then
        continue
    fi
    time=`printf "%04d" $time`

    if ! test -e small_t"$time".dat
    then
        echo "small_t$time.dat is not found." >&2
        continue
    fi

    nf=`awk '{if(NR==1) print NF}' small_t"$time".dat`
    if test $nf = $nfmesh
    then
        echo "small_t$time.dat already exists." >&2
        continue
    fi

    awk '{printf("%+e %+e %+e %+e %+e %+e\n", $1, $2, 0., 0., $3, $4);}' small_t"$time".dat > tmp.dat
    n=`head -n1 small_t"$time".dat | awk '{print NF}'`
    if test $n -eq $nfmesh
    then
        rm -f tmp.dat
    else
        mv tmp.dat small_t"$time".dat
    fi
done
