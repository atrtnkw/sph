if test $# -ne 3
then
    echo "sh $0 <tbgn> <tend> <dtsp>" >&2
    exit
fi

tbgn=$1
tend=$2
dtsp=$3

alias python3='~/work/yt-conda/bin/python3'
python3 ~/git-sph/tool.imbh/plot1d.py $tbgn $tend $dtsp

for time in $(seq $tbgn $dtsp $tend)
do
    if test $time -eq $tend
    then
        continue
    fi
    time=`printf "%04d" $time`

    if ! test -e mesh_t"$time".dat
    then
        echo "mesh_t$time.dat is not found." >&2
        continue
    fi

    awk '{for(i=0;i<17;i++) printf(" %+e", $(2*i+1));printf("\n");}' mesh_t"$time".dat > tmp.dat
    n=`head -n1 mesh_t"$time".dat | awk '{print NF}'`
    if test $n -eq 20
    then
        rm -f tmp.dat
    else
        mv tmp.dat mesh_t"$time".dat
    fi
done
