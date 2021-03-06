if test $# -ne 3
then
    echo "sh $0 <tbgn> <tend> <dtsp>" >&2
    exit
fi

tbgn=$1
tend=$2
dtsp=$3
tdir=`dirname $0`
nfmesh=19

alias python3='~/work/yt-conda/bin/python3'
#python3 ~/git-sph/tool.imbh/plot1d.py $tbgn $tend $dtsp
python3 "$tdir"/plot1d.py $tbgn $tend $dtsp

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

    nf=`awk '{if(NR==1) print NF}' mesh_t"$time".dat`
    if test $nf = $nfmesh
    then
        echo "mesh_t$time.dat already exists." >&2
        continue
    fi

    awk '{for(i=0;i<nf-2;i++){if(i==0||i==1){printf(" %+e", 0.);}printf(" %+e", $(2*i+1));}printf("\n");}' nf=$nfmesh mesh_t"$time".dat \
        > tmp.dat
    n=`head -n1 mesh_t"$time".dat | awk '{print NF}'`
    if test $n -eq 20
    then
        rm -f tmp.dat
    else
        mv tmp.dat mesh_t"$time".dat
    fi
done
