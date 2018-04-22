if test $# -ne 4
then
    echo "sh $0 <idir> <tbgn> <tend> <dtsp>" >&2
    exit
fi

idir=$1
tbgn=$2
tend=$3
dtsp=$4

python ~/git-sph/tool.cdet/plot1d.py $idir $tbgn $tend $dtsp

for time in $(seq -f "%04g" $tbgn $dtsp $tend)
do
    if test $time -eq $tend
    then
        continue
    fi

    awk '{for(i=0;i<20;i++) printf(" %+e", $(2*i+1));printf("\n");}' mesh_t"$time".dat > tmp.dat
    n=`head -n1 mesh_t"$time".dat | awk '{print NF}'`
    if test $n -eq 20
    then
        rm -f tmp.dat
    else
        mv tmp.dat mesh_t"$time".dat
    fi
done
