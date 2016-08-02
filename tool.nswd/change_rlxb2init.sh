if test $# -ne 2
then
    echo "sh $0 <idir> <ofile>" >&2
    exit
fi

idir=$1
ofile=$2

awk -f ~/git-sph/tool.hgas/change_t2t.awk "$idir"/final.dat > "$ofile".data
awk -f ~/git-sph/tool.hgas/change_t2t_bhns.awk "$idir"/bhns_final.dat > "$ofile".bhns
