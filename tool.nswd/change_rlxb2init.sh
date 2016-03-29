if test $# -ne 1
then
    echo "sh $0 <idir>"
    exit
fi

idir=$1

awk -f ~/git-sph/tool.hgas/change_t2t.awk "$idir"/final.dat \
    | awk -f ~/git-sph/tool.hgas/zero_alphu.awk \
    > "$idir"/bns-wd0.6.atc0.data
awk -f ~/git-sph/tool.hgas/change_t2t_bhns.awk "$idir"/bhns_final.dat \
    > "$idir"/bns-wd0.6.atc0.bhns
