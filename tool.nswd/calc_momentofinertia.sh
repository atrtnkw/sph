if test $# -ne 3
then
    echo "sh $0 <idir> <tbgn> <tend>" >&2
    exit
fi

idir=$1
tbgn=$2
tend=$3

tdir="$HOME"/git-sph

for time in $(seq $tbgn 1 $tend)
do
    ptim=`printf "%04d" $time`
    sfile="$idir"/sph_t"$ptim".dat
    bfile="$idir"/bhns_t"$ptim".dat

    cat $sfile $bfile \
        | awk 'BEGIN{cx=cy=cz=0.;;m=0.;}{cx+=$3*$4;cy+=$3*$5;cz+=$3*$6;m+=$3;}END{print cx/m, cy/m, cz/m;}' \
        > tmp.cm
    cx=`awk '{print $1;}' tmp.cm`
    cy=`awk '{print $2;}' tmp.cm`
    cz=`awk '{print $3;}' tmp.cm`
    awk -f "$tdir"/tool.nswd/calc_momentofinertia.awk time=$time cx=$cx cy=$cy cz=$cz $sfile $bfile
done

rm -f tmp.*
