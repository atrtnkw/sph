if test $# -ne 3
then
    echo "sh $0 <idir> <tbgn> <tend>" >&2
    exit
fi

idir=$1
tbgn=$2
tend=$3

for time in $(seq $tbgn 1 $tend)
do
    ptim=`printf "%04d" $time`
    sfile="$idir"/sph_t"$ptim".dat
    awk 'BEGIN{m=cx=cy=cz=0.0;}{w=$3;cx+=w*$4;cy+=w*$5;cz+=w*$6;m+=w;}END{print cx/m, cy/m, cz/m}' \
        $sfile > tmp.center
    cx=`awk '{print $1;}' tmp.center`
    cy=`awk '{print $2;}' tmp.center`
    cz=`awk '{print $3;}' tmp.center`
    printf "%4d" $time
    awk 'BEGIN{s2=0.;n=0;}{s2+=($4-cx)**2+($5-cy)**2+($6-cz)**2;n++;}END{printf(" %+e", sqrt(s2/n))}' cx=$cx cy=$cy cz=$cz $sfile
    printf "\n"
done

rm -f tmp.center
rm -f tmp.r2
