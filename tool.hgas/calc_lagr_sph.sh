if test $# -ne 3
then
    echo "sh $0 <idir> <tbgn> <tend>" >&2
    exit
fi

idir=$1
tbgn=$2
tend=$3

mlag="1 2 5 10 20 50 80 90 95 98 99"

for time in $(seq $tbgn 1 $tend)
do
    ptim=`printf "%04d" $time`
    sfile="$idir"/sph_t"$ptim".dat
    nptcl=`wc -l $sfile | awk '{print $1}'`
    awk 'BEGIN{m=cx=cy=cz=0.0;}{w=$3;cx+=w*$4;cy+=w*$5;cz+=w*$6;m+=w;}END{print cx/m, cy/m, cz/m}' \
        $sfile > tmp.center
    cx=`awk '{print $1;}' tmp.center`
    cy=`awk '{print $2;}' tmp.center`
    cz=`awk '{print $3;}' tmp.center`
    awk '{print ($4-cx)**2+($5-cy)**2+($6-cz)**2;}' cx=$cx cy=$cy cz=$cz $sfile | sort -g > tmp.r2
    printf "%4d" $time
    for ilag in $mlag
    do
        nline=`echo "$nptcl * $ilag * 0.01" | bc -l`
        awk 'BEGIN{flag=0;}{if(NR==int(nline)){printf(" %+e", sqrt($1));flag=1}}END{if(flag==0){printf(" %+e", 0.0);}}' nline=$nline tmp.r2
    done
    printf "\n"
done

rm -f tmp.center
rm -f tmp.r2
