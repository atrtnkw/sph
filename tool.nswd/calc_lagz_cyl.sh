if test $# -ne 4
then
    echo "sh $0 <idir> <nptcl> <tbgn> <tend>"
    exit
fi

idir=$1
nptcl=$2
tbgn=$3
tend=$4
mlag="10 20 30 40 50 60 70 80 90"

for time in $(seq $tbgn 1 $tend)
do
    ptim=`printf "%04d" $time`
    sfile="$idir"/sph_t"$ptim".dat
    dnptcl=`wc -l $sfile | awk '{print nptcl - $1}' nptcl=$nptcl`
    bfile="$idir"/bhns_t"$ptim".dat
    bx=`awk '{print $4;}' $bfile`
    by=`awk '{print $5;}' $bfile`
    bz=`awk '{print $6;}' $bfile`
    awk '{print ($6-bz)**2;}' bx=$bx by=$by bz=$bz $sfile | sort -g > tmp.r2
    printf "%4d" $time
    for ilag in $mlag
    do
        nline=`echo "$nptcl * $ilag * 0.01 - $dnptcl" | bc -l`
        awk 'BEGIN{flag=0;}{if(NR==int(nline)){printf(" %+e", sqrt($1));flag=1}}END{if(flag==0){printf(" %+e", 0.0);}}' nline=$nline tmp.r2
    done
    printf "\n"
done

rm -f tmp.*
