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
    bfile="$idir"/bhns_t"$ptim".dat
    if ! test -e $sfile
    then
        echo "File $sfile is not found"
        exit
    fi
    if ! test -e $bfile
    then
        echo "File $bfile is not found"
        exit
    fi
    bx=`awk '{print $4;}' $bfile`
    by=`awk '{print $5;}' $bfile`
    bz=`awk '{print $6;}' $bfile`
    awk 'BEGIN{cx=cy=cz=0.0;m=0.0;}
{w=$3;cx+=w*$4;cy+=w*$5;cz+=w*$6;m+=w;}
END{printf("%5d %+e\n", time, sqrt((cx/m-bx)**2+(cy/m-by)**2+(cz/m-bz)**2));}' \
    time=$time bx=$bx by=$by bz=$bz $sfile
done
