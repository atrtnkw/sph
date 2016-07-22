if test $# -ne 2
then
    echo "sh $0 <idir> <time>" >&2
    exit
fi

idir=$1
time=$2

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
END{print sqrt((cx/m-bx)**2+(cy/m-by)**2+(cz/m-bz)**2);}' bx=$bx by=$by bz=$bz $sfile
