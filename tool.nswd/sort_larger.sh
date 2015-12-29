if test $# -ne 3
then
    echo "sh $0 <idir> <time(????)> <ni>"
    exit
fi

idir=$1
time=$2
nint=$3

bfile="$idir"/bhns_t"$time".dat
sfile="$idir"/sph_t"$time".dat

xbh=`awk '{print $4}' $bfile`
ybh=`awk '{print $5}' $bfile`
zbh=`awk '{print $6}' $bfile`
nnow=`wc -l $sfile | awk '{print $1}'`

awk '{print sqrt(($4-xbh)**2+($5-ybh)**2)}' xbh=$xbh ybh=$ybh zbh=$zbh $sfile \
    | sort -g \
    | awk '{print $1, (NR+nint-nnow)/nint}' nint=$nint nnow=$nnow
