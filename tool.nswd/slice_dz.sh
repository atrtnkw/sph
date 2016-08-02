if test $# -ne 3
then
    echo "sh $0 <idir> <time:????> <dz>" >&2
    exit
fi

idir=$1
time=$2
delz=$3

bhns="$idir"/bhns_t"$time".dat
sph="$idir"/sph_t"$time".dat

xc=`awk '{print $4;}' $bhns`
yc=`awk '{print $5;}' $bhns`
zc=`awk '{print $6;}' $bhns`

awk '{$4-=xc;$5-=yc;$6-=zc;if(-dz<$6&&$6<dz) print $0;}' dz=$delz xc=$xc yc=$yc zc=$zc $sph
