if test $# -ne 2
then
    echo "sh $0 <idir> <time:????>" >&2
    exit
fi

idir=$1
time=$2

wbin=`ruby -e 'puts 10**0.05'`
bhns="$idir"/bhns_t"$time".dat
sph="$idir"/sph_t"$time".dat
xc=`awk '{print $4;}' $bhns`
yc=`awk '{print $5;}' $bhns`
zc=`awk '{print $6;}' $bhns`

rad0=`ruby -e 'puts 10**8'`
for i in $(seq 1 1 40)
do
    rad1=`echo "scale=5; $rad0 * $wbin" | bc`

    awk 'BEGIN{m=0;n=0;}{$4-=xc;$5-=yc;$6-=zc;if(r0**2<$4**2+$5**2&&$4**2+$5**2<r1**2){m+=$3;n++}}END{printf("%+e %+e %+e %8d\n", sqrt(r0*r1), m/(4.*3.1415927*(r1**2-r0**2)), m, n);}' \
        r0=$rad0 r1=$rad1 xc=$xc yc=$yc zc=$zc $sph

    rad0=$rad1
done

