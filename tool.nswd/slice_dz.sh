time=0200
idir=run
delz=5e7

bhns="$idir"/bhns_t"$time".dat
sph="$idir"/sph_t"$time".dat

xc=`awk '{print $4;}' $bhns`
yc=`awk '{print $5;}' $bhns`
zc=`awk '{print $6;}' $bhns`

awk '{$4-=xc;$5-=yc;$6-=zc;if(-dz<$6&&$6<dz) print $0;}' dz=$delz xc=$xc yc=$yc zc=$zc $sph
