if test $# -ne 6
then
    echo "sh $0 <BHmass[Msun]> <smallerWD> <rwd> <ri/rt> <beta> <ofile>" >&2
    exit
fi

bhmas=$1
ifile=$2
radius=$3
irirt=$4
ibeta=$5
ofile=$6

echo "sh $0 $bhmas $ifile $radius $irirt $ibeta $ofile" >&2
echo "sh $0 $bhmas $ifile $radius $irirt $ibeta $ofile" > "$ofile".log

bneps=1e6
gravc=0.000000066738480
msn=`echo "1.9891 * 10^33" | bc`
rwd=`echo "10^9" | bc`

m1=`echo "$bhmas * $msn" | bc`
m2=`awk 'BEGIN{mtot=0.0;}{mtot+=$3}END{printf("%lf\n", mtot);}' $ifile`
#rm=`awk 'BEGIN{r2max=0.0;}{r2=$4**2+$5**2+$6**2;if(r2>r2max)r2max=r2;}END{printf("%lf\n", sqrt(r2max));}' $ifile`
#rm=350000000
#rm=1100000000
rm=$radius
rt=`echo "scale=5; 1.2 * 10^11 * e(1/3.*l($m1/(1000000.*$msn))) * ($rm/$rwd) / e(1/3.*l($m2/(0.6*$msn)))" | bc -l`

rp=`echo "scale=5; $rt / $ibeta" | bc`
ri=`echo "scale=5; $rt * $irirt" | bc`
#hz=`echo "scale=5; e(0.5*l(2.*$gravc*$rp*($m1+$m2)))" | bc -l`
hz=`echo "scale=5; e(0.5*l(2.*$gravc*$rp*($m1)))" | bc -l`
vy=`echo "scale=5; $hz / $ri" | bc -l`
#vx=`echo "scale=5; e(0.5*l(2.*$gravc*($m1+$m2)/$ri-$vy*$vy))" | bc -l`
vx=`echo "scale=5; e(0.5*l(2.*$gravc*($m1)/$ri-$vy*$vy))" | bc -l`
v1=`echo "scale=5; e(0.5*l($vx*$vx+$vy*$vy))" | bc -l`

echo "0" > dummy

#e0=`awk '{print 1.+(-vy)*(+hz)/(g*(m1+m2))}' g=$gravc hz=$hz vx=$vx vy=$vy m1=$m1 m2=$m2 dummy`
e0=`awk '{print 1.+(-vy)*(+hz)/(g*(m1))}' g=$gravc hz=$hz vx=$vx vy=$vy m1=$m1 dummy`
#e1=`awk '{print   -(+vx)*(+hz)/(g*(m1+m2))}' g=$gravc hz=$hz vx=$vx vy=$vy m1=$m1 m2=$m2 dummy`
e1=`awk '{print   -(+vx)*(+hz)/(g*(m1))}' g=$gravc hz=$hz vx=$vx vy=$vy m1=$m1 dummy`
ee=`awk '{print sqrt(e0**2+e1**2);}' e0=$e0 e1=$e1 dummy`

cost=`awk '{print -e1}' e0=$e0 e1=$e1 dummy`
sint=`awk '{print -e0}' e0=$e0 e1=$e1 dummy`

dpx=`awk '{print cost * (-ri)}' cost=$cost sint=$sint ri=$ri dummy`
dpy=`awk '{print sint * (-ri)}' cost=$cost sint=$sint ri=$ri dummy`
dvx=`awk '{print cost * (+vx) - sint * (-vy)}' cost=$cost sint=$sint vx=$vx vy=$vy dummy`
dvy=`awk '{print sint * (+vx) + cost * (-vy)}' cost=$cost sint=$sint vx=$vx vy=$vy dummy`

printf "rt: %+e rp: %+e\n" $rt $rp >&2
printf "px: %+e py: %+e vx: %+e vy: %+e\n" $dpx $dpy $dvx $dvy >&2
awk '{printf("p/v: %+e\n", sqrt(px**2+py**2)/sqrt(vx**2+vy**2));}' px=$dpx py=$dpy vx=$dvx vy=$dvy \
    dummy >&2
echo "Caution! PW potential!" >&2
printf "rt: %+e rp: %+e\n" $rt $rp >> "$ofile".log
printf "px: %+e py: %+e vx: %+e vy: %+e\n" $dpx $dpy $dvx $dvy >> "$ofile".log
awk '{printf("p/v: %+e\n", sqrt(px**2+py**2)/sqrt(vx**2+vy**2));}' px=$dpx py=$dpy vx=$dvx vy=$dvy \
    dummy >> "$ofile".log

awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1, 0, $3, $4+dpx, $5+dpy, $6, $7+dvx, $8+dvy, $9, $13, $14, 0., $17, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44);}' dpx=$dpx dpy=$dpy dvx=$dvx dvy=$dvy $ifile > "$ofile".data

echo "$bhmas 1" > "$ofile".imbh

rm -f dummy

