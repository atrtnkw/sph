if test $# -ne 4
then
    echo "sh $0 <BHmass[Msun]> <smallerWD> <beta> <ofile>"
    exit
fi

bhmas=$1
ifile=$2
ibeta=$3
ofile=$4

bneps=1e6
gravc=0.000000066738480
msn=`echo "1.9891 * 10^33" | bc`
rwd=`echo "10^9" | bc`

m1=`echo "$bhmas * $msn" | bc`
m2=`awk 'BEGIN{mtot=0.0;}{mtot+=$3}END{printf("%lf\n", mtot);}' $ifile`
rm=`awk 'BEGIN{r2max=0.0;}{r2=$4**2+$5**2+$6**2;if(r2>r2max)r2max=r2;}END{printf("%lf\n", sqrt(r2max));}' $ifile`
rt=`echo "scale=5; 1.2 * 10^11 * e(1/3.*l($m1/(1000000.*$msn))) * ($rm/$rwd) / e(1/3.*l($m2/(0.6*$msn)))" | bc -l`

rp=`echo "scale=5; $rt / $ibeta" | bc`
ri=`echo "scale=5; $rt * 3." | bc`
vy=`echo "scale=5; e(0.5*l(2.*$gravc*$rp*($m1+$m2))) / $ri" | bc -l`
vx=`echo "scale=5; e(0.5*l(2.*$gravc*($m1+$m2)/$ri-$vy*$vy))" | bc -l`

echo "0" > dummy

r1=`awk '{print + ri * m2 / (m1 + m2)}' ri=$ri m1=$m1 m2=$m2 dummy`
r2=`awk '{print - ri * m1 / (m1 + m2)}' ri=$ri m1=$m1 m2=$m2 dummy`
vx1=`awk '{print - vx * m2 / (m1 + m2)}' vx=$vx m1=$m1 m2=$m2 dummy`
vx2=`awk '{print + vx * m1 / (m1 + m2)}' vx=$vx m1=$m1 m2=$m2 dummy`
vy1=`awk '{print + vy * m2 / (m1 + m2)}' vy=$vy m1=$m1 m2=$m2 dummy`
vy2=`awk '{print - vy * m1 / (m1 + m2)}' vy=$vy m1=$m1 m2=$m2 dummy`

printf "rtrp: %+e %+e\n" $rt $rp > "$ofile".log
printf "IMBH: %+e %+e %+e\n" $r1 $vx1 $vy1 >> "$ofile".log
printf "WD:   %+e %+e %+e\n" $r2 $vx2 $vy2 >> "$ofile".log

awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e\n", 0, 0, mns, x, 0., 0., vx, vy, 0., eps);}' mns=$m1 x=$r1 vx=$vx1 vy=$vy1 eps=$bneps dummy > "$ofile".bhns
awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1, 1, $3, $4+dx, $5, $6, $7+dvx, $8+dvy, $9, $13, $14, $15, $17, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44);}' dx=$r2 dvx=$vx2 dvy=$vy2 $ifile > "$ofile".data

cp "$ofile".bhns  "$ofile".atc0.bhns
awk -f ~/git-sph/tool.hgas/zero_alphu.awk "$ofile".data > "$ofile".atc0.data

rm -f dummy

