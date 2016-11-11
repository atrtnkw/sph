if test $# -ne 7
then
    echo "sh $0 <NSmass[Msun]> <snapWD> <rwd> <sep/Rlobe> <rp> <ecc> <ofile>" >&2
    exit
fi

echo "0" > dummy

nsmas=$1
small=$2
rwd=`awk '{printf("%lf\n", rwd);}' rwd=$3 dummy`
clobe=$4
rp=`awk '{printf("%lf\n", rp);}' rp=$5 dummy`
ecc=$6
ofile=$7

bneps=1e6
gravc=0.000000066738480

echo "sh $0 $nsmas $small $rwd $clobe $rp $ecc $ofile" >&2
echo "sh $0 $nsmas $small $rwd $clobe $rp $ecc $ofile" > "$ofile".log

m1=`echo "$nsmas * 1.9891 * 10^33" | bc -l`
m2=`awk 'BEGIN{mtot=0.0;}{mtot+=$3}END{printf("%lf\n", mtot);}' $small`
a0=`echo "$rp / (1.0 - $ecc)" | bc -l`

qq=`echo "$m2 / $m1" | bc -l`
qq13=`echo "e(1./3.*l($qq))" | bc -l`
qq23=`echo "e(2./3.*l($qq))" | bc -l`

ri=`echo "$clobe * $rwd * (0.6 * $qq23 + l(1. + $qq13)) /  (0.49 * $qq23)" | bc -l`
ra=`echo "$a0 * (1 + $ecc)" | bc -l`
bl=`echo "$ri > $ra" | bc -l`
if test $bl -eq 1
then
    echo "!!!!!!!!!!!!!!!!!" >&2
    echo "ERROR!!" >&2
    echo "Initial binary seperation is too distant!" >&2
    printf "ri: %+e ra: %+e\n" $ri $ra >&2
    echo "!!!!!!!!!!!!!!!!!" >&2
    exit
fi
vi=`echo "sqrt(2.0 * $gravc * ($m1 + $m2) * (1/$ri - 1/(2.0*$a0)))" | bc -l`
nn=`echo "sqrt($gravc * ($m1 + $m2) / ($a0 * $a0 * $a0))" | bc -l`
cosu=`echo "(1. - $ri / $a0) / $ecc" | bc -l`
sinu=`echo "- sqrt(1. - $cosu * $cosu)" | bc -l`

pxi=`echo "$a0 * ($cosu - $ecc)" | bc -l`
pyi=`echo "$a0 * sqrt(1. - $ecc * $ecc) * $sinu" | bc -l`
vxi=`echo "- $a0 * $a0 * $nn * $sinu / $ri" | bc -l`
vyi=`echo "$a0 * $a0 * $nn * sqrt(1. - $ecc * $ecc) * $cosu / $ri" | bc -l`

px1=`echo "- $pxi * $qq / (1.0 + $qq)" | bc -l`
py1=`echo "- $pyi * $qq / (1.0 + $qq)" | bc -l`
vx1=`echo "- $vxi * $qq / (1.0 + $qq)" | bc -l`
vy1=`echo "- $vyi * $qq / (1.0 + $qq)" | bc -l`
px2=`echo "- $pxi * -1. / (1.0 + $qq)" | bc -l`
py2=`echo "- $pyi * -1. / (1.0 + $qq)" | bc -l`
vx2=`echo "- $vxi * -1. / (1.0 + $qq)" | bc -l`
vy2=`echo "- $vyi * -1. / (1.0 + $qq)" | bc -l`

prd=`echo "2.0 * 3.141 * sqrt($a0^3 / ($gravc * ($m1 + $m2)))" | bc -l`

printf "Mass ratio: %+e\n" $qq           >&2
printf "Pericenter distance: %+e\n" $rpp >&2
printf "Binary separation: %+e\n" $ri    >&2
printf "Binary velocity: %+e\n" $vi      >&2
printf "Binary period: %+e\n" $prd       >&2
printf "Mass ratio: %+e\n" $qq           >> "$ofile".log
printf "Pericenter distance: %+e\n" $rpp >> "$ofile".log
printf "Binary separation: %+e\n" $ri    >> "$ofile".log
printf "Binary velocity: %+e\n" $vi      >> "$ofile".log
printf "Binary period: %+e\n" $prd       >> "$ofile".log

awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e\n", 0, 0, mns, px, py, 0., vx, vy, 0., eps);}' mns=$m1 px=$px1 py=$py1 vx=$vx1 vy=$vy1 eps=$bneps dummy > "$ofile".bhns
awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1, 1, $3, $4+px, $5+py, $6, $7+vx, $8+vy, $9, $13, $14, $15, $17, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44);}' px=$px2 py=$py2 vx=$vx2 vy=$vy2 $small > "$ofile".data

rm -f dummy

