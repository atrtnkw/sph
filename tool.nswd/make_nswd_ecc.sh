if test $# -ne 5
then
    echo "sh $0 <NSmass[Msun]> <snapWD> <rp/1e9cm> <ecc> <ofile>" >&2
    exit
fi

nsmas=$1
small=$2
rp=$3
ecc=$4
ofile=$5
bneps=1e6
gravc=0.000000066738480

echo "sh $0 $nsmas $small $rp $ecc $ofile" >&2
echo "sh $0 $nsmas $small $rp $ecc $ofile" > "$ofile".log

m1=`echo "$nsmas * 1.9891 * 10^33" | bc -l`
m2=`awk 'BEGIN{mtot=0.0;}{mtot+=$3}END{printf("%lf\n", mtot);}' $small`
a0=`echo "$rp * 10^9 / (1.0 - $ecc)" | bc -l`
r0=`echo "$a0 * (1.0 + $ecc)" | bc -l`
v0=`echo "sqrt($gravc * ($m1 + $m2) / $a0 * (1.0 - $ecc) / (1.0 + $ecc))" | bc -l`

qq=`echo "$m2 / $m1" | bc -l`

r1=`echo "$r0 * $qq / (1.0 + $qq)" | bc -l`
r2=`echo "$r0 * -1. / (1.0 + $qq)" | bc -l`
v1=`echo "$v0 * $qq / (1.0 + $qq)" | bc -l`
v2=`echo "$v0 * -1. / (1.0 + $qq)" | bc -l`

rpp=`echo "$rp * 10^9" | bc -l`
prd=`echo "2.0 * 3.141 * sqrt($a0^3 / ($gravc * ($m1 + $m2)))" | bc -l`

printf "Mass ratio: %+e\n" $qq           >&2
printf "Pericenter distance: %+e\n" $rpp >&2
printf "Binary separation: %+e\n" $r0    >&2
printf "Binary period: %+e\n" $prd    >&2
printf "Mass ratio: %+e\n" $qq           >> "$ofile".log
printf "Pericenter distance: %+e\n" $rpp >> "$ofile".log
printf "Binary separation: %+e\n" $r0    >> "$ofile".log
printf "Binary period: %+e\n" $prd       >> "$ofile".log

echo "0" > dummy

awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e\n", 0, 0, mns, dr, 0., 0., 0., dv, 0., eps);}' mns=$m1 dr=$r1 dv=$v1 eps=$bneps dummy > "$ofile".bhns
awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1, 1, $3, $4+dr, $5, $6, $7, $8+dv, $9, $13, $14, $15, $17, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44);}' dr=$r2 dv=$v2 somg=$somg $small > "$ofile".data

rm -f dummy

