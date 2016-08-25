if test $# -ne 6
then
    echo "sh $0 <NSmass[Msun]> <snapWD> <rWD> <sep/Rlobe> <ofile> <sync:0/nonr:1>" >&2
    exit
fi

nsmas=$1
small=$2
rc=$3
clobe=$4
ofile=$5
rflag=$6
bneps=1e6
gravc=0.000000066738480

echo "sh $0 $nsmas $small $rc $clobe $ofile $rflag" >&2
echo "sh $0 $nsmas $small $rc $clobe $ofile $rflag" > "$ofile".log

m1=`echo "$nsmas * 1.9891 * 10^33" | bc`
m2=`awk 'BEGIN{mtot=0.0;}{mtot+=$3}END{printf("%lf\n", mtot);}' $small`

qq=`echo "scale=10; $m2 / $m1" | bc`

echo "0" > dummy

a0=`awk '{print clobe * rc * (0.6 * q**(2./3.) + log(1.+q**(1./3.))) / (0.49 * q**(2./3.))}' \
    clobe=$clobe rc=$rc q=$qq dummy`

r1=`awk '{print + a0 * q / (1. + q)}' a0=$a0 q=$qq dummy`
r2=`awk '{print - a0     / (1. + q)}' a0=$a0 q=$qq dummy`

#v0=0.0
#v1=`awk '{print + v0 * q / (1. + q)}' v0=$v0 q=$qq dummy`
#v2=`awk '{print - v0     / (1. + q)}' v0=$v0 q=$qq dummy`

somg=0.0
if test $rflag -eq 1
then
    pa0=`printf "%f\n" $a0`
    somg=`echo "scale=10; $gravc * ($m1 + $m2) / $pa0^3" | bc -l`
fi

echo "Mass ratio: $qq" >&2
echo "Secondary radius: $rc" >&2
echo "Binary separation: $a0" >&2
echo "Sync=0 or Nrot=1: $rflag" >&2
echo "Mass ratio: $qq" > "$ofile".log
echo "Secondary radius: $rc" >> "$ofile".log
echo "Binary separation: $a0" >> "$ofile".log
echo "Sync=0 or Nrot=1: $rflag" >> "$ofile".log

awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e\n", 0, 0, mns, dr, 0., 0., 0., dv, 0., eps);}' mns=$m1 dr=$r1 dv=$v1 eps=$bneps dummy > "$ofile".bhns
#awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1+di, 1, $3, $4+dr, $5, $6, $7, $8+dv, $9, $13, $14, $15, $17, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44);}' di=$n1 dr=$r2 dv=$v2 $small > "$ofile".data
awk '{dvx=(-($5-0.)*somg); dvy=(-(-($4-dr)*somg)); printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1+di, 1, $3, $4+dr, $5, $6, $7+dvx, $8+dvy, $9, $13, $14, $15, $17, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44);}' di=$n1 dr=$r2 dv=$v2 somg=$somg $small > "$ofile".data

rm -f dummy

