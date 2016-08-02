if test $# -ne 5
then
    echo "sh $0 <NSmass[Msun]> <snapWD> <rWD> <sep/Rlobe> <ofile>" >&2
    exit
fi

nsmas=$1
small=$2
rc=$3
clobe=$4
ofile=$5
bneps=1e6
gravc=0.000000066738480

m1=`echo "$nsmas * 1.9891 * 10^33" | bc`
m2=`awk 'BEGIN{mtot=0.0;}{mtot+=$3}END{printf("%lf\n", mtot);}' $small`

qq=`echo "scale=10; $m2 / $m1" | bc`

#rc=`awk '{print sqrt($4**2+$5**2+$6**2)}' $small | sort -g | tail -n1`

echo "0" > dummy

a0=`awk '{print clobe * rc * (0.6 * q**(2./3.) + log(1.+q**(1./3.))) / (0.49 * q**(2./3.))}' \
    clobe=$clobe rc=$rc q=$qq dummy`

r1=`awk '{print + a0 * q / (1. + q)}' a0=$a0 q=$qq dummy`
r2=`awk '{print - a0     / (1. + q)}' a0=$a0 q=$qq dummy`

v0=0.0
v1=`awk '{print + v0 * q / (1. + q)}' v0=$v0 q=$qq dummy`
v2=`awk '{print - v0     / (1. + q)}' v0=$v0 q=$qq dummy`

n2=`wc -l $small | tail -n1 | awk '{print $1}'`

echo "Mass ratio: $qq" >&2
echo "Secondary radius: $rc" >&2
echo "Binary separation: $a0" >&2
echo "Relative velocity: $v0" >&2
echo "Mass ratio: $qq" > "$ofile".log
echo "Secondary radius: $rc" >> "$ofile".log
echo "Binary separation: $a0" >> "$ofile".log
echo "Relative velocity: $v0" >> "$ofile".log

awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e\n", 0, 0, mns, dr, 0., 0., 0., dv, 0., eps);}' mns=$m1 dr=$r1 dv=$v1 eps=$bneps dummy > "$ofile".bhns
awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1+di, 1, $3, $4+dr, $5, $6, $7, $8+dv, $9, $13, $14, $15, $17, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44);}' di=$n1 dr=$r2 dv=$v2 $small > "$ofile".data

rm -f dummy

