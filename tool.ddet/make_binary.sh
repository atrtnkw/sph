if test $# -ne 5
then
    echo "sh $0 <largerWD> <smallerWD> <sep/Rlobe> <rcomp/cm> <ofile>" >&2
    exit
fi

large=$1
small=$2
clobe=$3
rcomp=$4
ftype=$5
gravc=0.000000066738480

m1=`awk 'BEGIN{mtot=0.0;}{mtot+=$3}END{printf("%lf\n", mtot);}' $large`
m2=`awk 'BEGIN{mtot=0.0;}{mtot+=$3}END{printf("%lf\n", mtot);}' $small`

qq=`echo "scale=10; $m2 / $m1" | bc`

#rc=`awk '{print sqrt($4**2+$5**2+$6**2)}' $small | sort -g | tail -n1`
rc=$rcomp

echo "0" > dummy

a0=`awk '{print clobe * rc * (0.6 * q**(2./3.) + log(1.+q**(1./3.))) / (0.49 * q**(2./3.))}' \
    clobe=$clobe rc=$rc q=$qq dummy`

r1=`awk '{print + a0 * q / (1. + q)}' a0=$a0 q=$qq dummy`
r2=`awk '{print - a0     / (1. + q)}' a0=$a0 q=$qq dummy`

v0=`awk '{print sqrt(g * (m1 + m2) / a);}' g=$gravc m1=$m1 m2=$m2 a=$a0 dummy`
v1=`awk '{print + v0 * q / (1. + q)}' v0=$v0 q=$qq dummy`
v2=`awk '{print - v0     / (1. + q)}' v0=$v0 q=$qq dummy`

np=`wc -l $large $small | tail -n1 | awk '{print $1}'`
n1=`wc -l $large        | tail -n1 | awk '{print $1}'`

echo "Mass ratio: $qq"        >&2
echo "Secondary radius: $rc"  >&2
echo "Binary separation: $a0" >&2
echo "Relative velocity: $v0" >&2
echo "sh $0 $1 $2 $3 $4 $5"    > "$ftype".log
echo "Mass ratio: $qq"        >> "$ftype".log
echo "Secondary radius: $rc"  >> "$ftype".log
echo "Binary separation: $a0" >> "$ftype".log
echo "Relative velocity: $v0" >> "$ftype".log

awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1+di, 0, $3, $4+dr, $5, $6, $7, $8+dv, $9, $10, $11, $12, $13, $14+dX, $15+dX, $16+dX, $17+dX, $18+dX, $19+dX, $20+dX, $21+dX, $22+dX, $23+dX, $24+dX, $25+dX, $26+dX);}' di=0   dr=$r1 dv=$v1 $large dX=1e-10 > "$ftype".data
awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1+di, 1, $3, $4+dr, $5, $6, $7, $8+dv, $9, $13, $14, $15, $17, $32+dX, $33+dX, $34+dX, $35+dX, $36+dX, $37+dX, $38+dX, $39+dX, $40+dX, $41+dX, $42+dX, $43+dX, $44+dX);}' di=$n1 dr=$r2 dv=$v2 $small dX=1e-10 >> "$ftype".data

rm -f dummy


