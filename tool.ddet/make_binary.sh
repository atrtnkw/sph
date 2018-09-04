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

if ! test -e $large
then
    echo "Error: Not found $large" >&2
    exit
fi

if ! test -e $small
then
    echo "Error: Not found $small" >&2
    exit
fi

touch $ftype.log
if ! test -e $ftype.log
then
    echo "Error: Not found $ftype.log" >&2
    exit
fi

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

#awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1+di, 0, $3, $4+dr, $5, $6, $7, $8+dv, $9, $10, $11, $12, $13, $14+dX, $15+dX, $16+dX, $17+dX, $18+dX, $19+dX, $20+dX, $21+dX, $22+dX, $23+dX, $24+dX, $25+dX, $26+dX);}' di=0   dr=$r1 dv=$v1 $large dX=1e-10 > "$ftype".data
awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1+di, 0, $3, $4+dr, $5, $6, $7, $8+dv, $9, $10, $11, $12, $13, $14+1e-10, $15+1e-10, $16+1e-10, $17+1e-10, $18+1e-10, $19+1e-10, $20+1e-10, $21+1e-10, $22+1e-10, $23+1e-10, $24+1e-10, $25+1e-10, $26+1e-10);}' di=0 dr=$r1 dv=$v1 $large > "$ftype".data
awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1+di, 1, $3, $4+dr, $5, $6, $7, $8+dv, $9, $13, $14, $15, $17, $32+1e-10, $33+1e-10, $34+1e-10, $35+1e-10, $36+1e-10, $37+1e-10, $38+1e-10, $39+1e-10, $40+1e-10, $41+1e-10, $42+1e-10, $43+1e-10, $44+1e-10);}' di=$n1 dr=$r2 dv=$v2 $small dX=1e-10 >> "$ftype".data

rm -f dummy


