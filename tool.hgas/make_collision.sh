m1=1.2
m2=1.2
mass=`echo "$m1+$m2" | bc -l`

if test $# -ne 4
then
    echo "sh $0 <ifile[OnlyTheSameWD]> <InitialSeparation[cm]> <PericenterDistance[cm]> <VelocityAtInfinity[cm/s]>" >&2
    exit
fi

file=$1
rini=$2
rpri=$3
vinf=$4

tdir=`dirname $0`

touch dummy

vini=`awk 'BEGIN{msun=1.989e33;grav=6.674e-8;}END{v0=(rini**2*vinf**2+2.*grav*mass*msun*rini)**0.5/rini; print v0;}' \
rini=$rini rpri=$rpri vinf=$vinf mass=$mass dummy`
sini=`awk 'BEGIN{msun=1.989e33;grav=6.674e-8;}END{s0=((rpri**2*vinf**2+2.*grav*mass*msun*rpri)/(rini**2*vinf**2+2.*grav*mass*msun*rini))**0.5; print s0;}' \
rini=$rini rpri=$rpri vinf=$vinf mass=$mass dummy`
r1x=`awk 'END{r1x=+m2/mass*rini; print r1x}' mass=$mass m1=$m1 m2=$m2 rini=$rini vini=$vini sini=$sini dummy`
v1x=`awk 'END{cini=(1.-sini**2)**0.5;v1x=-m2/mass*vini*cini; print v1x}' mass=$mass m1=$m1 m2=$m2 rini=$rini vini=$vini sini=$sini dummy`
v1y=`awk 'END{cini=(1.-sini**2)**0.5;v1y=+m2/mass*vini*sini; print v1y}' mass=$mass m1=$m1 m2=$m2 rini=$rini vini=$vini sini=$sini dummy`
r2x=`awk 'END{r2x=-m1/mass*rini; print r2x}' mass=$mass m1=$m1 m2=$m2 rini=$rini vini=$vini sini=$sini dummy`
v2x=`awk 'END{cini=(1.-sini**2)**0.5;v2x=+m1/mass*vini*cini; print v2x}' mass=$mass m1=$m1 m2=$m2 rini=$rini vini=$vini sini=$sini dummy`
v2y=`awk 'END{cini=(1.-sini**2)**0.5;v2y=-m1/mass*vini*sini; print v2y}' mass=$mass m1=$m1 m2=$m2 rini=$rini vini=$vini sini=$sini dummy`

echo "m1=$m1 [Msun]" >&2
echo "m2=$m2 [Msun]" >&2
echo "ri=$rini [cm]" >&2
echo "rp=$rpri [cm]" >&2
echo "vinf=$vinf [cm/s]" >&2
echo "vi=$vini [cm/s]" >&2
echo "sintheta=$sini" >&2
printf "r1x=%+e v1x=%+e v1y=%+e\n" $r1x $v1x $v1y >&2
printf "r2x=%+e v2x=%+e v2y=%+e\n" $r2x $v2x $v2y >&2

n1=`wc -l $file | awk '{print $1}'`

awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1+di, 0, $3, $4+drx, $5, $6, $7+dvx, $8+dvy, $9, $13, $14, $15, $17, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44);}' di=0 drx=$r1x dvx=$v1x dvy=$v1y $file
awk '{printf("%10d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1+di, 1, $3, $4+drx, $5, $6, $7+dvx, $8+dvy, $9, $13, $14, $15, $17, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44);}' di=$n1 drx=$r2x dvx=$v2x dvy=$v2y $file

rm -f dummy
