if test $# -ne 2
then
    echo "sh $0 <file> <fhe>"
    exit
fi

file=$1
frac=$2

awk '{if($2==0) print $0;}' $file \
    | awk 'BEGIN{x=y=z=0.0;m=0.0;n=0;}
{x+=$3*$4;y+=$3*$5;z+=$3*$6;m+=$3;n++;}
END{printf("%+.16e %+.16e %+.16e %d\n", x/m, y/m, z/m, n);}' \
    > tmp.xc
x0=`awk '{print $1}' tmp.xc`
y0=`awk '{print $2}' tmp.xc`
z0=`awk '{print $3}' tmp.xc`
n0=`awk '{print $4}' tmp.xc`

awk '{if($2==1) print $0;}' $file \
    | awk 'BEGIN{x=y=z=0.0;m=0.0;n=0;}
{x+=$3*$4;y+=$3*$5;z+=$3*$6;m+=$3;n++;}
END{printf("%+.16e %+.16e %+.16e %d\n", x/m, y/m, z/m, n);}' \
    > tmp.xc
x1=`awk '{print $1}' tmp.xc`
y1=`awk '{print $2}' tmp.xc`
z1=`awk '{print $3}' tmp.xc`
n1=`awk '{print $4}' tmp.xc`

nhe0=`echo "scale=0; $n0 * $frac" | bc`
nhe0=`printf "%.f\n" $nhe0`
nhe1=`echo "scale=0; $n1 * $frac" | bc`
nhe1=`printf "%.f\n" $nhe1`

awk '{if($2==0){$4-=dx;$5-=dy;$6-=dz;printf("%10d %+e\n", $1, $4**2+$5**2+$6**2);}}' \
    dx=$x0 dy=$y0 dz=$z0 $file \
    | sort -g -k 2 | tail -n"$nhe0" | sort -n | awk '{print $1}'

awk '{if($2==1){$4-=dx;$5-=dy;$6-=dz;printf("%10d %+e\n", $1, $4**2+$5**2+$6**2);}}' \
    dx=$x1 dy=$y1 dz=$z1ls $file \
    | sort -g -k 2 | tail -n"$nhe1" | sort -n | awk '{print $1}'

rm -f tmp.*

