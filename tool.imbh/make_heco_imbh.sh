if test $# -ne 3
then
    echo "sh $0 <InputFile> <FractionOfHe> <FractionOfMixing>"
    exit
fi

ifile=$1
frche=$2
frcmx=$3
dtool=~/git-sph/tool.hgas

awk 'BEGIN{x=y=z=0.;n=0;}{x+=$4;y+=$5;z+=$6;n++;}END{print x/n, y/n, z/n, n;}' $ifile > tmp.cen
x=`awk '{print $1;}' tmp.cen`
y=`awk '{print $2;}' tmp.cen`
z=`awk '{print $3;}' tmp.cen`
n=`awk '{print $4;}' tmp.cen`

ncrit=`echo "scale=16; $n * (1. - $frche / $frcmx)" | bc`

r2crit=`awk '{print ($4-x)**2+($5-y)**2+($6-z)**2;}' x=$x y=$y z=$z $ifile \
    | sort -g \
    | awk '{if(NR>=ncrit) printf("%+e\n", $1);}' ncrit=$ncrit \
    | head -n1`

awk '{r2=($4-x)**2+($5-y)**2+($6-z)**2;if(r2>=r2crit){$14=frcmx;$15=0.5*(1-frcmx);$16=0.5*(1-frcmx);} print $0;}' \
    x=$x y=$y z=$z r2crit=$r2crit frcmx=$frcmx $ifile \
    | awk '{printf("%8d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26);}'

rm -f tmp.*

