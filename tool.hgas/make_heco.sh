if test $# -ne 2
then
    echo "sh $0 <InputFile> <FractionOfHe>"
    exit
fi

ifile=$1
frche=$2
dtool=~/git-sph/tool.hgas

awk -f "$dtool"/change_t2t.awk $ifile > tmp.data

awk 'BEGIN{x=y=z=0.;n=0;}{x+=$4;y+=$5;z+=$6;n++;}END{print x/n, y/n, z/n, n;}' tmp.data > tmp.cen
x=`awk '{print $1;}' tmp.cen`
y=`awk '{print $2;}' tmp.cen`
z=`awk '{print $3;}' tmp.cen`
n=`awk '{print $4;}' tmp.cen`

ncrit=`echo "$n * (1. - $frche)" | bc`

r2crit=`awk '{print ($4-x)**2+($5-y)**2+($6-z)**2;}' x=$x y=$y z=$z tmp.data \
    | sort -g \
    | awk '{if(NR>=ncrit) printf("%+e\n", $1);}' ncrit=$ncrit \
    | head -n1`

awk '{r2=($4-x)**2+($5-y)**2+($6-z)**2;if(r2>=r2crit){$14=1.;$15=0.;$16=0.;} print $0;}' \
    x=$x y=$y z=$z r2crit=$r2crit tmp.data \
    | awk '{printf("%8d %2d %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26);}'

rm -f tmp.*

