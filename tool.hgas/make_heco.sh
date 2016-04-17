if test $# -ne 3
then
    echo "sh $0 <InputFile> <FractionOfHeOverall> <MassFractionOfHePerParticle>"
    exit
fi

ifile=$1
hefrc=$2
hemfrc=$3

nptcl=`wc -l $ifile | awk '{print $1}'`
nhelium=`echo "$nptcl * $hefrc / $hemfrc" | bc -l | awk '{printf("%d\n", $1);}'`
rmin=`awk '{print $4**2+$5**2+$6**2;}' $ifile | sort -g \
    | tail -n"$nhelium" | head -n1 | awk '{print sqrt($1);}'`

awk '{r2=$4**2+$5**2+$6**2;if(r2>=rmin*rmin){co=0.5*(1.-he);$32=he;$33=co;$34=co;} print $0;}' \
    rmin=$rmin he=$hemfrc $ifile
