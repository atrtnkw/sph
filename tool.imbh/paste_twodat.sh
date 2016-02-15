if test $# -ne 4
then
    echo "sh $0 <idir> <odir> <time0> <time1>"
    exit
fi

idir=$1
odir=$2
tim0=`printf "%04d\n" $3`
tim1=`printf "%04d\n" $4`

ifile0="$idir"/sph_t"$tim0".dat
ifile1="$idir"/sph_t"$tim1".dat
ofile="$odir"/sph_t"$tim0"-"$tim1".dat

sort -n $ifile0 > tmp0.dat
sort -n $ifile1 > tmp1.dat

awk 'BEGIN{nl=0}{dn=$1-nl;for(il=0;il<dn;il++){printf("\n");}print $0;nl=$1+1;}' tmp0.dat > tmp0.cor
awk 'BEGIN{nl=0}{dn=$1-nl;for(il=0;il<dn;il++){printf("\n");}print $0;nl=$1+1;}' tmp1.dat > tmp1.cor

paste tmp0.cor tmp1.cor > $ofile

rm tmp?.*
