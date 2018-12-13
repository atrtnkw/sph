if test $# -ne 2
then
    echo "sh $0 <idir> <Z/Zsun>" >&2
    exit
fi

idir=$1
metl=$2
odir=data
zsun=0.012

if test -e $odir
then
    echo "Error: Directory $odir exists" >&2
    exit
fi

xne22=`echo "1.11 * $zsun" | bc -l`
xc12=`echo "(1 - $xne22) / 2" | bc -l`
xo16=$xc12

files=`ls "$idir"/id*.in`

if ! test -e $odir
then
    mkdir -p $odir
fi

for file in $files
do
    nc=`echo $file | wc | awk '{print $3}'`
    cb=`echo "$nc - 13" | bc`
    cf=`echo "$nc - 4"  | bc`
    id=`echo $file | cut -c $cb-$cf`
    echo -e ""   > input
    echo -e "6" >> input
    echo -e "9" >> input
    echo -e ""  >> input
    echo -e "$file" >> input
    echo -e "3" >> input
    echo -e "0. 0. 0. $xc12 0. 0. $xo16 0. $xne22 0. 0.0 0. 0. 0. 0. 0. 0. 0." >> input
    echo -e "data/id$id.out_" >> input
    echo -e "-1" >> input
    ./run < input
done
