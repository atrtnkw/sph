if test $# -ne 1
then
    echo "sh $0 <Z/Zsun>" >&2
    exit
fi

metl=$1
idir=init
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
    id=`echo $file | cut -c 8-17`
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
