if test $# -ne 2
then
    echo "sh $0 <snr?d> <N???>" >&2
    exit
fi

nmesh=$2
itype=`printf "%sN%03d" $1 $nmesh`
ifile1="$itype"E04.dat
ifile2="$itype"E14.dat
ifile3="$itype"E24.dat

filelist=" $ifile1 $ifile2 $ifile3"

for ifile in $filelist
do
    if ! test -e $ifile
    then
        echo "Error: $ifile is not found." >&2
        exit
    else
        echo "$ifile exists." >&2
    fi
done

tdir=`dirname $0`

echo "start dens..." >&2
otype=`printf "densN%03d" $nmesh`
sh "$tdir"/transfer1dDataFromTextToNumpyArray.sh $ifile1 4 $otype

exit

echo "start ni56..." >&2
otype=`printf "ni56N%03d" $nmesh`
sh "$tdir"/transfer1dDataFromTextToNumpyArray.sh $ifile1 5 $otype

AtomicNumber=1
for ifile in $filelist
do
    for icol in $(seq 6 1 15)
    do
        echo "start Z=$AtomicNumber..." >&2
        otype=`printf "z%03dN%03d" $AtomicNumber $nmesh`
        sh "$tdir"/transfer1dDataFromTextToNumpyArray.sh $ifile $icol $otype
        AtomicNumber=`echo "$AtomicNumber+1" | bc`
    done
done

