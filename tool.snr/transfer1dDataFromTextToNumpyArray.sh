if test $# -ne 2
then
    echo "sh $0 <ifile> <odir>" >&2
    exit
fi

ityp=$1
odir=$2

tdir=`dirname $0`
ftype=`basename --suffix=.sh $0`

alias python3='~/work/yt-conda/bin/python3'

ifile="$ityp"E04.dat
python3 "$tdir"/"$ftype".py $ifile 3 $odir/dens
python3 "$tdir"/"$ftype".py $ifile 4 $odir/ni56

for znum in $(seq 1 1 30)
do
    if test $znum -le 10
    then
        ifile="$ityp"E04.dat
        icol=`echo "$znum+4" | bc`
    elif test $znum -le 20
    then
        ifile="$ityp"E14.dat
        icol=`echo "$znum+4-10" | bc`
    else
        ifile="$ityp"E24.dat
        icol=`echo "$znum+4-20" | bc`
    fi
    pnum=`printf "%03d" $znum`
    python3 "$tdir"/"$ftype".py $ifile $icol $odir/z"$pnum"
done

