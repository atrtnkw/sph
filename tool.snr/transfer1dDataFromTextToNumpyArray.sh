if test $# -ne 3
then
    echo "sh $0 <ifile> <icolumn=1:x,2:y,3:z,4:dn...> <otype>" >&2
    exit
fi

ifile=$1
icolumn=`echo "$2-1" | bc`
otype=$3

tdir=`dirname $0`
ftype=`basename --suffix=.sh $0`

if test -e $otype.npy
then
    echo "Error: $otype.npy exists." >&2
    exit
fi

alias python3='~/work/yt-conda/bin/python3'
python3 "$tdir"/"$ftype".py $ifile $icolumn $otype
