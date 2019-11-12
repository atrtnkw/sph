if test $# -ne 1
then
    echo "sh $0 <ifile>"
    exit
fi

ifile=$1

if ! test -e $ifile
then
    echo "Error: $ifile is Not found." >&2
    exit
fi

tdir=`dirname $0`
ftype=`basename --suffix=.sh $0`

alias python3='~/work/yt-conda/bin/python3'
python3 "$tdir"/"$ftype".py $ifile

