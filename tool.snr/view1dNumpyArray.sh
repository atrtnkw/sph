if test $# -ne 2
then
    echo "sh $0 <ifile> <nmesh>"
    exit
fi

ifile=$1
nmesh=$2

if ! test -e $ifile
then
    echo "Error: $ifile is Not found." >&2
    exit
fi

tdir=`dirname $0`
ftype=`basename --suffix=.sh $0`

alias python3='~/work/yt-conda/bin/python3'
if ! test -e hoge.dat1
then
    python3 "$tdir"/"$ftype".py $ifile > hoge.dat1
fi

if ! test -e hoge.dat2
then
    awk 'BEGIN{RS=OFS;}{print $1}' hoge.dat1 \
        | awk '{if(NR!=1&&NF!=0) print $0}' \
        > hoge.dat2
fi

nline=`wc -l hoge.dat2 | awk '{print $1}'`
nmesh2=`echo "$nmesh * $nmesh         " | bc`
nmesh3=`echo "$nmesh * $nmesh * $nmesh" | bc`
if test $nline -eq $nmesh2
then
    awk '{i=int((NR-1)%nmesh);j=int((NR-1)/nmesh);printf("%4d %4d %4d %+e\n", i, j, 0, $1);}' nmesh=$nmesh hoge.dat2
elif test $nline -eq $nmesh3
then
    awk '{i=(NR-1)%nmesh;j=int((NR-1)/nmesh)%nmesh;k=int((NR-1)/(nmesh*nmesh));printf("%4d %4d %4d %+e\n", i, j, k, $1)}' \
        nmesh=$nmesh hoge.dat2
else
    echo "Error: nmesh is wrong." >&2
    exit
fi

rm -f hoge.*
