if test $# -ne 2
then
    echo "sh $0 <idir> <odir>"
    exit
fi

idir=$1
odir=$2

if ! test -e $odir
then
    mkdir $odir
fi

idlist=`awk '{if($1==0.) print $2;}' "$idir"/debug*.log | sort -n`

for id in $idlist
do
    pid=`printf "%010d" $id`
    awk '{if($2==id) print $1, $3, $4}' id=$id "$idir"/debug*.log > "$odir"/history.dat."$pid"
done
