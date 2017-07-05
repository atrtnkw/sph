if test $# -ne 1
then
    echo "sh $0 <idir>"
    exit
fi

idir=$1

if ! test -e list
then
    mkdir list
fi

awk '{printf("%16.10f %+e %+e\n", $1, $3, $4) > sprintf("list/co.dat.%010d", $2);}' \
    "$idir"/debug_p*.log

exit

lfile=`ls list/`
for file in $lfile
do
    sort -g $lfile > tmp.sort
    mv tmp.sort > $lfile
done

rm -f tmp.sort
