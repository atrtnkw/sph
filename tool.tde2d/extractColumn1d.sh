if test $# -ne 2
then
    echo "sh $0 <idir> <col>" >&2
    exit
fi
idir=$1
ncol=$2
tdir=`dirname $0`

filelist=`ls $idir/mesh*.dat`

for file in $filelist
do
    awk -f "$tdir"/extractColumn1d.awk col=$ncol $file > $file.col"$ncol"
done
