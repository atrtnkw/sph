if test $# -ne 3
then
    echo "sh $0 <idir> <dx> <col>" >&2
    exit
fi
idir=$1
dxxx=$2
ncol=$3
tdir=`dirname $0`

filelist=`ls $idir/mesh*.dat`

for file in $filelist
do
    awk -f "$tdir"/extractColumn1d.awk col=$ncol dx=$dxxx $file > $file.col"$ncol"
done
