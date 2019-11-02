if test $# -ne 3
then
    echo "sh $0 <idir> <xcol> <dx>" >&2
    exit
fi
idir=$1
pxxx=$2
dxxx=$3
tdir=`dirname $0`

filelist=`ls $idir/mesh*.dat`

for file in $filelist
do
    awk -f "$tdir"/extractColumn1d.awk col=-1 xcol=$pxxx dx=$dxxx $file > $file.x"$pxxx"
done
