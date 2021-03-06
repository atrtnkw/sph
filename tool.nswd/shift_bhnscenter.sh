if test $# -ne 2
then
    echo "sh $0 <idir> <time>" >&2
    exit
fi

idir=$1
time=$2

ptim=`printf "%04d" $time`
sfile="$idir"/sph_t"$ptim".dat
tfile="$idir"/sph_t"$ptim"_*_i000000.dat
bfile="$idir"/bhns_t"$ptim".dat
#if ! test -e $sfile
if ! test -e $sfile -o $tfile
then
    echo "File $sfile is not found" >&2
    exit
fi
if ! test -e $bfile
then
    echo "File $bfile is not found" >&2
    exit
fi
bx=`awk '{print $4;}' $bfile`
by=`awk '{print $5;}' $bfile`
bz=`awk '{print $6;}' $bfile`

#awk '{$4-=bx;$5-=by;$6-=bz;print $0;}' bx=$bx by=$by bz=$bz $sfile
awk '{$4-=bx;$5-=by;$6-=bz;print $0;}' bx=$bx by=$by bz=$bz $idir/sph_t"$ptim"*.dat
