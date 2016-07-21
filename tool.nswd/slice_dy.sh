if test $# -ne 4
then
    echo "sh $0 <HOME> <idir> <time=xxxx> <dy>" >&2
    exit
fi

gdir="$1"/git-sph
idir=$2
time=$3
dy=$4

sfile="$idir"/sph_t"$time".dat
bfile="$idir"/bhns_t"$time".dat

px=`awk '{print $4;}' $bfile`
py=`awk '{print $5;}' $bfile`
pz=`awk '{print $6;}' $bfile`

awk -f "$gdir"/tool.hgas/slice.awk py=$py dy=$dy $sfile \
    | awk '{$4-=px;$5-=py;$6-=pz; print $0;}' px=$px py=$py pz=$pz
