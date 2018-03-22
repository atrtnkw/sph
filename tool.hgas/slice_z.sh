if test $# -ne 6
then
    echo "sh $0 <idir> <z> <dz> <tbgn> <tend> <dtsp>"
    exit
fi

idir=$1
z=$2
dz=$3
tbgn=$4
tend=$5
dtsp=$6

for time in $(seq -f "%04g" $tbgn $dtsp $tend)
do
    ifile="$idir"/sph_t"$time".dat
    ofile="$idir"/sph_t"$time".dat.dz
    awk -f ~/git-sph/tool.hgas/slice.awk z=$z dz=$dz $ifile > $ofile
done
