if test $# -ne 3
then
    echo "sh $0 <ifile> <ofile> <starID<0/1>" >&2
    exit
fi

ifile=$1
ofile=$2
sid=$3

awk 'BEGIN{n=cx=cy=cz=0.}{if($2==sid){n++;cx+=$4;cy+=$5;cz+=$6;}}END{print cx/n,cy/n,cz/n}' \
    sid=$sid $ifile > tmp.cntr

cx=`awk '{print $1}' tmp.cntr`
cy=`awk '{print $2}' tmp.cntr`
cz=`awk '{print $3}' tmp.cntr`

echo "# $cx $cy $cz" >&2
echo "# $cx $cy $cz" > $ofile

awk '{if($2==sid){printf("%+e %10d\n", sqrt(($4-cx)**2+($5-cy)**2+($6-cz)**2), $1)}}' \
    sid=$sid cx=$cx cy=$cy cz=$cz $ifile \
    | sort -g >> $ofile

rm -f tmp.*
