if test $# -ne 3
then
    echo "sh $0 <ifile> <ofile> <starID<0/1>" >&2
    exit
fi

ifile=$1
ofile=$2
sid=$3

awk 'BEGIN{n=cx=cy=cz=0.}{if($2==sid){n++;cx+=$4;cy+=$5;cz+=$6;}}END{printf("%+.16e %+.16e %+.16e\n", cx/n,cy/n,cz/n)}' \
    sid=$sid $ifile > tmp.cntr

cx=`awk '{printf("%+.16e\n", $1)}' tmp.cntr`
cy=`awk '{printf("%+.16e\n", $2)}' tmp.cntr`
cz=`awk '{printf("%+.16e\n", $3)}' tmp.cntr`

printf "# %+.16e %+.16e %+.16e\n" $cx $cy $cz >&2
printf "# %+.16e %+.16e %+.16e\n" $cx $cy $cz > $ofile

awk '{if($2==sid){printf("%+.16e %10d\n", sqrt(($4-cx)**2+($5-cy)**2+($6-cz)**2), $1)}}' \
    sid=$sid cx=$cx cy=$cy cz=$cz $ifile \
    | sort -g >> $ofile

rm -f tmp.*
