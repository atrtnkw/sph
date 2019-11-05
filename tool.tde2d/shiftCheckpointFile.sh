if test $# -ne 3
then
    echo "sh $0 <nshift> <ibgn> <iend>" >&2
    exit
fi

nshift=$1
ibgn=$2
iend=$3

for torg in $(seq $ibgn 1 $iend)
do
    tmov=`echo "$torg + $nshift" | bc`
    ptorg=`printf "%04d" $torg`
    ptmov=`printf "%04d" $tmov`
    ifile=tde2d_hdf5_chk_$ptorg
    ofile=tde2d_hdf5_chk_$ptmov
    if test -e $ofile
    then
        echo "Error: $ofile exists." >&2
        exit
    fi
    mv $ifile $ofile
done
