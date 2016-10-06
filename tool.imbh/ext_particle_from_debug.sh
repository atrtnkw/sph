if test $# -ne 2
then
    echo "sh $0 <idir> <odir>"
    exit
fi

#idir=cowd1.2_pw.bh5e2.b03.00
idir=$1
odir=$2

idlist=`awk '{print $1}' $idir/id.list`
ofile="$odir"/id_all.dat
if test -e $ofile
then
    rm -f $ofile
fi
for id in $idlist
do
    pid=`printf "%09d" $id`
#    ofile="$idir"/id"$pid".dat
    olist="$idir"/id"$pid".list
#    awk '{if($2==id) print $0}' id=$id "$idir"/debug* | sort -g > $ofile
    awk '{if($2==id) print $0}' id=$id "$idir"/debug* | sort -g >> $ofile
    echo "" >> $ofile
#    awk '{print $1, $22, $17}' $ofile > $olist
done
