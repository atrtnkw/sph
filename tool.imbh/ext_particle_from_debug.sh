idir=cowd1.2_pw.bh5e2.b03.00

idlist=`awk '{print $1}' $idir/id.list`
for id in $idlist
do
    pid=`printf "%06d" $id`
    ofile="$idir"/id"$pid".dat
    olist="$idir"/id"$pid".list
    awk '{if($2==id) print $0}' id=$id "$idir"/debug* | sort -g > $ofile
    awk '{print $1, $22, $17}' $ofile > $olist
done
