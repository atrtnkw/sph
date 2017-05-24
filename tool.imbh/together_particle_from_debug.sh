if test $# -ne 2
then
    echo "sh $0 <idir_list> <odir>"
    exit
fi

idir_list=$1
odir=$2

if ! test -e $odir
then
    mkdir $odir
fi

for idir in $idir_list
do
    ifile="$idir"/id_all.dat
    ilist=`cat $idir/id.list`
    for id in $ilist
    do
        pid=`printf "%010d" $id`
        prevt="-1"
        if test -e "$odir"/histroy.dat."$pid"
        then
            prevt=`tail -n1 "$odir"/histroy.dat."$pid" | awk '{print $1;}'`
        fi    
        awk '{if($1>prevt&&$2==id)print $1, $22, $17, $33, $34, $35;}' \
            prevt=$prevt id=$id $ifile >> "$odir"/histroy.dat."$pid"
    done
done

for idir in $idir_list
do
    ilist=`cat $idir/id.list`
    for id in $ilist
    do
        pid=`printf "%010d" $id`
        name=`head -n1 "$odir"/histroy.dat."$pid" | awk '{if($4>0.1){print "he"}else{print "co"}}'`
        awk '{print $1, $2, $3;}' "$odir"/histroy.dat."$pid" > "$odir"/"$name".dat."$pid"
    done
    break
done

