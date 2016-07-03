if test $# -ne 1
then
    echo "sh $0 <idir>"
    exit
fi

idir=$1

for ndir in $(seq -f "%02g" 0 1 99)
do
    fdir="$idir"/t"$ndir"
    if ! test -e $fdir
    then
	break
    fi
    file="$fdir"/time.log
    if test -e "$idir"/unfy/time.log
    then
	step=`tail -n1 "$idir"/unfy/time.log | awk '{if($1=="time:") print $2}'`
    else
	step=-1
    fi
#    echo $file 1>&2
    awk '{if($2>step) print $0;}' step=$step $file
done
