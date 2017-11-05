if test $# -ne 1
then
    echo "sh $0 <ifile>" >&2
    exit
fi

ifile=$1

mtot=`awk 'BEGIN{m=0}{if($1>=0)m+=$2;}END{print m}' $ifile`
awk 'BEGIN{m=0}{if($1>=0){m+=$2;print $1, $2, m/mtot}}' mtot=$mtot $ifile
