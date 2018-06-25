if test $# -ne 1
then
    echo "sh $0 <ifile>" >&2 
    exit
fi

ifile=$1

awk '{if($1=="time:") print $0}' $ifile \
    | awk 'BEGIN{tsum=0.}{tsum+=$10;nstep=$2}END{printf("%+e\n", tsum/nstep);}'
