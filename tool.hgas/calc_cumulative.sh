if test $# -ne 2
then
    echo "sh $0 <ftype> <column>" >&2
    exit
fi

ftype=$1
col=$2

n=`cat "$ftype"*.dat | wc -l`

awk '{print $(col)}' col=$col "$ftype"*.dat \
    | sort -g \
    | awk '{print $1, NR/n}' n=$n
