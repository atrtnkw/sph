if test $# -ne 2
then
    echo "sh $0 <ifile> <n/1msun>" >&2
    exit
fi

ifile=$1
nmsun=$2

echo "# ifile=$ifile n/1msun=$nmsun"
awk '{if(NR==1){pn=$2;}else{printf("%4d %+e\n", $1, (pn-$2)/nsun);} pn=$2;}' nsun=$nmsun $ifile
