if test $# -ne 1
then
    echo "sh $0 <idir>"
    exit
fi

idir=$1

awk '{if($1=="###"&&NF>40){$1="";print $0;}}' "$idir"/t*/debug_p* \
    | sort -n \
    | awk 'BEGIN{pid=-1;}{if($1!=pid) print $0; pid=$1;}'
#    | sort -n -k 2 \
#    | awk 'BEGIN{pid=-1;}{if($2!=pid) print $0; pid=$2;}'
