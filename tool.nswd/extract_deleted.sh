idir=.

awk '{if($1=="###"&&NF>40) print $0;}' "$idir"/t*/debug_p* \
    | sort -n -k 2 \
    | awk 'BEGIN{pid=-1;}{if($2!=pid) print $0; pid=$2;}'
