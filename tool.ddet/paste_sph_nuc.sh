if test $# -ne 2
then
    echo "sh $0 <fnuc> <fsph>"
    exit
fi

fnuc=$1
fsph=$2

sort -n -k 37 $fnuc > tmp.fnuc
awk '{if(NR==1){printf("\n");}print $0;}' $fsph > tmp.fsph
paste tmp.fnuc tmp.fsph

rm tmp.fnuc tmp.fsph
