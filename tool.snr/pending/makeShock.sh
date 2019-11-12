if test $# -ne 3
then
    echo "sh $0 <tbgn> <tend> <dtsp>" >&2
    exit
fi

tbgn=$1
tend=$2
dtsp=$3

for time in $(seq $tbgn $dtsp $tend)
do
    ptime=`printf "%04d" $time`
    ifile=mesh_t"$ptime".dat
    echo $ifile >&2
    sort -g -k 4 $ifile | tail -n1 | awk '{print time, $0;}' time=$time
done
