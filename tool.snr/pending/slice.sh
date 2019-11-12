if test $# -ne 4
then
    echo "sh $0 <idir> <tbgn> <tend> <dt>" >&2
    exit
fi

idir=$1
tbgn=$2
tend=$3
dtsp=$4

for time in $(seq -f "%04g" $tbgn $dtsp $tend)
do
    dx=`head -n1 "$idir"/mesh_t"$time".dat | awk '{print $1}'`
    awk '{if($(col)==dx) print $0}' col=1 dx=$dx "$idir"/mesh_t"$time".dat \
        >  "$idir"/mesh_t"$time".datx
    awk '{if($(col)==dx) print $0}' col=2 dx=$dx "$idir"/mesh_t"$time".dat \
        >  "$idir"/mesh_t"$time".daty
    awk '{if($(col)==dx) print $0}' col=3 dx=$dx "$idir"/mesh_t"$time".dat \
        >  "$idir"/mesh_t"$time".datz
done
