if test $# -ne 6
then
    echo "sh $0 <idir> <x/y/z column> <x/y/z position> <tbgn> <tend> <dt>" >&2
    exit
fi

idir=$1
column=$2
dx=$3
tbgn=$4
tend=$5
dtsp=$6

if test $column = x
then
    col=1
elif  test $column = y
then
    col=2
elif  test $column = z
then
    col=3
else
    echo "Error: Not found such column $column" >&2
    exit
fi

for time in $(seq -f "%04g" $tbgn $dtsp $tend)
do
    awk '{if($(col)==dx) print $0}' col=$col dx=$dx "$idir"/mesh_t"$time".dat \
        >  "$idir"/mesh_t"$time".dat.d"$column"
done
