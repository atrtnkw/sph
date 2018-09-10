if test $# -ne 2
then
    echo "sh $0 <xxx.awk> <NDIM(1/3)>" >&2
    exit
fi

fawk=$1
ndim=$2

touch dummy

if test $ndim -eq 1
then
    awk -f $fawk xmax=2e19 dummy
elif test $ndim -eq 3
then
    awk -f $fawk xmax=2e19 ymax=2e19 zmax=2e19 dummy
else
    echo "Error: NDIM should be 1 or 3" >&2
    exit
fi

rm -f dummy
