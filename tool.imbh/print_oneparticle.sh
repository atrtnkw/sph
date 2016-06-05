if test $# -ne 2
then
    echo "sh $0 <idir> <ldir>" 1>&2
    exit
fi

idir=$1
ldir=$2

lfile="$ldir"/idlist.dat
list=`awk '{print $1}' $lfile`

for id in $list
do
    awk '{if($2==id) print $0;}' id=$id "$idir"/debug* | sort -g
    printf "\n"
done


