if test $# -ne 1
then
    echo "sh $0 <ifile>" >&2
    exit
fi

ifile=$1

nline=`wc -l $ifile | awk '{print $1}'`
radius=`tail -n1 $ifile | awk '{print $1}'`
tmass=`tail -n1 $ifile | awk '{print $6}'`

echo "$nline"
echo "1.4000E+00"
echo "$radius"
cat $ifile
echo $tmass
