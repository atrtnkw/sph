if test $# -ne 1
then
    echo "sh $0 <ifile>" >&2
    exit
fi

ifile=$1

nline=`wc -l $ifile | awk '{print $1}'`
remnd=`echo "$nline%800" | bc`
if ! test $remnd -eq 0
then
    echo "Error: Update is needed. See this source.\n" >&2
    exit
fi
pitch=`echo "$nline/800" | bc`

awk -f ~/git-sph/tool.imbh/reverseTime.awk pitch=$pitch $ifile
