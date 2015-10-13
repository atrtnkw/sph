if test $# -ne 2
then
    echo "sh $0 <????.data> <nfile>"
    exit
fi

tfile=$1
nfile=$2

ifile=$1.data
nptcl=`wc -l $ifile | awk '{print $1}'`

awk 'BEGIN{
i = 0;
}
{
if(NR == int(nptcl * i / nfile) + 1) {
filename = sprintf("%s_p%06d_i%06d.data", tfile, nfile, i);
i++;
}
print $0 > filename;
}' tfile=$tfile nfile=$nfile nptcl=$nptcl $ifile
