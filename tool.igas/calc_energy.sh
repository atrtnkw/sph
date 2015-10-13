if test $# -ne 1
then
    echo "sh $0 <ifile>a"
    exit
fi

file=$1

awk 'BEGIN{ek=ep=et=0.;n=0;}
{
ek += 0.5 * $3 * ($7**2 + $8**2 + $9**2);
ep += 0.5 * $3 * $25;
et +=       $3 * $13;
}
END{printf("%+e %+e %+e\n", ek, ep, et);}' $file
