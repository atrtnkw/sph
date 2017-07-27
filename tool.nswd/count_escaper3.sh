if test $# -ne 3
then
    echo "sh $0 <idir> <tbgn> <tend>" >&2
    exit
fi

idir=$1
tbgn=$2
tend=$3

echo "Directory $idir" >&2

for time in $(seq $tbgn 1 $tend)
do
    ptim=`printf "%04d" $time`
    if ! test -e "$idir"/sph_t"$ptim".dat -o -e "$idir"/sph_t"$ptim"*000000.dat
    then
	    continue
    fi
    cat "$idir"/sph_t"$ptim"*.dat |  awk 'BEGIN{n=0;ek=0.;}{e=0.5*($7**2+$8**2+$9**2)+$25+$13;if(e>0.){n++;ek+=e;}}END{printf("%4d %10d %+e\n", time, n, ek);}' time=$time
done
