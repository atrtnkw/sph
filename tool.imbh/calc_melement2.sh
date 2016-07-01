if test $# -ne 4
then
    echo "sh $0 <idir> <tbgn> <tend> <dtsp>" >&2
    exit
fi

idir=$1
tbgn=$2
tend=$3
dtsp=$4

echo "Directory $idir" >&2

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    awk -f ~/git-sph/tool.imbh/calc_melement.awk time=$time dtsp=$dtsp "$idir"/sph_t"$time"_*.dat
#    awk 'BEGIN{msun=1.989e33;m=0.;mf=0.;mn=0.;}
#{m+=$3*($43+$44);mf+=$3*$43;mn+=$3*$44}
#END{printf("%16.8f %+e %+e %+e\n", time/256., m/msun, mf/msun, mn/msun);}' \
#        time=$time "$idir"/sph_t"$time".dat
done
