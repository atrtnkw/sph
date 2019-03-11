if test $# -ne 5
then
    echo "sh $0 <bdir> <tbgn> <tend> <dtsp> <starid[0/1]>" >&2
    exit
fi

bdir=$1
tbgn=$2
tend=$3
dtsp=$4
stid=$5
mlag="1 5 10 20 30 40 50 60 70 80 90 95 99"

for time in $(seq $tbgn $dtsp $tend)
do
    ptim=`printf "%04d" $time`
    for tdir in $(seq -f "%02g" 0 1 99)
    do
        ifile="$bdir"/t"$tdir"/sph_t"$ptim"*000000.dat
        if test -e $ifile
        then
            break
        fi
    done
    idir="$bdir"/t"$tdir"
    awk 'BEGIN{dx=dy=dz=dd=0.;}{if($2==0){dx+=$16*$4;dy+=$16*$5;dz+=$16*$6;dd+=$16}}END{print dx/dd, dy/dd, dz/dd}' "$idir"/sph_t"$ptim"*.dat > tmp.dat
    xc=`awk '{print $1}' tmp.dat`
    yc=`awk '{print $2}' tmp.dat`
    zc=`awk '{print $3}' tmp.dat`

    nptcl=`awk '{if($2==stid)print $0;}' stid=$stid "$idir"/sph_t"$ptim"*.dat | wc -l | awk '{print $1}'`
    awk '{if($2==stid){print ($4-xc)**2+($5-yc)**2+($6-zc)**2;}}' stid=$stid xc=$xc yc=$yc zc=$zc "$idir"/sph_t"$ptim"*.dat \
        | sort -g \
        > tmp.r2

    printf "%5d" $time
    for frac in $mlag
    do
        iptcl=`echo "scale=0; $nptcl * $frac * 0.01" | bc | awk '{print int($1);}'`
        awk '{if(NR==iptcl) printf(" %+e", sqrt($1))}' iptcl=$iptcl tmp.r2 
    done
    printf "\n"

done

rm -f tmp.*
