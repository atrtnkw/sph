if test $# -ne 3
then
    echo "sh $0 <idir> <tbgn> <tend>" >&2
    exit
fi

idir=$1
tbgn=$2
tend=$3

tdir="$HOME"/git-sph

for time in $(seq $tbgn 1 $tend)
do
    ptim=`printf "%04d" $time`
    sfile="$idir"/sph_t"$ptim".dat
    bfile="$idir"/bhns_t"$ptim".dat

    cat $sfile $bfile \
        | awk 'BEGIN{cx=cy=cz=0.;vx=vy=vz=0.;ax=ay=az=0.;m=0.;}{cx+=$3*$4;cy+=$3*$5;cz+=$3*$6;vx+=$3*$7;vy+=$3*$8;vz+=$3*$9;ax+=$3*$10;ay+=$3*$11;az+=$3*$12;m+=$3;}END{print cx/m, cy/m, cz/m, vx/m, vy/m, vz/m, ax/m, ay/m, az/m;}' \
        > tmp.cm
    cx=`awk '{print $1;}' tmp.cm`
    cy=`awk '{print $2;}' tmp.cm`
    cz=`awk '{print $3;}' tmp.cm`
    vx=`awk '{print $4;}' tmp.cm`
    vy=`awk '{print $5;}' tmp.cm`
    vz=`awk '{print $6;}' tmp.cm`
    ax=`awk '{print $7;}' tmp.cm`
    ay=`awk '{print $8;}' tmp.cm`
    az=`awk '{print $9;}' tmp.cm`
    awk -f "$tdir"/tool.nswd/calc_d2quad.awk time=$time cx=$cx cy=$cy cz=$cz \
        vx=$vx vy=$vy vz=$vz ax=$ax ay=$ay az=$az $sfile $bfile
done

rm -f tmp.*
