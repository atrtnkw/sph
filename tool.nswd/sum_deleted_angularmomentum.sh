if test $# -ne 4
then
    echo "sh $0 <idir> <tdirbgn(txx)> <tdirend(tyy)> <tfin>"
    exit
fi

idir=$1
tbgn=$2
tend=$3
tfin=$4

for time in 0 $tfin
do
    sh ~/git-sph/tool.nswd/shift_bhnscenter.sh "$idir"/unfy $time \
	| awk 'BEGIN{angx=angy=angz=0.;}{ms=$3;px=$4;py=$5;pz=$6;vx=$7;vy=$8;vz=$9;
angx+=ms*(py*vz-pz*vy);angy+=ms*(pz*vx-px*vz);angz+=ms*(px*vy-py*vx)}END{print angx, angy, angz;}' \
    > tmp.ang
    if test $time -eq 0
    then
	angx0=`awk '{printf("%+e\n", $1);}' tmp.ang`
	angy0=`awk '{printf("%+e\n", $2);}' tmp.ang`
	angz0=`awk '{printf("%+e\n", $3);}' tmp.ang`
    else
	angx1=`awk '{printf("%+e\n", $1);}' tmp.ang`
	angy1=`awk '{printf("%+e\n", $2);}' tmp.ang`
	angz1=`awk '{printf("%+e\n", $3);}' tmp.ang`
    fi
done

if test -e tmp.sum
then
    rm -f tmp.sum
fi

tprev=0.
pangx=0.
pangy=0.
pangz=0.
for ii in $(seq -f "%02g" $tbgn $tend)
do
    awk -f ~/git-sph/tool.nswd/sum_deleted_angularmomentum.awk "$idir"/t"$ii"/debug_p* \
	| sort -g \
	| awk 'BEGIN{angx=angy=angz=0.;}{angx+=$3;angy+=$4;angz+=$5;printf("%18.10f %+e %+e %+e\n", $1, angx, angy, angz);}' tp=$tprev angx=$pangx angy=$pangy angz=$pangz >> tmp.sum
    tprev=`tail -n1 tmp.sum | awk '{print $1}'`
    pangx=`tail -n1 tmp.sum | awk '{print $2}'`
    pangy=`tail -n1 tmp.sum | awk '{print $3}'`
    pangz=`tail -n1 tmp.sum | awk '{print $4}'`
done

tail -n1 tmp.sum \
    | awk '{printf("# error: angx=%+.3e angy=%+.3e angz=%+.3e\n", ($2+x1-x0)/x0, ($3+y1-y0)/y0, ($4+z1-z0)/z0);}' x0=$angx0 y0=$angy0 z0=$angz0 x1=$angx1 y1=$angy1 z1=$angz1

awk '{printf("%16.10f %+e %+e %+e %+e %+e %+e\n", $1, $2, $3, $4, angx0, angy0, angz0);}' \
    angx0=$angx0 angy0=$angy0 angz0=$angz0 tmp.sum

rm tmp.ang tmp.sum

