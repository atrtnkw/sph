tbgn=0
tend=500
idir=r002k/run

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    bhfile="$idir"/bhns_t"$time".dat
    spfile="$idir"/sph_t"$time".dat

    x=`awk '{print $4}' $bhfile`
    y=`awk '{print $5}' $bhfile`
    z=`awk '{print $5}' $bhfile`

    awk 'BEGIN{r2min=1e30}{r2=($4-x)**2+($5-y)**2+($6-z)**2;if(r2<r2min){r2min=r2;}}END{print time, sqrt(r2min);}' x=$x y=$y z=$z time=$time $spfile

done
