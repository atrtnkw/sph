tbgn=0
tend=500

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    file=unfy/sph_t"$time".dat
    if test -e $file
    then
        ttemp=`expr $time + 0`
        printf "%d " $ttemp
        awk '{print $20, $17;}' $file | sort -g | tail -n1
    fi
done
