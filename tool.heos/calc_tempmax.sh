tbgn=355
tend=500

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    for sd in $(seq -f "%02g" 0 1 20)
    do
        file=t"$sd"/sph_t"$time".dat
#        file=snap/sph_t"$time".dat
        if test -e $file
        then
            ttemp=`expr $time + 0`
            printf "%d " $ttemp
            awk '{print $20, $17;}' $file | sort -g | tail -n1
            break
        fi
    done
done
