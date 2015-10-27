for time in $(seq 0 1 245)
do
    ptime=`printf "%04d" $time`
    ifile=unfy/sph_t"$ptime".dat
    printf "time: %4d\n" $time
    awk -f ref/print_smalltau.awk $ifile
    printf "\n"
done
