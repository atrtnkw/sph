tbgn=0
tend=300
idir=b1.10-1.00/s03.bsep
tcol=19

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    file="$idir"/sph_t"$time".dat
    if ! test -e $file
    then
        echo "Not found file \"$file\""
        exit
    fi
    sort -g -k $tcol $file \
        | tail -n1 \
        | awk '{printf("%4d %e\n", time, $(tcol));}' time=$time tcol=$tcol
done
