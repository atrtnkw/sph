tbgn=0
tend=500
idir=run
vs=1e9

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    awk 'BEGIN{n=0;}{v=sqrt($7**2+$8**2+$9**2);if(v>vs)n++;}END{print time, n;}' \
        time=$time vs=$vs "$idir"/sph_t"$time".dat
done
