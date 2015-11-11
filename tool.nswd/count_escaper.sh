tbgn=0
tend=500
idir=run

echo "Directory $idir" >&2

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    awk 'BEGIN{n=0;}{e=0.5*($7**2+$8**2+$9**2)+$25;if(e>0.)n++;}END{print time, n;}' \
        time=$time "$idir"/sph_t"$time".dat
done
