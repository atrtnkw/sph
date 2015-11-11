tbgn=0
tend=500
idir=run

echo "Directory $idir" >&2

for time in $(seq -f "%04g" $tbgn 1 $tend)
do
    wc -l "$idir"/sph_t"$time".dat | awk '{print time, $1}' time=$time
done
