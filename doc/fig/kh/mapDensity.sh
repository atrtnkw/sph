for time in 0000 0001 0002 0003
do
awk -f mapDensity.awk xmin=-0.5 xmax=0.5 ymin=-0.5 ymax=0.5 width=0.00390625 sph_t"$time".dat > sph_t"$time".map
done
