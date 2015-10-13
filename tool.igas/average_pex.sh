for time in $(seq -f "%04g" 0 1 55)
do
    ifile=pex.ns32.dt/sph_t"$time".dat
    ofile=pex.ns32.dt/profile_t"$time".dat
    awk -f awk/average_pex.awk width=0.01 $ifile > $ofile
done
