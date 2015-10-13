pbgn=0
pend=143
idir=b1.10-1.00/s03.bsep
time=0290

for proc in $(seq -f "%06g" $pbgn 1 $pend)
do
    ifile="$idir"/t"$time"_p"$proc".hexa
    ofile="$idir"/new_t"$time"_p"$proc".hexa
    awk '{if(NR==4){printf("40352c5e5bf00776\n");} print $0;}' $ifile > $ofile
done
