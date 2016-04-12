#ifile=data/r008k/cowd1.2_pw.bh5e2.b03.00/sph_t0128.dat
ifile=preliminary/sph_t0128.dat
tmin=2e9

echo "1" > dummy.hoge

################
# density bin
################
for idens in 1e6 2e6 5e6 1e7 2e7 5e7 1e8 2e8
do
    idensmin=`awk '{print idens * 0.9}' idens=$idens dummy.hoge`
    idensmax=`awk '{print idens * 1.1}' idens=$idens dummy.hoge`
    awk '{if($46>tmin&&(rmin<$47&&$47<rmax)) print $0;}' tmin=$tmin rmin=$idensmin rmax=$idensmax \
        $ifile > tmp.ifile
    nph=`wc -l tmp.ifile | awk '{print int($1/2.)}'`
    sort -g -k 46 tmp.ifile | tail -n$nph | head -n1
done

################
#temperature bin
################
itemp=8e9
itempmin=`awk '{print itemp * 0.99}' itemp=$itemp dummy.hoge`
itempmax=`awk '{print itemp * 1.01}' itemp=$itemp dummy.hoge`
awk '{if(tmin<$46&&$46<tmax) print $0;}' tmin=$itempmin tmax=$itempmax $ifile > tmp.ifile
nph=`wc -l tmp.ifile | awk '{print int($1*0.1)}'`
sort -g -k 47 tmp.ifile | tail -n$nph | head -n1

################
# density max
################
sort -g -k 47 $ifile | tail -n1

rm -f dummy.hoge tmp.ifile
