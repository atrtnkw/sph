bfile=`ls -r ../t*/bhns_t*.dat`
sfile=`ls -r ../t*/sph_t*.dat`

for file in $bfile $sfile
do
    ln -s $file .
done
