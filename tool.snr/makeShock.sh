files=`ls mesh_*.dat`

for file in $files
do
    sort -g -k 4 $file | tail -n1
done
