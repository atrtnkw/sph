idir=init
odir=data

files=`ls "$idir"/id*.in`

if ! test -e $odir
then
    mkdir -p $odir
fi

for file in $files
do
    id=`echo $file | cut -c 8-15`
    echo -e ""   > input
    echo -e "6" >> input
    echo -e "9" >> input
    echo -e ""  >> input
    echo -e "$file" >> input
    echo -e "3" >> input
    echo -e "0. 0. 0. 0.5 0. 0. 0.5 0. 0. 0. 0.0 0. 0. 0. 0. 0. 0. 0." >> input
    echo -e "data/id$id.out_" >> input
    echo -e "-1" >> input
    ./run < input
done
