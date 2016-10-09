#PBS -N mk1d
#PBS -l mppwidth=1
#PBS -j oe
#PBS -q short-a

NPARALLEL=1
NPROCESS=$NPARALLEL
fexe=run
file=input.list

dn_a=+2.0e-8
dn_b=+1.2e+8
vz_c=+1.6e+1
fhe=0.00
fcb=0.00
fox=0.60
fne=0.35
fmg=0.05
nhlf=1000
temp=+1.0e+7
zmax=+8.0e+7
leng=+4.0e+8
filetype=../sim1d/init/lowT_vz1

cd $PBS_O_WORKDIR

echo "$dn_a $dn_b $vz_c" > $file
echo "$fhe $fcb $fox $fne $fmg" >> $file
echo "$nhlf $temp $zmax $leng" >> $file
echo "$filetype" >>$file

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./"$fexe" $file
