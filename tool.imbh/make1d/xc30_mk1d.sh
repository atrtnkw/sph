#PBS -N mk1d
#PBS -l mppwidth=1
#PBS -j oe
##PBS -q short-a
#PBS -q bulk-a

NPARALLEL=1
NPROCESS=$NPARALLEL
fexe=run
file=input.list

flag=0
nhlf=10000
frvz=2
fhe=0.00
fcb=0.00
fox=0.60
fne=0.35
fmg=0.05
#fcb=0.50
#fox=0.50
#fne=0.00
#fmg=0.00
if test $flag = 0
then
    echo "lowT"
    if test $frvz = 1
    then
        vz_c=+1.6e+1    
        vzvz=1.0
    else
        vz_c=+2.4e+1
        vzvz=1.5
    fi
    dn_a=+2.0e-8
    dn_b=+1.2e+8
    zmax=+7.0e+7
    temp=lowT
elif test $flag = 1
then
    echo "hghT"
    if test $frvz = 1
    then
        vz_c=+3.2e+1
        vzvz=1.0
    else
        vz_c=+4.8e+1
        vzvz=1.5
    fi
    dn_a=+3.0e-8
    dn_b=+9.5e+7
    zmax=+6.0e+7
    temp=hghT
elif test $flag = 2
then
    echo "extT"
    if test $frvz = 1
    then
        vz_c=+4.0e+1
        vzvz=1.0
    else
        vz_c=+6.0e+1
        vzvz=1.5
    fi
    dn_a=+3.5e-8
    dn_b=+6.5e+7
    zmax=+4.0e+7
    temp=extT
elif test $flag = 3
then
    echo "hypT"
    if test $frvz = 1
    then
        vz_c=+8.0e+1
        vzvz=1.0
    else
        vz_c=+12.0e+1
        vzvz=1.5
    fi
    dn_a=+2.5e-8
    dn_b=+1.3e+7
    zmax=+2.5e+7
    temp=hypT
else
    echo "ultT"
    if test $frvz = 1
    then
        vz_c=+7.5e+1
        vzvz=1.0
    else
        vz_c=+11.25e+1
        vzvz=1.5
    fi
    dn_a=+2.0e-8
    dn_b=+1.3e+7
    zmax=+2.5e+7
    temp=ultT
fi
if test $nhlf = 10
then
    ndir=n1e1
elif test $nhlf = 100
then
    ndir=n1e2
elif test $nhlf = 1000
then
    ndir=n1e3
else
    ndir=n1e4
fi
#
filetype=../sim1d/"$ndir"/init/"$temp"_vz"$vzvz"
#filetype="$temp"_vz"$vzvz"
#filetype=../sim1d/"$ndir"/init/co_"$temp"_vz"$vzvz"
temp=+1.0e+7
leng=+1.0e+10

cd $PBS_O_WORKDIR

echo "$dn_a $dn_b $vz_c" > $file
echo "$fhe $fcb $fox $fne $fmg" >> $file
echo "$nhlf $temp $zmax $leng" >> $file
echo "$filetype" >>$file

export OMP_NUM_THREADS=`echo "$NPARALLEL / $NPROCESS" | bc`

aprun -n $NPROCESS -d $OMP_NUM_THREADS ./"$fexe" $file
