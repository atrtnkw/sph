#PBS -N pxy
#PBS -l mppwidth=24
#PBS -j oe
##PBS -q large-a
#PBS -q bulk-a
##PBS -q short-a
##PBS -q debug

cd $PBS_O_WORKDIR

mtot=1.0
nres=64
nxnx=512
xmin=-1e10
xmax=+1e10
idir=../../nswd/data/bns-wd1.0/r064k/run/unfy
tbgn=446
tend=500
dtim=1
##########
odir=snap

###############################
nproc=1
export OMP_NUM_THREADS=24
nptcl=`echo "$mtot * 10 * $nres * 1024" | bc`

if ! test -e $odir
then
    mkdir $odir
fi

for time in $(seq -f "%04g" $tbgn $dtim $tend)
do
    ifile="$idir"/sph_t"$time".dat
    bfile="$idir"/bhns_t"$time".dat
    ofile="$odir"/pxy_t"$time".dat
    xpos=`awk '{print $4}' $bfile`
    ypos=`awk '{print $5}' $bfile`
    zpos=`awk '{print $6}' $bfile`
    echo "$nptcl $xpos $ypos $zpos $nxnx $xmin $xmax" \
        | awk '{printf("%d %+e %+e %+e %d %+e %+e\n", $1, $2, $3, $4, $5, $6, $7);}' \
        > header.tmp
    aprun -n $nproc -d $OMP_NUM_THREADS ./run header.tmp $ifile $ofile
done
