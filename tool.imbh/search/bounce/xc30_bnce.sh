#PBS -N bounce
##PBS -l mppwidth=128
#PBS -l mppwidth=16
#PBS -j oe
##PBS -q large-b+
##PBS -q bulk-b+
#PBS -q short-b+
##PBS -q debug

## The number of processes must be 2^n.
## The number of "nnxxx" must be 2^n.

#NPROCESS=128
NPROCESS=16
export OMP_NUM_THREADS=1

itbgn=61
itend=61
idtsp=2
idir=../r001m/run.cowd0.60_h100-100_bh3e2_b05.00
odir=../r001m/run.cowd0.60_h100-100_bh3e2_b05.00/fitting.mach04.+-
#odir=../r002m/run.cowd0.60_h100-100_bh3e2_b05.00/fitting.chkmach04_r001m
nfile=192
#nfile=288
#
fflag=1
xmin0=-5e9
ymin0=-5e9
xmax=1e11
width=10e9
nnxxx=512
#

cd $PBS_O_WORKDIR

if ! test -e $odir
then
    mkdir $odir
fi

cp $0 $odir/xc30_bnce.sh

for itime in $(seq -f "%04g" $itbgn $idtsp $itend)
do
    for dnum in $(seq -f "%02g" 0 1 99)
    do
        pnfile=`printf "%06d" $nfile`
        if ! test -e "$idir"/t"$dnum"/sph_t"$itime"_p"$pnfile"_i000000.dat
        then
            echo "Not found"$idir"/t"$dnum"/sph_t"$itime"_p"$pnfile".dat"
            continue
        fi
        
        itype="$idir"/t"$dnum"/sph_t"$itime"
        otype="$odir"/t"$itime"
        
        echo "$itype"               > input.list
        echo "$fflag $nfile"       >> input.list
        echo "$xmin0 $ymin0 $xmax" >> input.list
        echo "$width $nnxxx"       >> input.list
        echo "$otype"              >> input.list
        
        aprun -n $NPROCESS -d $OMP_NUM_THREADS ./run input.list
        
        break
    done
done
