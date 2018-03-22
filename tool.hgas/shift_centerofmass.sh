if test $# -ne 4
then
    echo "sh $0 <idir> <tbgn> <tend> <dtsp>" 1>&2
    exit
fi

idir=$1
tbgn=$2
tend=$3
dtsp=$4

for time in $(seq -f "%04g" $tbgn $dtsp $tend)
do
    ifile="$idir"/sph_t"$time".dat
    ofile="$idir"/sph_t"$time".dat.sft
    if test -e $ifile
    then
        echo "File $ifile is found." 1>&2
    else
        echo "error: $ifile is not found." 1>&2
        continue
    fi
    awk -f ~/git-sph/tool.hgas/print_centerofmass.awk $ifile > tmp.com
    px=`awk '{if(NR==1) print $2}' tmp.com`
    py=`awk '{if(NR==1) print $3}' tmp.com`
    pz=`awk '{if(NR==1) print $4}' tmp.com`
    vx=`awk '{if(NR==2) print $2}' tmp.com`
    vy=`awk '{if(NR==2) print $3}' tmp.com`
    vz=`awk '{if(NR==2) print $4}' tmp.com`
    
    awk '{$4-=px;$5-=py;$6-=pz;$7-=vx;$8-=vy;$9-=vz; print $0;}' px=$px py=$py pz=$pz \
        vx=$vx vy=$vy vz=$vz $ifile > $ofile
done

rm -f tmp.com
