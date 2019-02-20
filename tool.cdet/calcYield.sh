if test $# -ne 3
then
    echo "sh $0 <idir> <nfile> <0:final,1:decayed>" >&2
    echo "Caution: Limited to 1Msun WD!" >&2
    exit
fi

idir=$1
nfile=$2
type=$3

if test $type -eq 0
then
    ftype=final
elif test $type -eq 1
then
    ftype=decayed
else
    echo "Error: argv[2] is only 0 or 1." >&2
    exit
fi

files=`ls "$idir"/id*.out_0000_z1.dat`

if test -e tmp.dat
then
    rm -f tmp.dat
fi

nline=-1

for n in $(seq 1 1 $nfile)
do
    pn=`printf "%010d" $n`
    dfile="$idir"/id"$pn".out_0000_z1.dat
    if ! test -e $dfile
    then
	echo "Error: Not fould $dfile" >&2
	exit
    fi

    dmass=`awk 'BEGIN{dR=6e6;msun=1.989e33}{if($1==1){printf("%+e\n", $4*4./3.*3.14*((dR*(n-0.5))**3-(dR*((n>1?n:1.5)-1.5))**3)/msun);}}' n=$n $dfile`

    fresult="$idir"/id"$pn".out_"$ftype".dat
    if test $nline -eq -1
    then
	nline=`wc -l $fresult | awk '{print $1}'`
    fi
    awk '{$3*=dmass; print $0;}' dmass=$dmass $fresult >> tmp.dat
done

awk -f $HOME/git-torch/tool/reduceTorchResult.awk nline=$nline type=$type tmp.dat \
    | awk '{printf("%10s %+e\n", $1, $2*nfile);}' nfile=$nfile

rm -f tmp.dat
