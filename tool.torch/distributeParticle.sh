if test $# -ne 2
then
    echo "sh $0 <idirs> <odir>" >&2
    exit
fi

idirs=$1
odir=$2
#idirs="t01 t99"
#odir=init.test

if test -e $odir
then
    echo "Error: Directory $odir exists" >&2
    exit
fi

mkdir -p $odir

for idir in $idirs
do
    awk -f ~/git-sph/tool.torch/distributeParticle.awk odir=$odir "$idir"/debug*.log 
done

files=`ls $odir/id*`

for file in $files
do
    sort -g $file > "$odir"/hoge
#    mv "$odir"/hoge $file
    awk 'BEGIN{tprev=-1;}{if($1<=tprev){next;}else{print $0;tprev=$1;}}' "$odir"/hoge > $file
done

