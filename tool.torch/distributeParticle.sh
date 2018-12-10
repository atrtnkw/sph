idir=t01
odir=init

if test -e init
then
    rm -r init
fi

mkdir -p init

awk -f ~/git-sph/tool.torch/distributeParticle.awk odir=$odir "$idir"/debug*.log 

files=`ls init/id*`

for file in $files
do
    sort -g $file > "$odir"/hoge
    mv "$odir"/hoge $file
done

