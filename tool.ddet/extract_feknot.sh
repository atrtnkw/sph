if test $# -ne 2
then
    echo "sh $0 <ifile> <otype>"
    exit
fi

ifile=$1
otype=$2

awk '{if(NR>1 && $25/$27<0.023 && $29/$27<0.029 && $27 != 0.){print $0;}}' $ifile > $otype.feknot
awk '{if(NR>1 && (!($25/$27<0.023 && $29/$27<0.029 && $27 != 0.))){print $0;}}' $ifile \
    > $otype.other
