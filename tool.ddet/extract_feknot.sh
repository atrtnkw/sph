if test $# -ne 2
then
    echo "sh $0 <ifile> <otype>"
    exit
fi

ifile=$1
otype=$2

awk '{if(NR>1 && $27 > 0.85 && (5.3e9 < $84 && $84 < 5.7e9) && $25/$27< 0.023 && $29/$27< 0.029){print $0;}}' $ifile \
    > $otype.feknot
awk '{if(NR>1 && $27 > 0.85 && (5.3e9 < $84 && $84 < 5.7e9) && $25/$27< 0.023 && $29/$27>=0.029){print $0;}}' $ifile \
    > $otype.onlycr
awk '{if(NR>1 && $27 > 0.85 && (5.3e9 < $84 && $84 < 5.7e9) && $25/$27>=0.023 && $29/$27< 0.029){print $0;}}' $ifile \
    > $otype.onlyni
awk '{if(NR>1 && ($27 <= 0.85 || ($84 <= 5.3e9 || 5.7e9 <= $84) || ($25/$27>=0.023 && $29/$27>=0.029))){print $0;}}' $ifile \
    > $otype.other
