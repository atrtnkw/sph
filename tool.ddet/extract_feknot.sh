if test $# -ne 2
then
    echo "sh $0 <ifile> <otype>"
    exit
fi

ifile=$1
otype=$2

awk '{if(NR>1 && $27 > 0.85 && (5.3e9 < $84 && $84 < 5.7e9) && $25/$27< 0.023 && $29/$27< 0.029){print $0;}}' $ifile \
    > $otype.fe_ytp
awk '{if(NR>1 && $27 > 0.85 && (5.3e9 < $84 && $84 < 5.7e9) && $25/$27< 0.023 && $29/$27>=0.029){print $0;}}' $ifile \
    > $otype.cr_ytp
awk '{if(NR>1 && $27 > 0.85 && (5.3e9 < $84 && $84 < 5.7e9) && $25/$27>=0.023 && $29/$27< 0.029){print $0;}}' $ifile \
    > $otype.ni_ytp
awk '{if(NR>1 && $27 > 0.85 && (5.3e9 < $84 && $84 < 5.7e9) && $25/$27>=0.023 && $29/$27>=0.029){print $0;}}' $ifile \
    > $otype.ng_ytp

awk '{if(NR>1 && $27 > 0.85 && ($84 <= 5.3e9 || 5.7e9 <= $84) && $25/$27< 0.023 && $29/$27< 0.029){print $0;}}' $ifile \
    > $otype.fe_ntp
awk '{if(NR>1 && $27 > 0.85 && ($84 <= 5.3e9 || 5.7e9 <= $84) && $25/$27< 0.023 && $29/$27>=0.029){print $0;}}' $ifile \
    > $otype.cr_ntp
awk '{if(NR>1 && $27 > 0.85 && ($84 <= 5.3e9 || 5.7e9 <= $84) && $25/$27>=0.023 && $29/$27< 0.029){print $0;}}' $ifile \
    > $otype.ni_ntp
awk '{if(NR>1 && $27 > 0.85 && ($84 <= 5.3e9 || 5.7e9 <= $84) && $25/$27>=0.023 && $29/$27>=0.029){print $0;}}' $ifile \
    > $otype.ng_ntp

awk '{if(NR>1 && $27 <= 0.85) print $0;}' $ifile > $otype.other
