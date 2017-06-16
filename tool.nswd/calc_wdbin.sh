if test $# -ne 1
then
    echo "sh $0 <time>"
    exit
fi

idir=unfy
time=`printf "%04d" $1`

nfile="$idir"/bhns_t"$time".dat
sfile="$idir"/sph_t"$time".dat

mns=`awk '{print $3;}' $nfile`
xns=`awk '{print $4;}' $nfile`
yns=`awk '{print $5;}' $nfile`
zns=`awk '{print $6;}' $nfile`

awk 'BEGIN{ek=ep=ei=0.;g=6.67259e-8;ksrh=1.936492;}
{
ek += 0.5*$3*($7**2+$8**2+$9**2);
r2  = ($4-xns)**2+($5-yns)**2+($6-zns)**2;
eni = -0.5*$3*g*mns/sqrt(r2)
eei = -g*$3*$3/($17/ksrh);
ep += 0.5*$3*$25 + eni - eei;
en += eni;
ei += $3*$13;
}
END{printf("eb_wd-ns: %+e eb_wd-only: %+e\n", ek+ep+ei, (ep-2*en)+ei);}' \
mns=$mns xns=$xns yns=$yns zns=$zns $sfile
