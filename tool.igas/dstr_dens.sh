if test $# -ne 1
then
    echo "sh $0 <ifile>a"
    exit
fi

file=$1

awk 'BEGIN{xc=yc=zc=0.;n=0;}{xc+=$4;yc+=$5;zc+=$6;n++;}END{print xc/n, yc/n, zc/n}' $file > tmp.xc
xc=`awk '{print $1}' tmp.xc`
yc=`awk '{print $2}' tmp.xc`
zc=`awk '{print $3}' tmp.xc`

echo "# r, rho, temperature, kernel "
awk '{printf("%e %e %e %e\n", sqrt(($4-xc)**2+($5-yc)**2+($6-zc)**2), $16, $19, $15);}' \
    xc=$xc yc=$yc zc=$zc $file \
    | sort -g

rm -f tmp.*
