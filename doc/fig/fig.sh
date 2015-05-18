gnuplot $1.gpl
awk -f $1.awk tmp.eps > tmp2.eps
mv tmp2.eps $1.eps
rm tmp.eps

convert -density 288x288 $1.eps $1.png
extractbb $1.png
