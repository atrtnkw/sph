if test $# -ne 1
then
    echo "sh $0 <mesh_txxxx.dat>" >&2
    exit
fi

file=$1

gnuplot<<EOF
file="$file"

set lmargin 7
set rmargin 2
set tmargin 2
set bmargin 3

set size 0.9
set size ratio 0.3

set multiplot

dy=0.30

sinv=1e-8

#====================================

set origin 0, dy*(+1.2)

set xrange [-2.5:2.5]
set xtic -3, 1, 3
set mxtic 4
set format x ""

set yrange [0:1.5]
set ytic -3, 0.5, 3
set mytic 5

set cbrange [0:8]
set cbtic 0, 2, 8
set mcbtic 1

pl \
file u (\$1*sinv):(\$2*sinv):(log10(\$3)) pal ps 0.6 pt 5

#====================================

set origin 0, dy*(+0.2)

set cbrange [6:10]
set cbtic 0, 1, 10
set mcbtic 1

pl \
file u (\$1*sinv):(\$2*sinv):(log10(\$5)) pal ps 0.6 pt 5

#====================================

set origin 0, dy*(-0.8)

set format x "%g"

set cbrange [-1:3]
set cbtic -10, 1, 10
set mcbtic 1

pl \
file u (\$1*sinv):(\$2*sinv):(\$4*1e-9) pal ps 0.6 pt 5


#====================================

#pause -1
pause 100000
EOF
