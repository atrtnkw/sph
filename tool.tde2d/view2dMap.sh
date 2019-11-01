if test $# -ne 1
then
    echo "sh $0 <mesh_txxxx.dat>" >&2
    exit
fi

file=$1

gnuplot<<EOF
file="$file"
PointSize=0.2

set lmargin 7
set rmargin 2
set tmargin 2
set bmargin 3

set size 1
set size ratio 0.2

set multiplot

dy=0.23
ndy=-1.5

sinv=1e-8

unset key

#====================================

set origin 0, dy*(3+ndy)

set xrange [-2.5:2.5]
set xtic -3, 1, 3
set mxtic 4
set format x ""

set yrange [0:1]
set ytic -3, 0.5, 3
set mytic 5

set cbrange [0:8]
set cbtic 0, 2, 8
set mcbtic 1
set cblabel "dens" offset 1,0

pl \
file u (\$1*sinv):(\$2*sinv):(log10(\$5)) pal ps PointSize pt 5

#====================================

set origin 0, dy*(2+ndy)

set cbrange [0:1]
set cbtic 0, 1, 1
set mcbtic 1
set cblabel "56Ni"

pl \
file u (\$1*sinv):(\$2*sinv):19 pal ps PointSize pt 5

#====================================

set origin 0, dy*(1+ndy)

set cbrange [-0.5:0.5]
set cbtic -0.5, 0.5, 0.5
set mcbtic 1
set cblabel "vels"

pl \
file u (\$1*sinv):(\$2*sinv):(\$3*1e-9) pal ps PointSize pt 5

#====================================

set origin 0, dy*(0+ndy)

set format x "%g"

set cblabel "velz"

pl \
file u (\$1*sinv):(\$2*sinv):(\$4*1e-9) pal ps PointSize pt 5


#====================================

#pause -1
pause 100000
EOF
