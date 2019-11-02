if test $# -ne 3
then
    echo "sh $0 <xmin> <xlen> <mesh_txxxx.dat>" >&2
    exit
fi

xmin=$1
xlen=$2
file=$3

gnuplot<<EOF
file = "$file"
xmin = $xmin
xlen = $xlen
xmax = xmin + xlen
PointSize=0.25
Pitch=1

set lmargin 7
set rmargin 2
set tmargin 2
set bmargin 3

set size 1
set size ratio 10*1e7/xlen

set multiplot

px=-0.45
dx=0.07*xlen/1e7
dy=0.23
ndy=-1.5

sinv=1e-8

unset key

#====================================

set origin px+dx*0, dy*0

set title "dens" offset 0, -1

set xrange [xmin*sinv:xmax*sinv]
set xtic xmin*sinv-xlen*sinv*0.5, xlen*sinv, xmax*sinv+xlen*sinv*0.5
set mxtic 1
set format x "%g"

set yrange [0:1]
set ytic -3, 0.5, 3
set mytic 5

unset colorbox
set cbrange [0:8]

pl \
file u (\$1*sinv):(\$2*sinv):(log10(\$5)) ev Pitch pal ps PointSize pt 5

#====================================

set origin px+dx*1, dy*0

set title "temp"

set format y ""

set cbrange [6:10]

pl \
file u (\$1*sinv):(\$2*sinv):(log10(\$6)) ev Pitch pal ps PointSize pt 5

#====================================

set origin px+dx*2, dy*0

set title "vels"

set cbrange [-0.5:0.5]

pl \
file u (\$1*sinv):(\$2*sinv):(\$3*1e-9) ev Pitch pal ps PointSize pt 5

#====================================

set origin px+dx*3, dy*0

set title "velz"

set cbrange [-1:3]

pl \
file u (\$1*sinv):(\$2*sinv):(\$4*1e-9) ev Pitch pal ps PointSize pt 5

#====================================

set origin px+dx*4, dy*0

set title "4He"

set format y ""

set cbrange [0:1]
set cbtic 0, 1, 1
set mcbtic 1

pl \
file u (\$1*sinv):(\$2*sinv):7 ev Pitch pal ps PointSize pt 5

#====================================

set origin px+dx*5, dy*0

set title "56Ni"

set format y ""

set cbrange [0:1]
set cbtic 0, 1, 1
set mcbtic 1

set label 1 file at xmax*sinv*1, 1 left

pl \
file u (\$1*sinv):(\$2*sinv):19 ev Pitch pal ps PointSize pt 5

#====================================

pause 100000
EOF
