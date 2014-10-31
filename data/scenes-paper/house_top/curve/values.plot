clear
reset
#set nokey
set terminal svg size 600,150 
set termoption dashed
set xrange[1:9]
set yrange[0:1.2]
#unset border
#unset tics
set style line 1 lt 1 lc rgb "black" lw 2
set style line 2 lt 1 lc rgb "black" lw 2
set style line 3 lt 3 lc rgb "slategray" lw 2

set style line 4 lt 1 lc rgb "black" lw 4

set key right bottom

set output "/home/nmellado/git/aron/globOpt/data/scenes-paper/house_top/curve/curves.svg"

plot "./values.csv" u 1:($2/0.220595941775) w lines ls 1 t "$\scriptstyle \\sigma / \\sigma_{GT}$", "./values.csv" u 1:3 w lines ls 2 t "F-measure (relations)"
