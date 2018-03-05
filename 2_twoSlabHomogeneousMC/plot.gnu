
set style line 1 lt rgb "cyan" lw 3 pt 6
set style line 2 lt rgb "red" lw 3 pt 6

set title "Homogeneous 2-material neutron flux (1e7 neutrons)"
set xlabel "Distance (cm)"
set ylabel "Normalized flux (neutrons per cm)"
show xlabel
show ylabel
show title
plot "collisions.txt" using 1:2:4 with yerrorbars ls 1 title "Collision based",\
 "collisions.txt" using 1:3:5 with yerrorbars ls 2 title "Distance based", \
 "collisions.txt" using 1:2 with lines ls 1 notitle,\
 "collisions.txt" using 1:3 with lines ls 2 notitle
set key right top
set term png 
set output 'plot.png'
replot
set term x11
