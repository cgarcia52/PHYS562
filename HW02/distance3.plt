# plot.plt
set term png
set output "r3.png"
set title " "
set grid
set autoscale
set xlabel "time "
set ylabel "distance of the spring "
plot "fort.33" using 1:2 title "fort.33"