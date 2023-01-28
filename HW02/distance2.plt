# plot.plt
set term png
set output "r2.png"
set title " "
set grid
set autoscale
set xlabel "time "
set ylabel "distance of the spring "
plot "fort.22" using 1:2 title "fort.22"