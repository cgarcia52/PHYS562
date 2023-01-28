# plot.plt
set term png
set output "theta.png"
set title " "
set grid
set autoscale
set xlabel "time "
set ylabel "angle "
plot "fort.1" using 1:2 title " ", "fort.2" using 1:2 title " "