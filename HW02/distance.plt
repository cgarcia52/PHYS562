# plot.plt
set term png
set output "r.png"
set title " "
set grid
set autoscale
set xlabel "time "
set ylabel "distance of the spring "
plot "fort.11" using 1:2 title "fort.11"