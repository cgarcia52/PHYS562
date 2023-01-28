# plot.plt
set term png
set output "double1.png"
set title " "
set grid
set autoscale
set xlabel "time "
set ylabel "distance of the spring "
plot "fort.66" using 1:5