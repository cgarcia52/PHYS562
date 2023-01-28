# plot.plt
set term png
set output "plot.png"
set title "planet orbit"
set grid
set xlabel " "
set ylabel " "
set autoscale
plot "fort.22" using 1:2 title "Flow 1", "fort.21" using 1:2 title "Flow 2", "fort.23" using 1:2 title "Flow 3", "fort.77" using 1:2 title "Flow 4"