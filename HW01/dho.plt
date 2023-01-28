# plot.plt
set term png
set output "plot.png"
set title "D Harmnonic Osciilator"
set grid
set autoscale
set xlabel " "
set ylabel " "
plot "fort.1" using 1:2 title "Flow 1", "fort.2" using 1:2 title "Flow 2", "fort.3" using 1:2 title "Flow 3" 