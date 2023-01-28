# plot.plt
set term png
set output "plot.png"
set title "Lotka-Volterra Equations"
set grid
set xlabel "Time"
set ylabel "Population"
plot "fort.17" using 1:2, "fort.19" using 1:2