# plot.plt
set term png
set output "plot1.png"
set title "Harmonic"
set grid
set xlabel " "
set ylabel " "
set autoscale
plot "hdata.txt" using 1:2 title "Flow 1", "hdatatwo.txt" using 1:2 title "Flow 2", "hdatathree.txt" using 1:2 title "Flow 3"