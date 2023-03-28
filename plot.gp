set xrange[0:0.7];
set xlabel "V"
set ylabel "P"
plot "carnot" using 3:2 with lines lw 3 title "Carnot", "stirling" using 3:2 with lines lw 3 title "Stirling", "ericsson" using 3:2 with lines lw 3 title "Ericsson"

pause -1
