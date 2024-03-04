set xrange[0:0.9];
set xlabel "V"
set ylabel "P"

array ARGV[ARGC]

do for [i=1:ARGC] {
	eval sprintf("ARGV[%d] = ARG%d", i, i);
}

plot for [i=1:ARGC] ARGV[i] using 3:2 with lines lw 3 title ARGV[i]

pause -1
