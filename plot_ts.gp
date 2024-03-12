set xrange[0:1000];
set xlabel "s"
set ylabel "T (K)"

array ARGV[ARGC]

do for [i=1:ARGC] {
	eval sprintf("ARGV[%d] = ARG%d", i, i);
}

plot for [i=1:ARGC] ARGV[i] using 5:4 with lines lw 3 title ARGV[i]

pause -1
