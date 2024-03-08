set xrange[0:0.9];
set xlabel "V"
set ylabel "P"

array ARGV[ARGC]

do for [i=1:ARGC] {
	eval sprintf("ARGV[%d] = ARG%d", i, i);
}

#set label "1" at 0.6, 300 point pointtype 7 pointsize 1
plot for [i=1:ARGC] ARGV[i] using 3:2 with lines lw 3 title ARGV[i]

pause -1
