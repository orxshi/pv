set xrange[0:7000];
set xlabel "V (cc)"
set ylabel "P (kPa)"

array ARGV[ARGC]

do for [i=1:ARGC] {
	eval sprintf("ARGV[%d] = ARG%d", i, i);
}

#set label "1" at 0.6, 300 point pointtype 7 pointsize 1
plot for [i=1:ARGC] ARGV[i] using ($3*10000):($2/1000) with lines lw 3 title ARGV[i]

pause -1
