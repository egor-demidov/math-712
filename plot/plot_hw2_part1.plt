#!/usr/bin/gnuplot

set encoding utf8
set term pngcairo size 600*1.5, 400*1.5 font "Helvetica,"
set output "hw2_part1.png"

set datafile separator ","
set key autotitle columnhead

set multiplot layout 2, 2 rowsfirst

set xlabel "x"
set ylabel "v"
set ytics 0.2

POS = "at graph 0.90,0.90 font ',20'"

set label 1 'a' @POS
set key left center
plot for [i=2:5] "sol_one_sided.csv" using 1:i with lines lw 2;
unset label 1
unset key

set label 1 'b' @POS
plot for [i=2:5] "sol_ghost_point.csv" using 1:i with lines lw 2;
unset label 1

set label 1 'c' @POS
plot for [i=2:5] "sol_exact.csv" using 1:i with lines lw 2;
unset label 1

set ytics 10
set logscale y 10
set logscale x 10
set format y "%.0s*10^{%T}"

set label 1 'd' @POS
set key left center
set xlabel "t"
set ylabel "L^2 norm of error"
plot for [i=2:3] "error.csv" using 1:i with points lw 2;
unset label 1
