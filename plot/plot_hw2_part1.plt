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
plot for [i=2:4] "hw2_part1_o1.csv" using 1:i with lines lw 2;
unset label 1
unset key

set label 1 'b' @POS
plot for [i=2:4] "hw2_part1_o2.csv" using 1:i with lines lw 2;
unset label 1

set label 1 'c' @POS
plot [0:1] cos(pi*x/2)*exp(-pi*pi/4*0.06) with lines lw 2,\
    cos(pi*x/2)*exp(-pi*pi/4*0.1) with lines lw 2,\
    cos(pi*x/2)*exp(-pi*pi/4*0.9) with lines lw 2;
unset label 1
