#!/usr/bin/gnuplot

set encoding utf8
set term pngcairo size 600*1.5/2.0, 400*1.5/2.0 font "Helvetica,"
set output "hw2_part2.png"

set datafile separator ","
#set key autotitle columnhead

set xlabel "log_2h"
set ylabel "log_2e"

unset key

plot "hw2_part2.csv" using 2:4 with points lw 2,\
    "hw2_part2.csv" using 2:6 with points lw 2;
