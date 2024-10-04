#!/usr/bin/gnuplot

set encoding utf8
set term pngcairo size 600*1.5/2.0, 400*1.5/2.0 font "Helvetica,"
set output "hw2_part3.png"

set xlabel ""

plot [10:50] 1-sin(2*pi/x)/(2*pi/x) with lines lw 2 title "e_2(N)",\
    [10:50] (1+1/3*sin(4*pi/x)/(4*pi/x)-4/3*sin(2*pi/x)/(2*pi/x)) with lines lw 2 title "e_4(N)";
