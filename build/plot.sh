#!/bin/bash
gnuplot -p << 'EOF'
set grid
set xlabel 'x'
set ylabel 'u(x,t)'
set title 'Нестационарный МКЭ: затухание решения'
plot 'result_t0000.txt' u 1:2 w lp pt 7 ps 0.8 lc 1 t 't=0.0', \
     'result_t0020.txt' u 1:2 w lp pt 7 ps 0.8 lc 2 t 't=0.2', \
     'result_t0040.txt' u 1:2 w lp pt 7 ps 0.8 lc 3 t 't=0.4', \
     'result_t0099.txt' u 1:2 w lp pt 7 ps 0.8 lc 4 t 't=1.0'
EOF