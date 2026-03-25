#!/bin/bash
gnuplot -p << 'EOF'
set grid
set xlabel 'x'
set ylabel 'u(x,t)'
set title 'Нестационарный МКЭ: затухание решения'
plot 'result_t0000.txt' u 1:2 w lp pt 7 ps 0.8 lc 1 t 't=0.0', \
     'result_t0002.txt' u 1:2 w lp pt 7 ps 0.8 lc 2 t 't=0.2', \
     'result_t0004.txt' u 1:2 w lp pt 7 ps 0.8 lc 3 t 't=0.4', \
     'result_t0006.txt' u 1:2 w lp pt 7 ps 0.8 lc 3 t 't=0.6', \
     'result_t0008.txt' u 1:2 w lp pt 7 ps 0.8 lc 3 t 't=0.8', \
     'result_t0009.txt' u 1:2 w lp pt 7 ps 0.8 lc 4 t 't=1.0'
EOF