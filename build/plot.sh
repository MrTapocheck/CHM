#!/bin/bash
gnuplot -p << 'EOF'
set grid
set xlabel 'x'
set ylabel 'u(x,t)'
set title 'Нестационарный МКЭ: затухание решения'
plot 'result_t0000.txt' u 1:2 w lp pt 7 ps 0.8 lc 1 t 't=0.5', \
     'result_t0001.txt' u 1:2 w lp pt 7 ps 0.8 lc 2 t 't=1.0', \
     'result_t0002.txt' u 1:2 w lp pt 7 ps 0.8 lc 3 t 't=1.5', \
     'result_t0003.txt' u 1:2 w lp pt 7 ps 0.8 lc 4 t 't=2.0'
EOF