#!/bin/bash
gnuplot -p << 'EOF'
set grid
set multiplot layout 3,1 title 'МКЭ: линейный базис' font ',12'
set title 'Решение'
plot 'plot.dat' u 1:2 w p pt 7 ps 1 lc 'red' t 'МКЭ', '' u 1:3 w l lw 2 lc 'blue' t 'Точное'
set title 'Абсолютная ошибка'
plot 'plot.dat' u 1:4 w lp pt 7 ps 0.7 lc 'purple' t '|ош|'
set title 'log_{10}(|ош|)'
set logscale y
plot 'plot.dat' u 1:(log10($4 + 1e-30)) w lp pt 7 ps 0.7 lc 'green' t 'log_{10}(|ош|)'
unset multiplot
EOF
