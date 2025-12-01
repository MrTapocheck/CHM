#!/bin/bash
gnuplot -p << 'EOF'
set grid
set multiplot layout 2,1 title 'МКЭ: линейный базис' font ',12'
set title 'Решение'
plot 'plot.dat' u 1:2 w p pt 7 ps 1 lc 'red' t 'МКЭ', \
     '' u 1:3 w l lw 2 lc 'blue' t 'Точное'
set title 'Абсолютная ошибка'
plot 'plot.dat' u 1:4 w lp pt 7 ps 0.8 lc 'purple' t '|ошибка|'
unset multiplot
EOF
