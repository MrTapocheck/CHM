#!/bin/bash
echo "Тест на порядок сходимости по времени"
echo "======================================"
echo "tau    | ошибка    | порядок"
echo "-------|-----------|---------"

# Фиксируем N=500 (мелкая сетка)
N=500
T=0.1  # момент времени для сравнения

# Точное решение в центре при t=0.1
u_exact=$(echo "scale=10; e(-3.1415926535^2*0.1)*sin(3.1415926535*0.5)" | bc -l)
echo "Точное решение в центре: $u_exact"

prev_err=0

for tau in 0.1 0.05 0.025 0.0125; do
    n_time=$(echo "scale=0; $T / $tau" | bc)
    
    # Запускаем расчёт
    ./fem 0 1 $N 1 1 1 > /dev/null 2>&1
    
    # Читаем амплитуду из последнего файла
    u_num=$(grep "0.50000000" result_t*.txt | tail -1 | awk '{print $2}')
    
    # Считаем ошибку
    err=$(echo "scale=10; sqrt(($u_num - $u_exact)^2)" | bc -l)
    
    # Считаем порядок
    if [ "$prev_err" != "0" ]; then
        order=$(echo "scale=4; l($prev_err / $err) / l(2)" | bc -l)
        printf "%.4f | %.2e | %.2f\n" $tau $err $order
    else
        printf "%.4f | %.2e | —\n" $tau $err
    fi
    
    prev_err=$err
done