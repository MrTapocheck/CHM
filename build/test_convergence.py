#!/usr/bin/env python3
import subprocess
import math
import os

print("Тест на порядок сходимости по времени")
print("=" * 50)
print(f"{'tau':<10} | {'ошибка':<12} | {'порядок':<10}")
print("-" * 50)

# Точное решение в центре при t=0.1
T = 0.1
x_center = 0.5
u_exact = math.exp(-math.pi**2 * T) * math.sin(math.pi * x_center)
print(f"Точное решение в центре (t={T}): {u_exact:.6e}\n")

prev_err = None

for tau in [0.1, 0.05, 0.025, 0.0125]:
    n_time = int(T / tau)
    
    # Запускаем расчёт
    cmd = f"./fem 0 1 500 1 1 1"
    subprocess.run(cmd, shell=True, capture_output=True)
    
    # Читаем последний файл
    result_file = f"result_t{n_time-1:04d}.txt"
    if not os.path.exists(result_file):
        print(f"❌ Файл {result_file} не найден")
        continue
    
    with open(result_file, 'r') as f:
        lines = f.readlines()
        # Ищем строку с центром (близко к 0.5)
        for line in lines[1:]:  # пропускаем заголовок
            parts = line.split()
            if len(parts) >= 2:
                x = float(parts[0])
                if abs(x - x_center) < 0.01:  # близко к центру
                    u_num = float(parts[1])
                    break
    
    # Считаем ошибку
    err = abs(u_num - u_exact)
    
    # Считаем порядок
    if prev_err is not None:
        order = math.log2(prev_err / err)
        print(f"{tau:<10.4f} | {err:<12.2e} | {order:<10.2f}")
    else:
        print(f"{tau:<10.4f} | {err:<12.2e} | {'—':<10}")
    
    prev_err = err

print("-" * 50)