#include "include/mke.h"
#include <unistd.h>
#include <cmath>

int main(int argc, char* argv[]) {
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        printf("Рабочая директория: %s\n", cwd);
    }

    auto get_arg = [&](int i, auto def) {
        return (argc > i) ? static_cast<decltype(def)>(atof(argv[i])) : def;
    };

    double a = get_arg(1, 0.0);
    double b = get_arg(2, 1.0);
    int n_points = get_arg(3, 6);
    bool uniform = get_arg(4, 1) != 0;
    int le = get_arg(5, 1);
    int re = get_arg(6, 1);
    basis = get_arg(7, 0) != 0;

    // 1. Генерация сетки
    if (generate_mesh(a, b, n_points, uniform, le, re)) {
        printf("❌ Ошибка генерации сетки\n");
        return 1;
    }
    printf("✅ Сетка [%g, %g], %d точек (%s)\n", a, b, n_points, uniform ? "равномерная" : "случайная");
    printf("Левое краевое - %d, Правое краевое - %d\n", le, re);
    printf("basis = %s\n", basis ? "кубический" : "линейный");
    printf("N = %d\n", N);

    // 2. Инициализация СЛАУ
    if (gen_mat()) return 1;

    // === НЕСТАЦИОНАРНЫЙ ЦИКЛ ПО ВРЕМЕНИ ===
    q_old = new double[N];

    // Начальное условие: u(x,0) = sin(πx)
    for (int i = 0; i < N; i++) {
        q_old[i] = sin(PI * node[i]);
    }

    // ПАРАМЕТРЫ ДЛЯ ТЕСТА (менять вручную):
    tau = 0.0125;
    n_time = 8;

    printf("🚀 Начинаю расчёт: tau=%.4f, n_time=%d, T=%.3f\n", tau, n_time, tau*n_time);

    for (int n = 0; n < n_time; n++) {
        double t = (n + 1) * tau;

        // Собираем матрицу A = (γ/τ)M + K
        clear_mat();
        gen_matrix_mass();   // ← ПОЛНАЯ масса с учётом GAMMA
        gen_matrix_zest();   // ← жёсткость

        // Умножаем ВСЁ на 1/τ (γ уже внутри массы!)
        for (int i = 0; i < N; i++) di[i] *= (1.0 / tau);
        for (int i = 0; i < ig[N]-1; i++) gg[i] *= (1.0 / tau);

        // Правая часть: b = (γ/τ) * M * q_old
        for (int i = 0; i < N; i++) f[i] = 0.0;
        
        for (int i = 0; i < kol_elem; i++) {
            double h = node[i+1] - node[i];
            double m00 = GAMMA * h / 3.0;
            double m01 = GAMMA * h / 6.0;
            double m11 = GAMMA * h / 3.0;
            
            f[i]     += (m00 * q_old[i] + m01 * q_old[i+1]) / tau;
            f[i+1]   += (m01 * q_old[i] + m11 * q_old[i+1]) / tau;
        }

        // Краевые условия
        apply_boundary_conditions();

        // Решаем
        if (solve()) {
            printf("❌ Ошибка на шаге %d\n", n);
            delete[] q_old;
            return 1;
        }

        // Сохраняем ТОЛЬКО последний шаг (в момент времени Т)
        if (n == n_time - 1) {
            char fname[64];
            sprintf(fname, "result_final.txt");
            output_solution(fname, t);  // ← ПЕРЕДАЁМ ВРЕМЯ!
            printf("💾 t=%.3f: u(0.5)=%.6f\n", t, q[N/2]);
        }

        // q -> q_old
        for (int i = 0; i < N; i++) {
            q_old[i] = q[i];
        }
    }

    delete[] q_old;
    printf("✅ Готово. Файл: result_final.txt\n");
    return 0;
}