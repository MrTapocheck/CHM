#include "include/mke.h"
#include <unistd.h>
#include <cmath>

int main(int argc, char* argv[]) {
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        printf("Рабочая директория: %s\n", cwd);
    }

    // Параметры по умолчанию
    double a = 0.0, b = 1.0;
    int n_points = 50;
    bool uniform = true;
    int le = 1, re = 1;
    basis = false;

    // Парсинг аргументов
    if (argc > 1) a = atof(argv[1]);
    if (argc > 2) b = atof(argv[2]);
    if (argc > 3) n_points = atoi(argv[3]);
    if (argc > 4) uniform = (atoi(argv[4]) != 0);
    if (argc > 5) le = atoi(argv[5]);
    if (argc > 6) re = atoi(argv[6]);
    if (argc > 7) basis = (atoi(argv[7]) != 0);

    // Генерация сетки
    if (generate_mesh(a, b, n_points, uniform, le, re)) {
        printf("❌ Ошибка генерации сетки\n");
        return 1;
    }
    printf("✅ Сетка [%g, %g], %d точек (%s)\n", a, b, n_points, uniform ? "равномерная" : "случайная");
    printf("Левое краевое - %d, Правое краевое - %d\n", le, re);
    printf("basis = %s\n", basis ? "кубический" : "линейный");
    printf("N = %d\n", N);

    // Инициализация СЛАУ
    if (gen_mat()) return 1;

    // Нестационарный цикл
    q_old = new double[N];

    // Начальное условие: u(x,0) = sin(πx)
    for (int i = 0; i < N; i++) {
        q_old[i] = sin(3.1415926535 * node[i]);
    }

    tau = 0.1;
    n_time = 10;

    for (int n = 0; n < n_time; n++) {
        double t = (n + 1) * tau;

        // Собираем матрицу A = (1/τ)M + K
        clear_mat();
        
        // ВОТ ОН — ВЫЗОВ МАССЫ! ЭТОТ СТРОКА КРИТИЧЕСКАЯ!
        gen_matrix_mass();
        
        gen_matrix_zest();

        // Умножаем всё на 1/τ
        for (int i = 0; i < N; i++) di[i] *= (1.0 / tau);
        for (int i = 0; i < ig[N]-1; i++) gg[i] *= (1.0 / tau);

        // Правая часть: b = (1/τ) * M_lumped * q_old
        for (int i = 0; i < N; i++) f[i] = 0.0;
        for (int i = 0; i < kol_elem; i++) {
            double h = node[i+1] - node[i];
            f[i]     += (h / 2.0 * q_old[i]) / tau;
            f[i+1]   += (h / 2.0 * q_old[i+1]) / tau;
        }

        // Краевые условия (ТОЛЬКО ОДИН РАЗ!)
        apply_boundary_conditions();

        // Решаем (ТОЛЬКО ОДИН РАЗ!)
        if (solve()) {
            printf("❌ Ошибка на шаге %d\n", n);
            delete[] q_old;
            return 1;
        }

        // Сохраняем
        if (n == 0 || n == 2 || n == 4 || n == 6 || n == 8 || n == 9) {
            char fname[64];
            sprintf(fname, "result_t%04d.txt", n);
            output_solution(fname);
            printf("💾 %s: u(0.5)=%.4f\n", fname, q[N/2]);
        }

        // q -> q_old
        for (int i = 0; i < N; i++) q_old[i] = q[i];
    }

    delete[] q_old;
    printf("✅ Готово. Файлы: result_t*.txt\n");
    return 0;
}