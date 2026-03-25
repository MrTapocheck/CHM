#include "include/mke.h"
#include <unistd.h>

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

    // 2. Сборка СЛАУ
    if (gen_mat()) return 1;
    printf("В main после gen_mat: N=%d\n", N);

    clear_mat();
    gen_matrix_zest();
    gen_right_vector();

    // 3. Краевые условия
    apply_boundary_conditions();

    // === НЕСТАЦИОНАРНЫЙ ЦИКЛ ПО ВРЕМЕНИ ===
    q_old = new double[N];

    // Начальное условие: u(x,0) = sin(pi*x)
    for (int i = 0; i < N; i++) {
        q_old[i] = sin(3.1415926535 * node[i]);
    }

    tau = 0.01;
    n_time = 100;

    printf("🚀 Начинаю нестационарный расчёт: %d шагов, tau=%.3f\n", n_time, tau);

    for (int n = 0; n < n_time; n++) {
        double t = (n + 1) * tau;

        // Собираем матрицу A = (1/tau)*I + K
        clear_mat();
        gen_matrix_zest();

        // Умножаем диагональ на 1/tau (M = I, т.к. gamma=0)
        for (int i = 0; i < N; i++) {
            di[i] += 1.0 / tau;
        }

        // Правая часть: b = q_old / tau
        for (int i = 0; i < N; i++) {
            f[i] = q_old[i] / tau;
        }

        // Краевые условия
        apply_boundary_conditions();

        // Решаем
        if (solve()) {
            printf("❌ Ошибка на шаге %d\n", n);
            delete[] q_old;
            return 1;
        }

        // Сохраняем результат каждые 20 шагов
        if (n % 20 == 0 || n == n_time - 1) {
            char fname[64];
            sprintf(fname, "result_t%04d.txt", n);
            output_solution(fname);
            printf("💾 Сохранено: %s (t=%.2f, u(0.5)=%.4e)\n", fname, t, q[N/2]);
        }

        // q -> q_old
        for (int i = 0; i < N; i++) {
            q_old[i] = q[i];
        }
    }

    delete[] q_old;
    printf("✅ Нестационарный расчёт завершён. Файлы: result_t*.txt\n");
    return 0;
}