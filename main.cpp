#include "include/mke.h"
#include <unistd.h>

// Параметры: ./fem [a] [b] [n_points] [uniform] [basis]
int main(int argc, char* argv[]) {
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        printf("Рабочая директория: %s\n", cwd);
    }

    double a = 0.0, b = 1.0;
    int n_points = 6;
    bool uniform = true;

    if (argc >= 5) {
        a = atof(argv[1]);
        b = atof(argv[2]);
        n_points = atoi(argv[3]);
        uniform = (atoi(argv[4]) != 0);
        if (argc >= 6) basis = (atoi(argv[5]) != 0);
    }

    // 1. Генерация сетки
    if (generate_mesh(a, b, n_points, uniform)) {
        printf("❌ Ошибка генерации сетки\n");
        return 1;
    }
    printf("✅ Сетка [%g, %g], %d точек (%s)\n", a, b, n_points, uniform ? "равномерная" : "случайная");
    
    printf("basis = %s\n", basis ? "кубический" : "линейный");
    printf("N = %d\n", N);

    // 2. Сборка СЛАУ
    if (gen_mat()) return 1;
    printf("В main после gen_mat: N=%d\n", N);

    printf("До clear_mat: di[0]=%e, di[1]=%e\n", di[0], di[1]);
    clear_mat();
    printf("После clear_mat: di[0]=%e, di[1]=%e\n", di[0], di[1]);
    gen_matrix_mass();
    gen_matrix_zest();
    gen_right_vector();

    // 3. Краевые условия
    apply_boundary_conditions();

    // 4. Решение
    if (solve()) {
        printf("❌ Ошибка решения\n");
        return 1;
    }
    printf("После solve: q[0]=%e, q[1]=%e, q[2]=%e\n", q[0], q[1], q[2]);

    // 5. Вывод
    output_solution("result.txt");
    save_plot_data("plot.dat");
    generate_plot_script("plot.sh");

    printf("✅ Готово. Результат: result.txt, график: ./plot.sh\n");
    return 0;
}