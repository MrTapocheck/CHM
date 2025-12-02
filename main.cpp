#include "include/mke.h"
#include <unistd.h>

// Параметры: ./fem [a] [b] [n_points] [uniform] [basis]
int main(int argc, char* argv[]) {
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        printf("Рабочая директория: %s\n", cwd);
    }




    auto get_arg = [&](int i, auto def) {
        return (argc > i) ? static_cast<decltype(def)>(atof(argv[i])) : def;
    };
    //по умолчанию a = 0 b = 1
    double a = get_arg(1, 0.0); 
    double b = get_arg(2, 1.0);
    int n_points = get_arg(3, 6); //5 точек
    bool uniform = get_arg(4, 1) != 0;//сетка равномерная по умолчанию
    //краевые по умолчанию первые
    int le = get_arg(5, 1);
    int re = get_arg(6, 1);
    basis = get_arg(7, 0) != 0; //базис по умолчанию линейный

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

    printf("До clear_mat: di[0]=%e, di[1]=%e\n", di[0], di[1]);
    clear_mat();
    printf("После clear_mat: di[0]=%e, di[1]=%e\n", di[0], di[1]);
    // gen_matrix_mass();
    gen_matrix_zest();
    gen_right_vector();

    // 3. Краевые условия
    printf("До BC: di[0]=%e, f[0]=%e\n", di[0], f[0]);
    apply_boundary_conditions();
    printf("После BC: di[0]=%e, f[0]=%e\n", di[0], f[0]);


    // printf("N=%d, ig[N]=%d, gg size=%d\n", N, ig[N], ig[N]-1);
    // for (int i = 1; i < N; i++) {
    //     printf("Строка %d: ширина = %d\n", i, ig[i+1]-ig[i]);
    // }
    // 4. Решение
    if (solve()) {
        printf("❌ Ошибка решения\n");
        return 1;
    }
    // printf("N=%d, ig[N]=%d, gg size=%d\n", N, ig[N], ig[N]-1);
    // for (int i = 1; i < N; i++) {
    //     printf("Строка %d: ширина = %d\n", i, ig[i+1]-ig[i]);
    // }

    printf("После solve: q[0]=%e, q[1]=%e, q[2]=%e\n", q[0], q[1], q[2]);

    // 5. Вывод
    output_solution("result.txt");
    save_plot_data("plot.dat");
    generate_plot_script("plot.sh");

    printf("✅ Готово. Результат: result.txt, график: ./plot.sh\n");
    return 0;
}