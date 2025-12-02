#include "../include/mke.h"
#include <ctime>
#include <algorithm>

int generate_mesh(double a, double b, int n_points, bool uniform,
                  int left_edge, int right_edge) {
    if (n_points < 2 || a >= b) return 1;

    double* x = new double[n_points];
    if (uniform) {
        double h = (b - a) / (n_points - 1);
        for (int i = 0; i < n_points; i++) x[i] = a + i * h;
    } else {
        std::srand(std::time(nullptr) ^ (size_t)&a);
        x[0] = a; x[n_points - 1] = b;
        for (int i = 1; i < n_points - 1; i++) {
            x[i] = a + (double)std::rand() / RAND_MAX * (b - a);
        }
        std::sort(x, x + n_points);
    }

    // setka.txt
    FILE* f = fopen("setka.txt", "w");
    for (int i = 0; i < n_points; i++) fprintf(f, "%.12f\n", x[i]);
    fclose(f);

    // param.txt
    f = fopen("param.txt", "w");
    fprintf(f, "%d\n%d\n%d\n", n_points - 1, left_edge, right_edge);
    fclose(f);

    delete[] x;
    return 0;
}