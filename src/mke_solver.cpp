#include "../include/mke.h"

// Инициализация глобальных переменных
bool basis = false; //true - кубический, fasle - линейный
int kol_elem = 0;
int left_edge = 1;
int right_edge = 1;
int N = 0;

const double GAMMA = 1.0;
const double BETA = 1.0;

int* ig = nullptr;
double* di = nullptr;
double* gg = nullptr;
double* f = nullptr;
double* node = nullptr;
double* q = nullptr;

// Физика
double lambda(double r) { return r; }
double f_func(double r) {
    return basis ? (-24. * r * r * r + r * r * r * r) : (-3 + r);
}
double u_func(double r) {
    return basis ? r * r * r * r : r;
}
double du_func(double r) {
    return basis ? 4 * r * r * r : 1;
}

// Генерация матрицы и сетки
int gen_mat() {
    printf("gen_mat: kol_elem=%d\n", kol_elem);
    FILE* in;
    int i, j;

    in = fopen("param.txt", "r");
    if (!in) {
        perror("fopen param.txt");  // ← ПОКАЖЕТ, ГДЕ ИЩЕТ
        return 1;
    }    
//    if ((in = fopen("param.txt", "r")) == nullptr) return 1;
    if (fscanf(in, "%d", &kol_elem) != 1) { fclose(in); return 1; }
    if (fscanf(in, "%d", &left_edge) != 1) { fclose(in); return 1; }
    if (fscanf(in, "%d", &right_edge) != 1) { fclose(in); return 1; }
    fclose(in);

    N = basis ? (1 + kol_elem * 3) : (kol_elem + 1);

    ig = new int[N + 1]();
    if (!ig) return 1;

    ig[0] = 1; ig[1] = 1;
    if (basis) {
        for (i = 0, j = 2; i < kol_elem; i++) {
            ig[j] = ig[j-1] + 1; j++;
            ig[j] = ig[j-1] + 2; j++;
            ig[j] = ig[j-1] + 3; j++;
        }
    } else {
        for (i = 0, j = 2; i < kol_elem; i++) {
            ig[j] = ig[j-1] + 1; j++;
        }
    }

    di = new double[N](); f = new double[N](); q = new double[N]();
    gg = new double[ig[N] - 1](); node = new double[kol_elem + 1]();
    if (!di || !f || !q || !gg || !node) return 1;

    if ((in = fopen("setka.txt", "r")) == nullptr) return 1;
    for (i = 0; i <= kol_elem; i++) {
        if (fscanf(in, "%lf", &node[i]) != 1) { fclose(in); return 1; }
    }
    fclose(in);
    return 0;
}

// Очистка
void clear_mat() {
    for (int i = 0; i < N; i++) { f[i] = 0; di[i] = 0; }
    for (int i = 0; i < ig[N] - 1; i++) gg[i] = 0;
}

// Матрицы
void gen_matrix_mass() { /* ... твои формулы ... */ }
void gen_matrix_zest() { /* ... твои формулы ... */ }
void gen_right_vector() { /* ... твои формулы ... */ }

// Краевые условия
void left_edge_1() {
    if (basis) { di[0] = 1e12; f[0] = 1e12 * u_func(node[0]); }
    else { di[0] = 1e12; f[0] = 1e12 * u_func(node[0]); }
}
void right_edge_1() {
    int idx = basis ? (3 * kol_elem) : (N - 1);
    di[idx] = 1e12; f[idx] = 1e12 * u_func(node[kol_elem]);
}
void left_edge_2() { f[0] += -lambda(node[0]) * du_func(node[0]) * node[0] * node[0]; }
void right_edge_2() { f[N-1] += lambda(node[kol_elem]) * du_func(node[kol_elem]) * node[kol_elem] * node[kol_elem]; }
void left_edge_3() { di[0] += BETA; f[0] += (-lambda(node[0])*du_func(node[0])*node[0]*node[0])/BETA + u_func(node[0]); }
void right_edge_3() { di[N-1] += BETA; f[N-1] += (lambda(node[kol_elem])*du_func(node[kol_elem])*node[kol_elem]*node[kol_elem])/BETA + u_func(node[kol_elem]); }

void apply_boundary_conditions() {
    switch (left_edge) {
        case 1: left_edge_1(); break;
        case 2: left_edge_2(); break;
        case 3: left_edge_3(); break;
    }
    switch (right_edge) {
        case 1: right_edge_1(); break;
        case 2: right_edge_2(); break;
        case 3: right_edge_3(); break;
    }
}

// Решение
int LLT() { /* ... плотный Холецкий ... */ return 0; }
int gauss() { for (int i = 0; i < N; i++) f[i] = q[i]; return 0; }
int solve() {
    if (LLT()) return 1;
    if (gauss()) return 1;
    return 0;
}

// Вывод
int output_solution(const char* filename) {
    FILE* out = fopen(filename, "w");
    if (!out) return 1;
    fprintf(out, "x\tapprox\tanalytic\terror\n");
    for (int i = 0; i <= kol_elem; i++) {
        double approx = basis ? q[3*i] : q[i];
        double exact = u_func(node[i]);
        double err = fabs(approx - exact);
        fprintf(out, "%.8f\t%.8f\t%.8f\t%.2e\n", node[i], approx, exact, err);
    }
    fclose(out);
    return 0;
}