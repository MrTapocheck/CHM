#include "../include/mke.h"
#include <cmath>
#include <cstdio>
#include <cstring>



// Коэффициенты уравнения
double lambda(double r) { return 1.0; }  // теплопроводность
double f_func(double r) { return 0.0; }   // правая часть = 0 для точного решения

// Начальное условие (при t=0)
double u_func(double r) { 
    return sin(PI * r);  // u(x,0) = sin(πx)
}

double du_func(double r) { 
    return PI * cos(PI * r);  // производная для краевых условий 2-го рода
}

// Генерация матрицы и сетки
int gen_mat() {
    printf("gen_mat: kol_elem=%d\n", kol_elem);
    FILE* in;
    int i, j;
    
    in = fopen("param.txt", "r");
    if (!in) {
        perror("fopen param.txt");
        return 1;
    }
    if (fscanf(in, "%d", &kol_elem) != 1) { fclose(in); return 1; }
    if (fscanf(in, "%d", &left_edge) != 1) { fclose(in); return 1; }
    if (fscanf(in, "%d", &right_edge) != 1) { fclose(in); return 1; }
    fclose(in);
    
    printf("ПОСЛЕ чтения: kol_elem=%d\n", kol_elem);
    
    N = basis ? (1 + kol_elem * 3) : (kol_elem + 1);
    
    ig = new int[N + 1]();
    if (!ig) return 1;
    
    ig[0] = 1; ig[1] = 1;
    if (basis) { // кубические
        for (i = 0, j = 2; i < kol_elem; i++) {
            ig[j] = ig[j-1] + 1; j++;
            ig[j] = ig[j-1] + 2; j++;
            ig[j] = ig[j-1] + 3; j++;
        }
    } else { // линейные
        for (i = 0, j = 2; i < kol_elem; i++) {
            ig[j] = ig[j-1] + 1; j++;
        }
    }
    
    di = new double[N](); 
    f = new double[N](); 
    q = new double[N]();
    gg = new double[ig[N] - 1](); 
    node = new double[kol_elem + 1]();
    
    if (!di || !f || !q || !gg || !node) return 1;
    
    in = fopen("setka.txt", "r");
    if (!in) return 1;
    for (i = 0; i <= kol_elem; i++) {
        if (fscanf(in, "%lf", &node[i]) != 1) { fclose(in); return 1; }
    }
    fclose(in);
    
    printf("ig[N] = %d, размер gg = %d\n", ig[N], ig[N] - 1);
    return 0;
}

// Очистка матрицы
void clear_mat() {
    printf("Очистка: N=%d\n", N);
    for (int i = 0; i < N; i++) { 
        f[i] = 0.0; 
        di[i] = 0.0; 
    }
    for (int i = 0; i < ig[N] - 1; i++) {
        gg[i] = 0.0;
    }
}

// ПОЛНАЯ (консистентная) матрица массы
void gen_matrix_mass() {
    printf("=== gen_matrix_mass ===\n");
    for (int i = 0; i < kol_elem; i++) {
        double h = node[i+1] - node[i];
        
        // Диагональные элементы: ∫φ_i² dx = h/3
        di[i]     += GAMMA * h / 3.0;
        di[i+1]   += GAMMA * h / 3.0;
        
        // Внедиагональные элементы: ∫φ_i φ_{i+1} dx = h/6
        if (i+1 < N) {
            int idx = ig[i+2] - 2;
            if (idx >= 0 && idx < ig[N]-1) {
                gg[idx] += GAMMA * h / 6.0;
            }
        }
    }
}

// Матрица жёсткости
void gen_matrix_zest() {
    printf("=== gen_matrix_zest ===\n");
    for (int i = 0; i < kol_elem; i++) {
        double h = node[i + 1] - node[i];
        double k = lambda(0.5 * (node[i] + node[i+1])) / h;  // λ/h
        
        di[i]     += k;
        di[i + 1] += k;
        
        if (i + 1 < N) {
            int idx = ig[i + 2] - 2;
            if (idx >= 0 && idx < ig[N] - 1) {
                gg[idx] -= k;
            }
        }
    }
}

// Правая часть (для стационарной задачи)
void gen_right_vector() {
    printf("=== gen_right_vector ===\n");
    for (int i = 0; i < kol_elem; i++) {
        double h = node[i + 1] - node[i];
        double f0 = f_func(node[i]);
        double f1 = f_func(node[i + 1]);
        
        f[i]     += (f0 + f1) * h / 4.0;
        f[i + 1] += (f0 + f1) * h / 4.0;
    }
}

// Краевые условия
void left_edge_1() {
    di[0] = 1e12;
    f[0] = 1e12 * u_func(node[0]);  // u(0) = 0
}

void right_edge_1() {
    di[N-1] = 1e12;
    f[N-1] = 1e12 * u_func(node[kol_elem]);  // u(1) = 0
}

void left_edge_2() { 
    f[0] += -lambda(node[0]) * du_func(node[0]); 
}

void right_edge_2() { 
    f[N-1] += lambda(node[kol_elem]) * du_func(node[kol_elem]); 
}

void left_edge_3() { 
    di[0] += BETA; 
    f[0] += (-lambda(node[0])*du_func(node[0])*node[0]*node[0])/BETA + u_func(node[0]); 
}

void right_edge_3() { 
    di[N-1] += BETA; 
    f[N-1] += (lambda(node[kol_elem])*du_func(node[kol_elem])*node[kol_elem]*node[kol_elem])/BETA + u_func(node[kol_elem]); 
}

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

// LLT разложение
int LLT() {
    int i, j, k;
    double per;
    int a, b;
    
    di[0] = sqrt(di[0]);
    
    for (i = 1; i < N; i++) {
        a = i - ig[i+1] + ig[i];  // начало i-й строки
        
        for (j = 0; j < ig[i+1] - ig[i]; j++) {
            b = a + j - ig[a+j+1] + ig[a+j];  // начало j-й строки
            
            per = gg[ig[i] + j - 1];
            
            if (a < b) {  // i-я строка раньше j-й
                for (k = ig[a+j+1] - ig[a+j] - 1; k >= 0; k--) {
                    per -= gg[ig[a+j] + k - 1] * gg[ig[i] + b - a + k - 1];
                }
            } else {  // j-я строка раньше i-й
                for (k = a - b; k < ig[a+j+1] - ig[a+j]; k++) {
                    per -= gg[ig[i] + k - 1 - (a - b)] * gg[ig[a+j] + k - 1];
                }
            }
            
            gg[ig[i] + j - 1] = per / di[a + j];
        }
        
        // Диагональ
        per = di[i];
        for (k = 0; k < ig[i+1] - ig[i]; k++) {
            per -= gg[ig[i] + k - 1] * gg[ig[i] + k - 1];
        }
        di[i] = sqrt(per);
        
        if (di[i] == 0) {
            printf("LLT: нулевой диагональный элемент в строке %d\n", i);
            return 1;
        }
    }
    
    // Прямой ход: L y = f
    double* y = new double[N];
    for (i = 0; i < N; i++) {
        double sum = f[i];
        for (j = 0; j < ig[i+1] - ig[i]; j++) {
            int col = i - (ig[i+1] - ig[i]) + j;
            sum -= gg[ig[i] + j - 1] * y[col];
        }
        y[i] = sum / di[i];
    }
    
    // Обратный ход: L^T x = y
    for (i = N - 1; i >= 0; i--) {
        double sum = y[i];
        for (j = 0; j < ig[i+1] - ig[i]; j++) {
            int col = i - (ig[i+1] - ig[i]) + j;
            sum -= gg[ig[i] + j - 1] * q[col];
        }
        q[i] = sum / di[i];
    }
    
    delete[] y;
    return 0;
}

// Метод Гаусса
int gauss() {
    int i, j;
    
    // Прямой ход: L * y = f
    q[0] = f[0] / di[0];
    for (i = 1; i < N; i++) {
        q[i] = f[i];
        for (j = 0; j < ig[i+1] - ig[i]; j++) {
            int col = i - (ig[i+1] - ig[i]) + j;
            q[i] -= gg[ig[i] + j - 1] * q[col];
        }
        q[i] = q[i] / di[i];
    }
    
    // Копируем q → f для обратного хода
    for (i = 0; i < N; i++) {
        f[i] = q[i];
    }
    
    // Обратный ход: L^T * x = y
    for (i = N - 1; i >= 0; i--) {
        q[i] = f[i] / di[i];
        for (j = 0; j < ig[i+1] - ig[i]; j++) {
            int col = i - (ig[i+1] - ig[i]) + j;
            f[col] -= q[i] * gg[ig[i] + j - 1];
        }
    }
    
    return 0;
}

int solve() {
    if (LLT()) return 1;
    if (gauss()) return 1;
    return 0;
}

// Вывод решения с точным решением для времени t
int output_solution(const char* filename, double t) {
    FILE* out = fopen(filename, "w");
    if (!out) return 1;
    
    fprintf(out, "x\tapprox\tanalytic\terror\n");
    
    for (int i = 0; i <= kol_elem; i++) {
        double approx = basis ? q[3*i] : q[i];
        
        // ТОЧНОЕ РЕШЕНИЕ для уравнения теплопроводности
        // с начальным условием u(x,0) = sin(πx) и краями u(0)=u(1)=0
        double exact = exp(-PI*PI * t) * sin(PI * node[i]);
        
        double err = fabs(approx - exact);
        fprintf(out, "%.8f\t%.8f\t%.8f\t%.2e\n", node[i], approx, exact, err);
    }
    
    fclose(out);
    return 0;
}