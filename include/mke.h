#ifndef MKE_H
#define MKE_H

#include <cstdio>
#include <cstdlib>
#include <cmath>


extern double* q_old;
extern double tau;
extern int n_time;

// ---------- Глобальные параметры ----------
extern bool basis;           // true — кубический, false — линейный
extern int kol_elem;         // количество КЭ
extern int left_edge;        // тип КУ слева
extern int right_edge;       // тип КУ справа
extern int N;                // размер СЛАУ

extern const double GAMMA;
extern const double BETA;

// ---------- Глобальные массивы ----------
extern int* ig;
extern double* di;
extern double* gg;
extern double* f;
extern double* node;
extern double* q;

// ---------- Прототипы функций ----------
// solver
int gen_mat();
void clear_mat();
void gen_matrix_mass();
void gen_matrix_zest();
void gen_right_vector();
int LLT();
int gauss();
int solve();
int output_solution(const char* filename);

// mesh
int generate_mesh(double a, double b, int n_points, bool uniform,
                  int left_edge = 1, int right_edge = 1);
// plot
int save_plot_data(const char* filename);
int generate_plot_script(const char* script_name);

// physics
double lambda(double r);
double f_func(double r);
double u_func(double r);
double du_func(double r);

// boundary conditions
void apply_boundary_conditions();

// initial condition
double u0_func(double x);

// time-dependent right-hand side
void gen_right_vector_time(double t);
void apply_boundary_conditions_time(double t);

int output_solution(const char* filename, double t);
const double PI = 3.141592653589793;
#endif // MKE_H