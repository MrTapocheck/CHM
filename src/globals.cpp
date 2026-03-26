#include "../include/mke.h"

// ЕДИНСТВЕННОЕ определение всех глобальных переменных
bool basis = false;
int kol_elem = 0;
int left_edge = 1;
int right_edge = 1;
int N = 0;

const double GAMMA = 1.0;  // ← КРИТИЧЕСКИ ВАЖНО: γ = 1 для точного решения e^(-π²t)
const double BETA = 1.0;

int* ig = nullptr;
double* di = nullptr;
double* gg = nullptr;
double* f = nullptr;
double* node = nullptr;
double* q = nullptr;

double* q_old = nullptr;
double tau = 0.1;      // шаг по времени
int n_time = 1;        // число шагов