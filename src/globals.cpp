#include "../include/mke.h"

// ЕДИНСТВЕННОЕ определение всех глобальных переменных
bool basis = false;
int kol_elem = 0;
int left_edge = 1;
int right_edge = 1;
int N = 0;

const double GAMMA = 0.0;
const double BETA = 1.0;
int* ig = nullptr;
double* di = nullptr;
double* gg = nullptr;
double* f = nullptr;
double* node = nullptr;
double* q = nullptr;
