#pragma once
#include "const.h"

double enthalpy_next(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double a1, double Pre);
double enthalpy(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double a1, double Pre);