#pragma once
#include "const.h"

void metropolis_ads_next(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double a1, double beta);
void metropolis_ads_next_exchange(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double a1, double beta);
void metropolis_ads(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double a1, double beta);
void metropolis_position_next(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double a1, double beta);
void metropolis_position(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double a1, double beta);
void metropolis_volume_next(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double &a1, double Pre, double beta);
void metropolis_volume(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double &a1, double Pre, double beta);
void metropolis_entropic_ads_next(int sigma[N], double position[N][2], int n_b[N][2*edge], 
	double alpha, double k, double mu, double a1, double &E, double e_max, double e_min, double bin_wid,
	int count_E[bin_E], double g_E[bin_E], double delta_g_E);
void metropolis_entropic_position_next(int sigma[N], double position[N][2], int n_b[N][2*edge], 
	double alpha, double k, double a1, double &E, double e_max, double e_min, double bin_wid,
	int count_E[bin_E], double g_E[bin_E], double delta_g_E);
void metropolis_entropic_volume_next(int sigma[N], double position[N][2], int n_b[N][2*edge], 
	double alpha, double k, double mu, double &a1, double Pre, 
	double &E, double e_max, double e_min, double bin_wid,
	int count_E[bin_E], double g_E[bin_E], double delta_g_E);
void metropolis_multi_ads_next(int sigma[N], double position[N][2], int n_b[N][2*edge], 
	double alpha, double k, double mu, double a1, double &E, double e_max, double e_min, double bin_wid,
	int count_E[bin_E], double g_E[bin_E]);
void metropolis_multi_position_next(int sigma[N], double position[N][2], int n_b[N][2*edge], 
	double alpha, double k, double a1, double &E, double e_max, double e_min, double bin_wid,
	int count_E[bin_E], double g_E[bin_E]);
void metropolis_multi_volume_next(int sigma[N], double position[N][2], int n_b[N][2*edge], 
	double alpha, double k, double mu, double &a1, double Pre, 
	double &E, double e_max, double e_min, double bin_wid,
	int count_E[bin_E], double g_E[bin_E]);