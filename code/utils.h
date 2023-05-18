#pragma once
#include "const.h"

void make_network(int n_b[N][2*edge]);
int num_ads(int sigma[N]);
double lattice_distance(double x1, double y1, double x2, double y2);
double susceptibility(double m_sq_ave,double m_ave_sq,double beta);
double specific_heat(double E_sq_ave,double E_ave,double beta);
void local_potential(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu,
	double a1, double p_plaquette[N], double p_plaquette2[N], double p_plaquette3[N], double f_bond[N][4]);
double local_potential2(double posix[4], double posiy[4], int n_b[N][2*edge], double alpha, double k);
double local_potential2_next(double posix[4], double posiy[4], int n_b[N][2*edge], double alpha, double k);
void force_average(double f_bond[N][4], double f_ave[2]);
void force_correlation(double f_bond[N][4], double f_correlation[L/2][4]);
void particle_correlation(int sigma[N], double ads_correlation[L/2]);
void plaquette_volume(double position[N][2], int n_b[N][2*edge], double a1, double p_vol[N]);
void contact_list(int sigma[N], int n_b[N][2*edge], int contact_f[N][edge], int contact_v[N][edge]);
int site_percolation_f(int sigma[N], int n_b[N][2*edge], int cluster_num_f[N], int &max_cluster_size_f, int contact_f[N][edge]);
void cluster_ratio_f(int sigma[N], int n_b[N][2*edge], int cluster_num_f[N], double c_ratio_f[2], int contact_f[N][edge]);
int site_percolation_v(int sigma[N], int n_b[N][2*edge], int cluster_num_v[N], int &max_cluster_size_v, int contact_v[N][edge]);
void cluster_ratio_v(int sigma[N], int n_b[N][2*edge], int cluster_num_v[N], double c_ratio_v[2], int contact_v[N][edge]);
void cluster_ratio_hist_f(int sigma[N], int n_b[N][2*edge], int cluster_num_f[N], int contact_f[N][edge], int size_hist_f[N-1], double ratio_hist_f[N-1]);
void cluster_ratio_hist_v(int sigma[N], int n_b[N][2*edge], int cluster_num_v[N], int contact_v[N][edge], int size_hist_v[N-1], double ratio_hist_v[N-1]);