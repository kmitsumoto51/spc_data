#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include "mt19937ar.h"
#include "const.h"
#include "sim.h"
#include "hami.h"
#include "utils.h"

using namespace std;

void metropolis_ads_next(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double a1, double beta){
	int sigma_b;
	double delta_H;
	double r;
	double *posix;
	double *posiy;
	posix = new double[4];
	posiy = new double[4];
	double size_l = L*a1;
	double l_p_2;
	for(int i=0; i<N; i++){
		sigma_b = sigma[i];
		sigma[i] = (sigma_b+1)%2;
		posix[0] = position[i][0];
		posiy[0] = position[i][1];
		posiy[1] = position[n_b[i][0]][1];
		posix[3] = position[n_b[i][1]][0];
		if(i%L == L-1){
			posix[1] = position[n_b[i][0]][0] + size_l;
			posix[2] = position[n_b[i][4]][0] + size_l;
			if(i/L == L-1){
				posiy[2] = position[n_b[i][4]][1] + size_l;
				posiy[3] = position[n_b[i][1]][1] + size_l;
			}else{
				posiy[2] = position[n_b[i][4]][1];
				posiy[3] = position[n_b[i][1]][1];
			}
		}else{
			posix[1] = position[n_b[i][0]][0];
			posix[2] = position[n_b[i][4]][0];
			if(i/L == L-1){
				posiy[2] = position[n_b[i][4]][1] + size_l;
				posiy[3] = position[n_b[i][1]][1] + size_l;
			}else{
				posiy[2] = position[n_b[i][4]][1];
				posiy[3] = position[n_b[i][1]][1];
			}
		}
		l_p_2 = local_potential2_next(posix, posiy, n_b, alpha, k);
		delta_H = (sigma[i] - sigma_b)*(l_p_2 - mu);
		r=genrand_real2();
		if(r >= exp(-beta*delta_H)){
			sigma[i] = sigma_b;
		}
	}
	delete[] posix;
	delete[] posiy;
}

void metropolis_ads_next_exchange(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double a1, double beta){
	int sigma_b;
	double delta_H;
	double r;
	double *posix;
	double *posiy;
	posix = new double[4];
	posiy = new double[4];
	double *posix_j;
	double *posiy_j;
	posix_j = new double[4];
	posiy_j = new double[4];
	double size_l = L*a1;
	double l_p_2, l_p_2_j;
	int p_j;
	for(int i0=0; i0<N; i0++){
		int i = (int)(N*genrand_real2());
		p_j = n_b[i][(int)(edge*genrand_real2())];
		if(sigma[i] != sigma[p_j]){
			sigma_b = sigma[i];
			sigma[i] = sigma[p_j];
			sigma[p_j] = sigma_b;

			posix[0] = position[i][0];
			posiy[0] = position[i][1];
			posiy[1] = position[n_b[i][0]][1];
			posix[3] = position[n_b[i][1]][0];
			if(i%L == L-1){
				posix[1] = position[n_b[i][0]][0] + size_l;
				posix[2] = position[n_b[i][4]][0] + size_l;
				if(i/L == L-1){
					posiy[2] = position[n_b[i][4]][1] + size_l;
					posiy[3] = position[n_b[i][1]][1] + size_l;
				}else{
					posiy[2] = position[n_b[i][4]][1];
					posiy[3] = position[n_b[i][1]][1];
				}
			}else{
				posix[1] = position[n_b[i][0]][0];
				posix[2] = position[n_b[i][4]][0];
				if(i/L == L-1){
					posiy[2] = position[n_b[i][4]][1] + size_l;
					posiy[3] = position[n_b[i][1]][1] + size_l;
				}else{
					posiy[2] = position[n_b[i][4]][1];
					posiy[3] = position[n_b[i][1]][1];
				}
			}

			posix_j[0] = position[p_j][0];
			posiy_j[0] = position[p_j][1];
			posiy_j[1] = position[n_b[p_j][0]][1];
			posix_j[3] = position[n_b[p_j][1]][0];
			if(i%L == L-1){
				posix_j[1] = position[n_b[p_j][0]][0] + size_l;
				posix_j[2] = position[n_b[p_j][4]][0] + size_l;
				if(i/L == L-1){
					posiy_j[2] = position[n_b[p_j][4]][1] + size_l;
					posiy_j[3] = position[n_b[p_j][1]][1] + size_l;
				}else{
					posiy_j[2] = position[n_b[p_j][4]][1];
					posiy_j[3] = position[n_b[p_j][1]][1];
				}
			}else{
				posix_j[1] = position[n_b[p_j][0]][0];
				posix_j[2] = position[n_b[p_j][4]][0];
				if(i/L == L-1){
					posiy_j[2] = position[n_b[p_j][4]][1] + size_l;
					posiy_j[3] = position[n_b[p_j][1]][1] + size_l;
				}else{
					posiy_j[2] = position[n_b[p_j][4]][1];
					posiy_j[3] = position[n_b[p_j][1]][1];
				}
			}
			l_p_2 = local_potential2_next(posix, posiy, n_b, alpha, k);
			l_p_2_j = local_potential2_next(posix_j, posiy_j, n_b, alpha, k);
			delta_H = (sigma[i] - sigma[p_j])*(l_p_2 - l_p_2_j);
			r=genrand_real2();
			if(r >= exp(-beta*delta_H)){
				sigma[p_j] = sigma[i];
				sigma[i] = sigma_b;

			}
		}
	}
	delete[] posix;
	delete[] posiy;
	delete[] posix_j;
	delete[] posiy_j;
}

void metropolis_ads(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double a1, double beta){
	int sigma_b;
	double delta_H;
	double r;
	double *posix;
	double *posiy;
	posix = new double[4];
	posiy = new double[4];
	double size_l = L*a1;
	double l_p_2;
	for(int i=0; i<N; i++){
		sigma_b = sigma[i];
		sigma[i] = (sigma_b+1)%2;
		posix[0] = position[i][0];
		posiy[0] = position[i][1];
		posiy[1] = position[n_b[i][0]][1];
		posix[3] = position[n_b[i][1]][0];
		if(i%L == L-1){
			posix[1] = position[n_b[i][0]][0] + size_l;
			posix[2] = position[n_b[i][4]][0] + size_l;
			if(i/L == L-1){
				posiy[2] = position[n_b[i][4]][1] + size_l;
				posiy[3] = position[n_b[i][1]][1] + size_l;
			}else{
				posiy[2] = position[n_b[i][4]][1];
				posiy[3] = position[n_b[i][1]][1];
			}
		}else{
			posix[1] = position[n_b[i][0]][0];
			posix[2] = position[n_b[i][4]][0];
			if(i/L == L-1){
				posiy[2] = position[n_b[i][4]][1] + size_l;
				posiy[3] = position[n_b[i][1]][1] + size_l;
			}else{
				posiy[2] = position[n_b[i][4]][1];
				posiy[3] = position[n_b[i][1]][1];
			}
		}
		l_p_2 = local_potential2(posix, posiy, n_b, alpha, k);
		delta_H = (sigma[i] - sigma_b)*(l_p_2 - mu);
		r=genrand_real2();
		if(r >= exp(-beta*delta_H)){
			sigma[i] = sigma_b;
		}
	}
	delete[] posix;
	delete[] posiy;
}

void metropolis_position_next(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double a1, double beta){
	double l_d1,l_d2,l_d3,l_d4,l_d5,l_d6,l_d7,l_d8;
	double l_d1_new,l_d2_new,l_d3_new,l_d4_new,l_d5_new,l_d6_new,l_d7_new,l_d8_new;
	double radi;
	double theta;
	double posix_new, posiy_new;
	double delta_E;
	double r;
	double posix1, posix2, posix3, posix4, posix5, posix6, posix7, posix8;
	double posiy1, posiy2, posiy3, posiy4, posiy5, posiy6, posiy7, posiy8;
	double size_l = L*a1;
	double l_1 = 1+alpha;
	double l_2 = sqrt(2)*(1+alpha);
	for(int i=0; i<N; i++){
		radi = 0.1*sqrt(genrand_real2());
		theta = 2*M_PI*genrand_real2();
		posix_new = position[i][0] + radi*cos(theta);
		posiy_new = position[i][1] + radi*sin(theta);
		posiy1 = position[n_b[i][0]][1];
		posix2 = position[n_b[i][1]][0];
		posiy3 = position[n_b[i][2]][1];
		posix4 = position[n_b[i][3]][0];
		if(i%L == 0){
			posix1 = position[n_b[i][0]][0];
			posix3 = position[n_b[i][2]][0] - size_l;
			posix5 = position[n_b[i][4]][0];
			posix6 = position[n_b[i][5]][0] - size_l;
			posix7 = position[n_b[i][6]][0] - size_l;
			posix8 = position[n_b[i][7]][0];
			if(i/L == 0){
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1] - size_l;
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1] - size_l;
				posiy8 = position[n_b[i][7]][1] - size_l;
			}else if(i/L == L-1){
				posiy2 = position[n_b[i][1]][1] + size_l;
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1] + size_l;
				posiy6 = position[n_b[i][5]][1] + size_l;
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}else{
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}
		}else if(i%L == L-1){
			posix1 = position[n_b[i][0]][0] + size_l;
			posix3 = position[n_b[i][2]][0];
			posix5 = position[n_b[i][4]][0] + size_l;
			posix6 = position[n_b[i][5]][0];
			posix7 = position[n_b[i][6]][0];
			posix8 = position[n_b[i][7]][0] + size_l;
			if(i/L == 0){
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1] - size_l;
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1] - size_l;
				posiy8 = position[n_b[i][7]][1] - size_l;
			}else if(i/L == L-1){
				posiy2 = position[n_b[i][1]][1] + size_l;
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1] + size_l;
				posiy6 = position[n_b[i][5]][1] + size_l;
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}else{
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}
		}else{
			posix1 = position[n_b[i][0]][0];
			posix3 = position[n_b[i][2]][0];
			posix5 = position[n_b[i][4]][0];
			posix6 = position[n_b[i][5]][0];
			posix7 = position[n_b[i][6]][0];
			posix8 = position[n_b[i][7]][0];
			if(i/L == 0){
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1] - size_l;
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1] - size_l;
				posiy8 = position[n_b[i][7]][1] - size_l;
			}else if(i/L == L-1){
				posiy2 = position[n_b[i][1]][1] + size_l;
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1] + size_l;
				posiy6 = position[n_b[i][5]][1] + size_l;
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}else{
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}
		}
		int flg = 0;
		if(posix_new > posix1 || posix_new < posix3){
			flg = 1;
		}
		if(posix_new > posix5 || posix_new < posix6 || posix_new < posix7 || posix_new > posix8){
			flg = 1;
		}
		if(posiy_new > posiy2 || posiy_new < posiy4){
			flg = 1;
		}
		if(posiy_new > posiy5 || posiy_new > posiy6 || posiy_new < posiy7 || posiy_new < posiy8){
			flg = 1;
		}
		if(flg == 0){
			l_d1 = lattice_distance(position[i][0], position[i][1], posix1, posiy1);
			l_d2 = lattice_distance(position[i][0], position[i][1], posix2, posiy2);
			l_d3 = lattice_distance(position[i][0], position[i][1], posix3, posiy3);
			l_d4 = lattice_distance(position[i][0], position[i][1], posix4, posiy4);
			l_d5 = lattice_distance(position[i][0], position[i][1], posix5, posiy5);
			l_d6 = lattice_distance(position[i][0], position[i][1], posix6, posiy6);
			l_d7 = lattice_distance(position[i][0], position[i][1], posix7, posiy7);
			l_d8 = lattice_distance(position[i][0], position[i][1], posix8, posiy8);
			l_d1_new = lattice_distance(posix_new, posiy_new, posix1, posiy1);
			l_d2_new = lattice_distance(posix_new, posiy_new, posix2, posiy2);
			l_d3_new = lattice_distance(posix_new, posiy_new, posix3, posiy3);
			l_d4_new = lattice_distance(posix_new, posiy_new, posix4, posiy4);
			l_d5_new = lattice_distance(posix_new, posiy_new, posix5, posiy5);
			l_d6_new = lattice_distance(posix_new, posiy_new, posix6, posiy6);
			l_d7_new = lattice_distance(posix_new, posiy_new, posix7, posiy7);
			l_d8_new = lattice_distance(posix_new, posiy_new, posix8, posiy8);
			// cout << i << " " << l_d1 << " " << l_d2 << endl;
			delta_E = -(l_d1_new-l_d1) + 0.5*(l_d1_new*l_d1_new-l_d1*l_d1)
					  -(l_d2_new-l_d2) + 0.5*(l_d2_new*l_d2_new-l_d2*l_d2)
					  -(l_d3_new-l_d3) + 0.5*(l_d3_new*l_d3_new-l_d3*l_d3)
					  -(l_d4_new-l_d4) + 0.5*(l_d4_new*l_d4_new-l_d4*l_d4)
					  -sqrt(2)*(l_d5_new-l_d5) + 0.5*(l_d5_new*l_d5_new-l_d5*l_d5)
					  -sqrt(2)*(l_d6_new-l_d6) + 0.5*(l_d6_new*l_d6_new-l_d6*l_d6)
					  -sqrt(2)*(l_d7_new-l_d7) + 0.5*(l_d7_new*l_d7_new-l_d7*l_d7)
					  -sqrt(2)*(l_d8_new-l_d8) + 0.5*(l_d8_new*l_d8_new-l_d8*l_d8)
					  -sigma[i]*k*
					  (0.5*l_1*(l_d1_new-l_d1) - 0.25*(l_d1_new*l_d1_new-l_d1*l_d1)
					  +0.5*l_1*(l_d2_new-l_d2) - 0.25*(l_d2_new*l_d2_new-l_d2*l_d2)
					  +l_2*(l_d5_new-l_d5) - 0.5*(l_d5_new*l_d5_new-l_d5*l_d5))
					  -sigma[n_b[i][2]]*k*
					  (0.5*l_1*(l_d2_new-l_d2) - 0.25*(l_d2_new*l_d2_new-l_d2*l_d2)
					  +0.5*l_1*(l_d3_new-l_d3) - 0.25*(l_d3_new*l_d3_new-l_d3*l_d3)
					  +l_2*(l_d6_new-l_d6) - 0.5*(l_d6_new*l_d6_new-l_d6*l_d6))
					  -sigma[n_b[i][6]]*k*
					  (0.5*l_1*(l_d3_new-l_d3) - 0.25*(l_d3_new*l_d3_new-l_d3*l_d3)
					  +0.5*l_1*(l_d4_new-l_d4) - 0.25*(l_d4_new*l_d4_new-l_d4*l_d4)
					  +l_2*(l_d7_new-l_d7) - 0.5*(l_d7_new*l_d7_new-l_d7*l_d7))
					  -sigma[n_b[i][3]]*k*
					  (0.5*l_1*(l_d1_new-l_d1) - 0.25*(l_d1_new*l_d1_new-l_d1*l_d1)
					  +0.5*l_1*(l_d4_new-l_d4) - 0.25*(l_d4_new*l_d4_new-l_d4*l_d4)
					  +l_2*(l_d8_new-l_d8) - 0.5*(l_d8_new*l_d8_new-l_d8*l_d8));
			r=genrand_real2();
			if(r < exp(-beta*delta_E)){
				position[i][0] = posix_new;
				position[i][1] = posiy_new;
			}
		}
	}
	double posi_ori_x = position[0][0];
	double posi_ori_y = position[0][1];
	for(int i=0; i<N; i++){
		position[i][0] -= posi_ori_x;
		position[i][1] -= posi_ori_y;
	}
}

void metropolis_position(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double a1, double beta){
	double l_d1,l_d2,l_d3,l_d4;
	double l_d1_new,l_d2_new,l_d3_new,l_d4_new;
	double radi;
	double theta;
	double posix_new, posiy_new;
	double delta_E;
	double r;
	double posix1, posix2, posix3, posix4, posix5, posix6, posix7, posix8;
	double posiy1, posiy2, posiy3, posiy4, posiy5, posiy6, posiy7, posiy8;
	double size_l = L*a1;
	double l_1 = 1+alpha;
	for(int i=0; i<N; i++){
		radi = 0.1*sqrt(genrand_real2());
		theta = 2*M_PI*genrand_real2();
		posix_new = position[i][0] + radi*cos(theta);
		posiy_new = position[i][1] + radi*sin(theta);
		posiy1 = position[n_b[i][0]][1];
		posix2 = position[n_b[i][1]][0];
		posiy3 = position[n_b[i][2]][1];
		posix4 = position[n_b[i][3]][0];
		if(i%L == 0){
			posix1 = position[n_b[i][0]][0];
			posix3 = position[n_b[i][2]][0] - size_l;
			posix5 = position[n_b[i][4]][0];
			posix6 = position[n_b[i][5]][0] - size_l;
			posix7 = position[n_b[i][6]][0] - size_l;
			posix8 = position[n_b[i][7]][0];
			if(i/L == 0){
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1] - size_l;
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1] - size_l;
				posiy8 = position[n_b[i][7]][1] - size_l;
			}else if(i/L == L-1){
				posiy2 = position[n_b[i][1]][1] + size_l;
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1] + size_l;
				posiy6 = position[n_b[i][5]][1] + size_l;
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}else{
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}
		}else if(i%L == L-1){
			posix1 = position[n_b[i][0]][0] + size_l;
			posix3 = position[n_b[i][2]][0];
			posix5 = position[n_b[i][4]][0] + size_l;
			posix6 = position[n_b[i][5]][0];
			posix7 = position[n_b[i][6]][0];
			posix8 = position[n_b[i][7]][0] + size_l;
			if(i/L == 0){
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1] - size_l;
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1] - size_l;
				posiy8 = position[n_b[i][7]][1] - size_l;
			}else if(i/L == L-1){
				posiy2 = position[n_b[i][1]][1] + size_l;
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1] + size_l;
				posiy6 = position[n_b[i][5]][1] + size_l;
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}else{
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}
		}else{
			posix1 = position[n_b[i][0]][0];
			posix3 = position[n_b[i][2]][0];
			posix5 = position[n_b[i][4]][0];
			posix6 = position[n_b[i][5]][0];
			posix7 = position[n_b[i][6]][0];
			posix8 = position[n_b[i][7]][0];
			if(i/L == 0){
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1] - size_l;
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1] - size_l;
				posiy8 = position[n_b[i][7]][1] - size_l;
			}else if(i/L == L-1){
				posiy2 = position[n_b[i][1]][1] + size_l;
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1] + size_l;
				posiy6 = position[n_b[i][5]][1] + size_l;
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}else{
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}
		}
		int flg = 0;
		if(posix_new > posix1 || posix_new < posix3){
			flg = 1;
		}
		if(posix_new > posix5 || posix_new < posix6 || posix_new < posix7 || posix_new > posix8){
			flg = 1;
		}
		if(posiy_new > posiy2 || posiy_new < posiy4){
			flg = 1;
		}
		if(posiy_new > posiy5 || posiy_new > posiy6 || posiy_new < posiy7 || posiy_new < posiy8){
			flg = 1;
		}
		if(flg == 0){
			l_d1 = lattice_distance(position[i][0], position[i][1], posix1, posiy1);
			l_d2 = lattice_distance(position[i][0], position[i][1], posix2, posiy2);
			l_d3 = lattice_distance(position[i][0], position[i][1], posix3, posiy3);
			l_d4 = lattice_distance(position[i][0], position[i][1], posix4, posiy4);
			l_d1_new = lattice_distance(posix_new, posiy_new, posix1, posiy1);
			l_d2_new = lattice_distance(posix_new, posiy_new, posix2, posiy2);
			l_d3_new = lattice_distance(posix_new, posiy_new, posix3, posiy3);
			l_d4_new = lattice_distance(posix_new, posiy_new, posix4, posiy4);
			// cout << i << " " << l_d1 << " " << l_d2 << endl;
			delta_E = -(l_d1_new-l_d1) + 0.5*(l_d1_new*l_d1_new-l_d1*l_d1)
					  -(l_d2_new-l_d2) + 0.5*(l_d2_new*l_d2_new-l_d2*l_d2)
					  -(l_d3_new-l_d3) + 0.5*(l_d3_new*l_d3_new-l_d3*l_d3)
					  -(l_d4_new-l_d4) + 0.5*(l_d4_new*l_d4_new-l_d4*l_d4)
					  -sigma[i]*k*
					  (0.5*l_1*(l_d1_new-l_d1) - 0.25*(l_d1_new*l_d1_new-l_d1*l_d1)
					  +0.5*l_1*(l_d2_new-l_d2) - 0.25*(l_d2_new*l_d2_new-l_d2*l_d2))
					  -sigma[n_b[i][2]]*k*
					  (0.5*l_1*(l_d2_new-l_d2) - 0.25*(l_d2_new*l_d2_new-l_d2*l_d2)
					  +0.5*l_1*(l_d3_new-l_d3) - 0.25*(l_d3_new*l_d3_new-l_d3*l_d3))
					  -sigma[n_b[i][6]]*k*
					  (0.5*l_1*(l_d3_new-l_d3) - 0.25*(l_d3_new*l_d3_new-l_d3*l_d3)
					  +0.5*l_1*(l_d4_new-l_d4) - 0.25*(l_d4_new*l_d4_new-l_d4*l_d4))
					  -sigma[n_b[i][3]]*k*
					  (0.5*l_1*(l_d1_new-l_d1) - 0.25*(l_d1_new*l_d1_new-l_d1*l_d1)
					  +0.5*l_1*(l_d4_new-l_d4) - 0.25*(l_d4_new*l_d4_new-l_d4*l_d4));
			r=genrand_real2();
			if(r < exp(-beta*delta_E)){
				position[i][0] = posix_new;
				position[i][1] = posiy_new;
			}
		}
	}
	double posi_ori_x = position[0][0];
	double posi_ori_y = position[0][1];
	for(int i=0; i<N; i++){
		position[i][0] -= posi_ori_x;
		position[i][1] -= posi_ori_y;
	}
}

void metropolis_volume_next(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double &a1, double Pre, double beta){
	double H_b;
	double delta_H;
	double r = 2*genrand_real2()-1;
	double af_r = 1 + 0.01*r/a1;
	H_b = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);
	for(int j=0; j<N; j++){
		position[j][0] *= af_r;
		position[j][1] *= af_r;
	}
	a1 *= af_r;
	delta_H = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre) - H_b;
	r=genrand_real2();
	if(r >= exp(-beta*delta_H)){
		a1 /= af_r;
		for(int j=0; j<N; j++){
			position[j][0] /= af_r;
			position[j][1] /= af_r;
		}
	}
}

void metropolis_volume(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double &a1, double Pre, double beta){
	double H_b;
	double delta_H;
	double r = 2*genrand_real2()-1;
	double af_r = 1 + 0.01*r/a1;
	H_b = enthalpy(sigma, position, n_b, alpha, k, mu, a1, Pre);
	for(int j=0; j<N; j++){
		position[j][0] *= af_r;
		position[j][1] *= af_r;
	}
	a1 *= af_r;
	delta_H = enthalpy(sigma, position, n_b, alpha, k, mu, a1, Pre) - H_b;
	r=genrand_real2();
	if(r >= exp(-beta*delta_H)){
		a1 /= af_r;
		for(int j=0; j<N; j++){
			position[j][0] /= af_r;
			position[j][1] /= af_r;
		}
	}
}

void metropolis_entropic_ads_next(int sigma[N], double position[N][2], int n_b[N][2*edge], 
	double alpha, double k, double mu, double a1, double &E, double e_max, double e_min, double bin_wid,
	int count_E[bin_E], double g_E[bin_E], double delta_g_E){
	int sigma_b;
	double delta_H;
	double e_b;
	double e_a;
	int e_b_bin;
	int e_a_bin;
	double r;
	double *posix;
	double *posiy;
	posix = new double[4];
	posiy = new double[4];
	double size_l = L*a1;
	double l_p_2;
	for(int i=0; i<N; i++){
		sigma_b = sigma[i];
		sigma[i] = (sigma_b+1)%2;
		posix[0] = position[i][0];
		posiy[0] = position[i][1];
		posiy[1] = position[n_b[i][0]][1];
		posix[3] = position[n_b[i][1]][0];
		if(i%L == L-1){
			posix[1] = position[n_b[i][0]][0] + size_l;
			posix[2] = position[n_b[i][4]][0] + size_l;
			if(i/L == L-1){
				posiy[2] = position[n_b[i][4]][1] + size_l;
				posiy[3] = position[n_b[i][1]][1] + size_l;
			}else{
				posiy[2] = position[n_b[i][4]][1];
				posiy[3] = position[n_b[i][1]][1];
			}
		}else{
			posix[1] = position[n_b[i][0]][0];
			posix[2] = position[n_b[i][4]][0];
			if(i/L == L-1){
				posiy[2] = position[n_b[i][4]][1] + size_l;
				posiy[3] = position[n_b[i][1]][1] + size_l;
			}else{
				posiy[2] = position[n_b[i][4]][1];
				posiy[3] = position[n_b[i][1]][1];
			}
		}
		l_p_2 = local_potential2_next(posix, posiy, n_b, alpha, k);
		delta_H = (sigma[i] - sigma_b)*(l_p_2 - mu);
		e_b = E/(double)N;
		e_a = (E+delta_H)/(double)N;
		e_b_bin = (int)((e_b-e_min)/bin_wid);
		if(e_a > e_max){
			sigma[i] = sigma_b;
			count_E[e_b_bin] += 1;
			g_E[e_b_bin] += delta_g_E;
		}else{
			e_a_bin = (int)((e_a-e_min)/bin_wid);
			r=genrand_real2();
			if(r >= exp(-g_E[e_a_bin] + g_E[e_b_bin]) || e_a_bin < 100){
				sigma[i] = sigma_b;
				count_E[e_b_bin] += 1;
				g_E[e_b_bin] += delta_g_E;
			}else{
				count_E[e_a_bin] += 1;
				g_E[e_a_bin] += delta_g_E;
				E += delta_H;
			}
		}
	}
	delete[] posix;
	delete[] posiy;
}

void metropolis_entropic_position_next(int sigma[N], double position[N][2], int n_b[N][2*edge], 
	double alpha, double k, double a1, double &E, double e_max, double e_min, double bin_wid,
	int count_E[bin_E], double g_E[bin_E], double delta_g_E){
	double l_d1,l_d2,l_d3,l_d4,l_d5,l_d6,l_d7,l_d8;
	double l_d1_new,l_d2_new,l_d3_new,l_d4_new,l_d5_new,l_d6_new,l_d7_new,l_d8_new;
	double radi;
	double theta;
	double posix_new, posiy_new;
	double delta_E;
	double e_b;
	double e_a;
	int e_b_bin;
	int e_a_bin;
	double r;
	double posix1, posix2, posix3, posix4, posix5, posix6, posix7, posix8;
	double posiy1, posiy2, posiy3, posiy4, posiy5, posiy6, posiy7, posiy8;
	double size_l = L*a1;
	double l_1 = 1+alpha;
	double l_2 = sqrt(2)*(1+alpha);
	for(int i=0; i<N; i++){
		radi = 0.1*sqrt(genrand_real2());
		theta = 2*M_PI*genrand_real2();
		posix_new = position[i][0] + radi*cos(theta);
		posiy_new = position[i][1] + radi*sin(theta);
		posiy1 = position[n_b[i][0]][1];
		posix2 = position[n_b[i][1]][0];
		posiy3 = position[n_b[i][2]][1];
		posix4 = position[n_b[i][3]][0];
		if(i%L == 0){
			posix1 = position[n_b[i][0]][0];
			posix3 = position[n_b[i][2]][0] - size_l;
			posix5 = position[n_b[i][4]][0];
			posix6 = position[n_b[i][5]][0] - size_l;
			posix7 = position[n_b[i][6]][0] - size_l;
			posix8 = position[n_b[i][7]][0];
			if(i/L == 0){
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1] - size_l;
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1] - size_l;
				posiy8 = position[n_b[i][7]][1] - size_l;
			}else if(i/L == L-1){
				posiy2 = position[n_b[i][1]][1] + size_l;
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1] + size_l;
				posiy6 = position[n_b[i][5]][1] + size_l;
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}else{
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}
		}else if(i%L == L-1){
			posix1 = position[n_b[i][0]][0] + size_l;
			posix3 = position[n_b[i][2]][0];
			posix5 = position[n_b[i][4]][0] + size_l;
			posix6 = position[n_b[i][5]][0];
			posix7 = position[n_b[i][6]][0];
			posix8 = position[n_b[i][7]][0] + size_l;
			if(i/L == 0){
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1] - size_l;
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1] - size_l;
				posiy8 = position[n_b[i][7]][1] - size_l;
			}else if(i/L == L-1){
				posiy2 = position[n_b[i][1]][1] + size_l;
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1] + size_l;
				posiy6 = position[n_b[i][5]][1] + size_l;
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}else{
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}
		}else{
			posix1 = position[n_b[i][0]][0];
			posix3 = position[n_b[i][2]][0];
			posix5 = position[n_b[i][4]][0];
			posix6 = position[n_b[i][5]][0];
			posix7 = position[n_b[i][6]][0];
			posix8 = position[n_b[i][7]][0];
			if(i/L == 0){
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1] - size_l;
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1] - size_l;
				posiy8 = position[n_b[i][7]][1] - size_l;
			}else if(i/L == L-1){
				posiy2 = position[n_b[i][1]][1] + size_l;
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1] + size_l;
				posiy6 = position[n_b[i][5]][1] + size_l;
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}else{
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}
		}
		int flg = 0;
		if(posix_new > posix1 || posix_new < posix3){
			flg = 1;
		}
		if(posix_new > posix5 || posix_new < posix6 || posix_new < posix7 || posix_new > posix8){
			flg = 1;
		}
		if(posiy_new > posiy2 || posiy_new < posiy4){
			flg = 1;
		}
		if(posiy_new > posiy5 || posiy_new > posiy6 || posiy_new < posiy7 || posiy_new < posiy8){
			flg = 1;
		}
		if(flg == 0){
			l_d1 = lattice_distance(position[i][0], position[i][1], posix1, posiy1);
			l_d2 = lattice_distance(position[i][0], position[i][1], posix2, posiy2);
			l_d3 = lattice_distance(position[i][0], position[i][1], posix3, posiy3);
			l_d4 = lattice_distance(position[i][0], position[i][1], posix4, posiy4);
			l_d5 = lattice_distance(position[i][0], position[i][1], posix5, posiy5);
			l_d6 = lattice_distance(position[i][0], position[i][1], posix6, posiy6);
			l_d7 = lattice_distance(position[i][0], position[i][1], posix7, posiy7);
			l_d8 = lattice_distance(position[i][0], position[i][1], posix8, posiy8);
			l_d1_new = lattice_distance(posix_new, posiy_new, posix1, posiy1);
			l_d2_new = lattice_distance(posix_new, posiy_new, posix2, posiy2);
			l_d3_new = lattice_distance(posix_new, posiy_new, posix3, posiy3);
			l_d4_new = lattice_distance(posix_new, posiy_new, posix4, posiy4);
			l_d5_new = lattice_distance(posix_new, posiy_new, posix5, posiy5);
			l_d6_new = lattice_distance(posix_new, posiy_new, posix6, posiy6);
			l_d7_new = lattice_distance(posix_new, posiy_new, posix7, posiy7);
			l_d8_new = lattice_distance(posix_new, posiy_new, posix8, posiy8);
			// cout << i << " " << l_d1 << " " << l_d2 << endl;
			delta_E = -(l_d1_new-l_d1) + 0.5*(l_d1_new*l_d1_new-l_d1*l_d1)
					  -(l_d2_new-l_d2) + 0.5*(l_d2_new*l_d2_new-l_d2*l_d2)
					  -(l_d3_new-l_d3) + 0.5*(l_d3_new*l_d3_new-l_d3*l_d3)
					  -(l_d4_new-l_d4) + 0.5*(l_d4_new*l_d4_new-l_d4*l_d4)
					  -sqrt(2)*(l_d5_new-l_d5) + 0.5*(l_d5_new*l_d5_new-l_d5*l_d5)
					  -sqrt(2)*(l_d6_new-l_d6) + 0.5*(l_d6_new*l_d6_new-l_d6*l_d6)
					  -sqrt(2)*(l_d7_new-l_d7) + 0.5*(l_d7_new*l_d7_new-l_d7*l_d7)
					  -sqrt(2)*(l_d8_new-l_d8) + 0.5*(l_d8_new*l_d8_new-l_d8*l_d8)
					  -sigma[i]*k*
					  (0.5*l_1*(l_d1_new-l_d1) - 0.25*(l_d1_new*l_d1_new-l_d1*l_d1)
					  +0.5*l_1*(l_d2_new-l_d2) - 0.25*(l_d2_new*l_d2_new-l_d2*l_d2)
					  +l_2*(l_d5_new-l_d5) - 0.5*(l_d5_new*l_d5_new-l_d5*l_d5))
					  -sigma[n_b[i][2]]*k*
					  (0.5*l_1*(l_d2_new-l_d2) - 0.25*(l_d2_new*l_d2_new-l_d2*l_d2)
					  +0.5*l_1*(l_d3_new-l_d3) - 0.25*(l_d3_new*l_d3_new-l_d3*l_d3)
					  +l_2*(l_d6_new-l_d6) - 0.5*(l_d6_new*l_d6_new-l_d6*l_d6))
					  -sigma[n_b[i][6]]*k*
					  (0.5*l_1*(l_d3_new-l_d3) - 0.25*(l_d3_new*l_d3_new-l_d3*l_d3)
					  +0.5*l_1*(l_d4_new-l_d4) - 0.25*(l_d4_new*l_d4_new-l_d4*l_d4)
					  +l_2*(l_d7_new-l_d7) - 0.5*(l_d7_new*l_d7_new-l_d7*l_d7))
					  -sigma[n_b[i][3]]*k*
					  (0.5*l_1*(l_d1_new-l_d1) - 0.25*(l_d1_new*l_d1_new-l_d1*l_d1)
					  +0.5*l_1*(l_d4_new-l_d4) - 0.25*(l_d4_new*l_d4_new-l_d4*l_d4)
					  +l_2*(l_d8_new-l_d8) - 0.5*(l_d8_new*l_d8_new-l_d8*l_d8));
			e_b = E/(double)N;
			e_a = (E+delta_E)/(double)N;
			e_b_bin = (int)((e_b-e_min)/bin_wid);
			if(e_a > e_max){
				count_E[e_b_bin] += 1;
				g_E[e_b_bin] += delta_g_E;
			}else{
				e_a_bin = (int)((e_a-e_min)/bin_wid);
				r=genrand_real2();
				if(r >= exp(-g_E[e_a_bin] + g_E[e_b_bin]) || e_a_bin < 100){
					count_E[e_b_bin] += 1;
					g_E[e_b_bin] += delta_g_E;
				}else{
					position[i][0] = posix_new;
					position[i][1] = posiy_new;
					count_E[e_a_bin] += 1;
					g_E[e_a_bin] += delta_g_E;
					E += delta_E;
				}
			}
		}
	}
	double posi_ori_x = position[0][0];
	double posi_ori_y = position[0][1];
	for(int i=0; i<N; i++){
		position[i][0] -= posi_ori_x;
		position[i][1] -= posi_ori_y;
	}
}

void metropolis_entropic_volume_next(int sigma[N], double position[N][2], int n_b[N][2*edge], 
	double alpha, double k, double mu, double &a1, double Pre, 
	double &E, double e_max, double e_min, double bin_wid,
	int count_E[bin_E], double g_E[bin_E], double delta_g_E){
	double H_b;
	double e_b;
	double e_a;
	int e_b_bin;
	int e_a_bin;
	double delta_H;
	double r = 2*genrand_real2()-1;
	double af_r = 1 + 0.01*r/a1;
	H_b = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);
	for(int j=0; j<N; j++){
		position[j][0] *= af_r;
		position[j][1] *= af_r;
	}
	a1 *= af_r;
	delta_H = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre) - H_b;
	e_b = E/(double)N;
	e_a = (E+delta_H)/(double)N;
	e_b_bin = (int)((e_b-e_min)/bin_wid);
	if(e_a > e_max){
		count_E[e_b_bin] += 1;
		g_E[e_b_bin] += delta_g_E;
		a1 /= af_r;
		for(int j=0; j<N; j++){
			position[j][0] /= af_r;
			position[j][1] /= af_r;
		}
	}else{
		e_a_bin = (int)((e_a-e_min)/bin_wid);
		r=genrand_real2();
		if(r >= exp(-g_E[e_a_bin] + g_E[e_b_bin]) || e_a_bin < 100){
			count_E[e_b_bin] += 1;
			g_E[e_b_bin] += delta_g_E;
			a1 /= af_r;
			for(int j=0; j<N; j++){
				position[j][0] /= af_r;
				position[j][1] /= af_r;
			}
		}else{
			count_E[e_a_bin] += 1;
			g_E[e_a_bin] += delta_g_E;
			E += delta_H;
		}
	}
}

void metropolis_multi_ads_next(int sigma[N], double position[N][2], int n_b[N][2*edge], 
	double alpha, double k, double mu, double a1, double &E, double e_max, double e_min, double bin_wid,
	int count_E[bin_E], double g_E[bin_E]){
	int sigma_b;
	double delta_H;
	double e_b;
	double e_a;
	int e_b_bin;
	int e_a_bin;
	double r;
	double *posix;
	double *posiy;
	posix = new double[4];
	posiy = new double[4];
	double size_l = L*a1;
	double l_p_2;
	for(int i=0; i<N; i++){
		sigma_b = sigma[i];
		sigma[i] = (sigma_b+1)%2;
		posix[0] = position[i][0];
		posiy[0] = position[i][1];
		posiy[1] = position[n_b[i][0]][1];
		posix[3] = position[n_b[i][1]][0];
		if(i%L == L-1){
			posix[1] = position[n_b[i][0]][0] + size_l;
			posix[2] = position[n_b[i][4]][0] + size_l;
			if(i/L == L-1){
				posiy[2] = position[n_b[i][4]][1] + size_l;
				posiy[3] = position[n_b[i][1]][1] + size_l;
			}else{
				posiy[2] = position[n_b[i][4]][1];
				posiy[3] = position[n_b[i][1]][1];
			}
		}else{
			posix[1] = position[n_b[i][0]][0];
			posix[2] = position[n_b[i][4]][0];
			if(i/L == L-1){
				posiy[2] = position[n_b[i][4]][1] + size_l;
				posiy[3] = position[n_b[i][1]][1] + size_l;
			}else{
				posiy[2] = position[n_b[i][4]][1];
				posiy[3] = position[n_b[i][1]][1];
			}
		}
		l_p_2 = local_potential2_next(posix, posiy, n_b, alpha, k);
		delta_H = (sigma[i] - sigma_b)*(l_p_2 - mu);
		e_b = E/(double)N;
		e_a = (E+delta_H)/(double)N;
		e_b_bin = (int)((e_b-e_min)/bin_wid);
		if(e_a >= e_max){
			sigma[i] = sigma_b;
			count_E[e_b_bin] += 1;
		}else{
			e_a_bin = (int)((e_a-e_min)/bin_wid);
			r=genrand_real2();
			if(r >= exp(-g_E[e_a_bin] + g_E[e_b_bin]) || e_a_bin < 100){
				sigma[i] = sigma_b;
				count_E[e_b_bin] += 1;
			}else{
				count_E[e_a_bin] += 1;
				E += delta_H;
			}
		}
	}
	delete[] posix;
	delete[] posiy;
}

void metropolis_multi_position_next(int sigma[N], double position[N][2], int n_b[N][2*edge], 
	double alpha, double k, double a1, double &E, double e_max, double e_min, double bin_wid,
	int count_E[bin_E], double g_E[bin_E]){
	double l_d1,l_d2,l_d3,l_d4,l_d5,l_d6,l_d7,l_d8;
	double l_d1_new,l_d2_new,l_d3_new,l_d4_new,l_d5_new,l_d6_new,l_d7_new,l_d8_new;
	double radi;
	double theta;
	double posix_new, posiy_new;
	double delta_E;
	double e_b;
	double e_a;
	int e_b_bin;
	int e_a_bin;
	double r;
	double posix1, posix2, posix3, posix4, posix5, posix6, posix7, posix8;
	double posiy1, posiy2, posiy3, posiy4, posiy5, posiy6, posiy7, posiy8;
	double size_l = L*a1;
	double l_1 = 1+alpha;
	double l_2 = sqrt(2)*(1+alpha);
	for(int i=0; i<N; i++){
		radi = 0.1*sqrt(genrand_real2());
		theta = 2*M_PI*genrand_real2();
		posix_new = position[i][0] + radi*cos(theta);
		posiy_new = position[i][1] + radi*sin(theta);
		posiy1 = position[n_b[i][0]][1];
		posix2 = position[n_b[i][1]][0];
		posiy3 = position[n_b[i][2]][1];
		posix4 = position[n_b[i][3]][0];
		if(i%L == 0){
			posix1 = position[n_b[i][0]][0];
			posix3 = position[n_b[i][2]][0] - size_l;
			posix5 = position[n_b[i][4]][0];
			posix6 = position[n_b[i][5]][0] - size_l;
			posix7 = position[n_b[i][6]][0] - size_l;
			posix8 = position[n_b[i][7]][0];
			if(i/L == 0){
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1] - size_l;
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1] - size_l;
				posiy8 = position[n_b[i][7]][1] - size_l;
			}else if(i/L == L-1){
				posiy2 = position[n_b[i][1]][1] + size_l;
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1] + size_l;
				posiy6 = position[n_b[i][5]][1] + size_l;
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}else{
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}
		}else if(i%L == L-1){
			posix1 = position[n_b[i][0]][0] + size_l;
			posix3 = position[n_b[i][2]][0];
			posix5 = position[n_b[i][4]][0] + size_l;
			posix6 = position[n_b[i][5]][0];
			posix7 = position[n_b[i][6]][0];
			posix8 = position[n_b[i][7]][0] + size_l;
			if(i/L == 0){
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1] - size_l;
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1] - size_l;
				posiy8 = position[n_b[i][7]][1] - size_l;
			}else if(i/L == L-1){
				posiy2 = position[n_b[i][1]][1] + size_l;
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1] + size_l;
				posiy6 = position[n_b[i][5]][1] + size_l;
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}else{
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}
		}else{
			posix1 = position[n_b[i][0]][0];
			posix3 = position[n_b[i][2]][0];
			posix5 = position[n_b[i][4]][0];
			posix6 = position[n_b[i][5]][0];
			posix7 = position[n_b[i][6]][0];
			posix8 = position[n_b[i][7]][0];
			if(i/L == 0){
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1] - size_l;
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1] - size_l;
				posiy8 = position[n_b[i][7]][1] - size_l;
			}else if(i/L == L-1){
				posiy2 = position[n_b[i][1]][1] + size_l;
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1] + size_l;
				posiy6 = position[n_b[i][5]][1] + size_l;
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}else{
				posiy2 = position[n_b[i][1]][1];
				posiy4 = position[n_b[i][3]][1];
				posiy5 = position[n_b[i][4]][1];
				posiy6 = position[n_b[i][5]][1];
				posiy7 = position[n_b[i][6]][1];
				posiy8 = position[n_b[i][7]][1];
			}
		}
		int flg = 0;
		if(posix_new > posix1 || posix_new < posix3){
			flg = 1;
		}
		if(posix_new > posix5 || posix_new < posix6 || posix_new < posix7 || posix_new > posix8){
			flg = 1;
		}
		if(posiy_new > posiy2 || posiy_new < posiy4){
			flg = 1;
		}
		if(posiy_new > posiy5 || posiy_new > posiy6 || posiy_new < posiy7 || posiy_new < posiy8){
			flg = 1;
		}
		if(flg == 0){
			l_d1 = lattice_distance(position[i][0], position[i][1], posix1, posiy1);
			l_d2 = lattice_distance(position[i][0], position[i][1], posix2, posiy2);
			l_d3 = lattice_distance(position[i][0], position[i][1], posix3, posiy3);
			l_d4 = lattice_distance(position[i][0], position[i][1], posix4, posiy4);
			l_d5 = lattice_distance(position[i][0], position[i][1], posix5, posiy5);
			l_d6 = lattice_distance(position[i][0], position[i][1], posix6, posiy6);
			l_d7 = lattice_distance(position[i][0], position[i][1], posix7, posiy7);
			l_d8 = lattice_distance(position[i][0], position[i][1], posix8, posiy8);
			l_d1_new = lattice_distance(posix_new, posiy_new, posix1, posiy1);
			l_d2_new = lattice_distance(posix_new, posiy_new, posix2, posiy2);
			l_d3_new = lattice_distance(posix_new, posiy_new, posix3, posiy3);
			l_d4_new = lattice_distance(posix_new, posiy_new, posix4, posiy4);
			l_d5_new = lattice_distance(posix_new, posiy_new, posix5, posiy5);
			l_d6_new = lattice_distance(posix_new, posiy_new, posix6, posiy6);
			l_d7_new = lattice_distance(posix_new, posiy_new, posix7, posiy7);
			l_d8_new = lattice_distance(posix_new, posiy_new, posix8, posiy8);
			// cout << i << " " << l_d1 << " " << l_d2 << endl;
			delta_E = -(l_d1_new-l_d1) + 0.5*(l_d1_new*l_d1_new-l_d1*l_d1)
					  -(l_d2_new-l_d2) + 0.5*(l_d2_new*l_d2_new-l_d2*l_d2)
					  -(l_d3_new-l_d3) + 0.5*(l_d3_new*l_d3_new-l_d3*l_d3)
					  -(l_d4_new-l_d4) + 0.5*(l_d4_new*l_d4_new-l_d4*l_d4)
					  -sqrt(2)*(l_d5_new-l_d5) + 0.5*(l_d5_new*l_d5_new-l_d5*l_d5)
					  -sqrt(2)*(l_d6_new-l_d6) + 0.5*(l_d6_new*l_d6_new-l_d6*l_d6)
					  -sqrt(2)*(l_d7_new-l_d7) + 0.5*(l_d7_new*l_d7_new-l_d7*l_d7)
					  -sqrt(2)*(l_d8_new-l_d8) + 0.5*(l_d8_new*l_d8_new-l_d8*l_d8)
					  -sigma[i]*k*
					  (0.5*l_1*(l_d1_new-l_d1) - 0.25*(l_d1_new*l_d1_new-l_d1*l_d1)
					  +0.5*l_1*(l_d2_new-l_d2) - 0.25*(l_d2_new*l_d2_new-l_d2*l_d2)
					  +l_2*(l_d5_new-l_d5) - 0.5*(l_d5_new*l_d5_new-l_d5*l_d5))
					  -sigma[n_b[i][2]]*k*
					  (0.5*l_1*(l_d2_new-l_d2) - 0.25*(l_d2_new*l_d2_new-l_d2*l_d2)
					  +0.5*l_1*(l_d3_new-l_d3) - 0.25*(l_d3_new*l_d3_new-l_d3*l_d3)
					  +l_2*(l_d6_new-l_d6) - 0.5*(l_d6_new*l_d6_new-l_d6*l_d6))
					  -sigma[n_b[i][6]]*k*
					  (0.5*l_1*(l_d3_new-l_d3) - 0.25*(l_d3_new*l_d3_new-l_d3*l_d3)
					  +0.5*l_1*(l_d4_new-l_d4) - 0.25*(l_d4_new*l_d4_new-l_d4*l_d4)
					  +l_2*(l_d7_new-l_d7) - 0.5*(l_d7_new*l_d7_new-l_d7*l_d7))
					  -sigma[n_b[i][3]]*k*
					  (0.5*l_1*(l_d1_new-l_d1) - 0.25*(l_d1_new*l_d1_new-l_d1*l_d1)
					  +0.5*l_1*(l_d4_new-l_d4) - 0.25*(l_d4_new*l_d4_new-l_d4*l_d4)
					  +l_2*(l_d8_new-l_d8) - 0.5*(l_d8_new*l_d8_new-l_d8*l_d8));
			e_b = E/(double)N;
			e_a = (E+delta_E)/(double)N;
			e_b_bin = (int)((e_b-e_min)/bin_wid);
			if(e_a >= e_max){
				count_E[e_b_bin] += 1;
			}else{
				e_a_bin = (int)((e_a-e_min)/bin_wid);
				r=genrand_real2();
				if(r >= exp(-g_E[e_a_bin] + g_E[e_b_bin]) || e_a_bin < 100){
					count_E[e_b_bin] += 1;
				}else{
					position[i][0] = posix_new;
					position[i][1] = posiy_new;
					count_E[e_a_bin] += 1;
					E += delta_E;
				}
			}
		}
	}
	double posi_ori_x = position[0][0];
	double posi_ori_y = position[0][1];
	for(int i=0; i<N; i++){
		position[i][0] -= posi_ori_x;
		position[i][1] -= posi_ori_y;
	}
}

void metropolis_multi_volume_next(int sigma[N], double position[N][2], int n_b[N][2*edge], 
	double alpha, double k, double mu, double &a1, double Pre, 
	double &E, double e_max, double e_min, double bin_wid,
	int count_E[bin_E], double g_E[bin_E]){
	double H_b;
	double e_b;
	double e_a;
	int e_b_bin;
	int e_a_bin;
	double delta_H;
	double r = 2*genrand_real2()-1;
	double af_r = 1 + 0.01*r/a1;
	H_b = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);
	for(int j=0; j<N; j++){
		position[j][0] *= af_r;
		position[j][1] *= af_r;
	}
	a1 *= af_r;
	delta_H = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre) - H_b;
	e_b = E/(double)N;
	e_a = (E+delta_H)/(double)N;
	e_b_bin = (int)((e_b-e_min)/bin_wid);
	if(e_a >= e_max){
		count_E[e_b_bin] += 1;
		a1 /= af_r;
		for(int j=0; j<N; j++){
			position[j][0] /= af_r;
			position[j][1] /= af_r;
		}
	}else{
		e_a_bin = (int)((e_a-e_min)/bin_wid);
		r=genrand_real2();
		if(r >= exp(-g_E[e_a_bin] + g_E[e_b_bin]) || e_a_bin < 100){
			count_E[e_b_bin] += 1;
			a1 /= af_r;
			for(int j=0; j<N; j++){
				position[j][0] /= af_r;
				position[j][1] /= af_r;
			}
		}else{
			count_E[e_a_bin] += 1;
			E += delta_H;
		}
	}
}