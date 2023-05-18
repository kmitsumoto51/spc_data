#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>
#include "const.h"
#include "mt19937ar.h"
#include "utils.h"
#include "hami.h"
#include "sim.h"

using namespace std;

const int num_phys = 6;

int count_E[bin_E];
int count_E2[bin_E];
double g_E[bin_E];
double physical_q[num_phys][bin_E];
int n_b[N][2*edge];
int sigma[N];
double position[N][2];
double a1;
int mu_in;
int Pre_in;
int alpha_in;
int k_in;
int seed;

void initial_conf(int sigma[N], double position[N][2], double alpha, double k, double &a1){
	int sum_sig=0;
	for(int i=0; i<N; i++){
		// sigma[i] = 1;
		// sum_sig += 1;
		double r = genrand_real2();
		if(r > 0.5){
			sigma[i] = 1;
			sum_sig += 1;
		}else{
			sigma[i] = 0;
		}
	}
	a1 = 1 + alpha*k*sum_sig/((1+k)*N);
	for(int i=0; i<N; i++){
		position[i][0] = (i%L)*a1;
			position[i][1] = (i/L)*a1;
	}
}

int main(int argc, char **argv){
	if(argc > 1){
		mu_in = atoi(argv[1]);
		Pre_in = atoi(argv[2]);
		alpha_in = atoi(argv[3]);
		k_in = atoi(argv[4]);
		seed = atoi(argv[5]);
	}else{
		return 1;
	}
	double mu = mu_in/100.0;
	double alpha;
	if(alpha_in == 67){
		alpha = 2.0/3.0;
	}else{
		alpha = alpha_in/100.0;
	}
	double Pre = Pre_in/100.0;
	double k = k_in/100.0;
	make_network(n_b);
	init_genrand(seed);
	initial_conf(sigma, position, alpha,k,a1);
	double e_max = 1.2;
	double e_min = -mu + 3*k*alpha*alpha/(1+k);
	double bin_wid = (e_max-e_min)/(double)bin_E;
	for(int i=0; i<bin_E; i++){
		g_E[i] = 0.0;
	}
	for(int i=0; i<bin_E; i++){
		count_E[i] = 0;
	}
	fstream f1;
	stringstream name1;
	name1 << "eq_physical" << "/" << L << "_" << mu_in << "_" << Pre_in 
	<< "_" << alpha_in << "_" << k_in << "_" << step1 << "_" << 1 << "_w_l.txt";
	f1.open(name1.str().c_str(), ios::in);
	for(int i=100; i<bin_E; i++){
		double aa, bb; 
		int cc, dd, ee;
		f1 >> cc >> aa >> bb >> dd >> ee >> g_E[i];
	}
	f1.close();
	double E, e_t;
	int Nads;
	double volu;
	E = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);
	Nads = num_ads(sigma);
	cout << Nads/(double)N << " " << E/(double)N << endl;
	for(int t=0; t<step1*500; t++){
		metropolis_multi_ads_next(sigma, position, n_b, alpha, k, mu, a1, E, e_max, e_min, bin_wid, count_E, g_E);
		// cout << "ads OK" << endl;
		for(int l=0; l<L; l++){
			metropolis_multi_position_next(sigma, position, n_b, alpha, k, a1, E, e_max, e_min, bin_wid, count_E, g_E);
			metropolis_multi_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, E, e_max, e_min, bin_wid, count_E, g_E);
			// cout << "vol OK" << endl;
		}
	}
	for(int i=100; i<bin_E; i++){
		if(count_E[i] != 0){
			g_E[i] += log(count_E[i]);
		}
	}
	for(int i=100; i<bin_E; i++){
		count_E[i] = 0;
		count_E2[i] = 0;
	}
	for(int i=0; i<num_phys; i++){
			for(int j=0; j<bin_E; j++){
					physical_q[i][j] = 0.0;
			}
	}
	for(int t=0; t<step1*500; t++){
		metropolis_multi_ads_next(sigma, position, n_b, alpha, k, mu, a1, 
		E, e_max, e_min, bin_wid, count_E, g_E);
		for(int l=0; l<L; l++){
			metropolis_multi_position_next(sigma, position, n_b, alpha, k, a1, 
				E, e_max, e_min, bin_wid, count_E, g_E);
			metropolis_multi_volume_next(sigma, position, n_b, alpha, k, mu, a1, 
				Pre, E, e_max, e_min, bin_wid, count_E, g_E);
		}
		Nads = num_ads(sigma);
		volu = a1*a1*L*L;
		e_t = E/(double)N;
		for(int j=100; j<bin_E; j++){
				if(e_t > e_min + bin_wid*j && e_t < e_min + bin_wid*(j+1)){
					count_E2[j] += 1;
					physical_q[0][j] += E;
					physical_q[1][j] += E*E;
					physical_q[2][j] += Nads;
					physical_q[3][j] += Nads*Nads;
					physical_q[4][j] += volu;
					physical_q[5][j] += volu*volu;
				}
			}
		}
		for(int j=100; j<bin_E; j++){
			if(count_E2[j] ==0){
				for(int i=0; i<num_phys; i++){
						physical_q[i][j] = 0.0;
				}
			}else{
				for(int i=0; i<num_phys; i++){
						physical_q[i][j] /= count_E2[j];
				}
			}
		}
		fstream f2;
		stringstream name2;
		name2 << "eq_physical" << "/" << L 
		<< "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_" << step1 << "_" << seed << "_multi.txt";
		f2.open(name2.str().c_str(), ios::out);
		for(int j=100; j<bin_E; j++){
			f2 << j << " " //1
				 << e_min + bin_wid*(j+0.5) << " " //2
				 << count_E[j] << " " //3
				 << g_E[j] << " " //4
				 << physical_q[0][j] << " " //5 E
				 << physical_q[1][j] << " " //6 E_sq
				 << physical_q[2][j] << " " //7 N_ads
				 << physical_q[3][j] << " " //8 N_ads*N_ads
				 << physical_q[4][j] << " " //9 vol
				 << physical_q[5][j] << endl; //10 vol_sq
		}
		f2.close();
}
