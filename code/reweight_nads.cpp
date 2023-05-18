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
const int bin_n = 100;

double e_in[bin_E];
int count_E[bin_E];
double count_E_re[bin_E];
int count_E_n[bin_E][bin_n];
double count_n_re[bin_n];
double energy_ave[bin_n];
double g_E[bin_E];
double physical_q[num_phys][bin_E];
double physical_q_ave[num_phys];
int mu_in;
int Pre_in;
int alpha_in;
int k_in;
int seed;

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
	double alpha = alpha_in/100.0;
    double Pre = Pre_in/100.0;
    double k = k_in/100.0;
	double E_ave, E_sq_ave,C;
	int Nads;
	double Nads_ave, Nads_sq_ave, Nads_ave_sq, chi;
	double volu;
	double volume_ave, volume_sq_ave, volume_ave_sq, kappa;
	double T, beta;
	double e_max = 1.0;
  	double e_min = -mu + 3*k*alpha*alpha/(1+k);
  	double bin_wid = (e_max-e_min)/(double)bin_E;
  	for(int j=0; j<bin_E; j++){
  		e_in[j] = 0;
  		count_E[j] = 0;
  		g_E[j] = 0;
  		for(int i=0; i<num_phys; i++){
  			physical_q[i][j] = 0;
  		}
  	}
  	fstream f2;
	stringstream name2;
	name2 << "eq_physical" << "/" << L 
	<< "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_" << step1 << "_" << seed << "_multi.txt";
	f2.open(name2.str().c_str(), ios::in);
	for(int j=50; j<bin_E; j++){
		int aa;
		f2 >> aa >> e_in[j] >> count_E[j] >> g_E[j]
			  >> physical_q[0][j] >> physical_q[1][j]
			  >> physical_q[2][j] >> physical_q[3][j]
			  >> physical_q[4][j] >> physical_q[5][j];
	}
	f2.close();
	for(int j=0; j<bin_E; j++){
		for(int j2=0; j2<bin_n; j2++){
			count_E_n[j][j2] = 0;
		}
	}
	fstream f3;
	stringstream name3;
	name3 << "eq_physical" << "/" << L 
	<< "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_" << step1 << "_" << seed << "_multi_nads.txt";
	f3.open(name3.str().c_str(), ios::in);
	for(int j=50; j<bin_E; j++){
		for(int j2=0; j2<bin_n; j2++){
			int aa, bb, cc;
			f3 >> aa >> bb >> cc >> count_E_n[j][j2];
		}
	}
	f3.close();
	double denominator;
	for(int T_int=0; T_int< 700; T_int++){
		T = 0.8 - T_int*0.001;
		beta = 1.0/T;
		double e_c = -5000;
		for(int j=50; j<bin_E; j++){
			if(e_c < -beta*e_in[j]*N + g_E[j]){
				e_c = -beta*e_in[j]*N + g_E[j];
			}
		}
		for(int j2=0; j2<bin_n; j2++){
			count_n_re[j2] = 0;
			energy_ave[j2] = 0;
		}
		denominator = 0.0;
		for(int j=50; j<bin_E; j++){
			for(int j2=0; j2<bin_n; j2++){
				denominator += count_E_n[j][j2]*exp(-beta*e_in[j]*N - e_c + g_E[j]);
				count_n_re[j2] += count_E_n[j][j2]*exp((-beta)*e_in[j]*N - e_c + g_E[j]);
				energy_ave[j2] += e_in[j]*count_E_n[j][j2]*exp((-beta)*e_in[j]*N - e_c + g_E[j]);
			}	
		}
		// cout << T << " " << denominator << endl;
		fstream f4;
	    stringstream name4;
	    name4 << "eq_physical/hist/" << L << "_" << mu_in
	    << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_" << (int)((T+0.00001)*10000) << "_" << step1 << "_" << seed << "_reweight_nads.txt";
	    f4.open(name4.str().c_str(), ios::out);
	    for(int j2=0; j2<bin_n; j2++){
	    	f4 << j2 << " " //1
	    	<< (j2 + 0.5)/(double)bin_n << " " //2
	    	<< count_n_re[j2]*bin_n/denominator << " " //3
	    	<< log(count_n_re[j2])+log(bin_n)-log(denominator) << " " //4
	    	<< log(count_n_re[j2]) << " " //5
	    	<< log(denominator) << " " //6
	    	<< energy_ave[j2]/count_n_re[j2] << endl; //7
	    }
	    f4.close();
	}
}