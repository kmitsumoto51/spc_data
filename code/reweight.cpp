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

double e_in[bin_E];
int count_E[bin_E];
double count_E_re[bin_E];
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
	fstream f1;
	stringstream name1;
	name1 << "eq_physical" << "/" << L 
	<< "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_" << step1 << "_" << seed << "_eq.txt";
	f1.open(name1.str().c_str(), ios::out);
	double denominator;
	for(int T_int=0; T_int< 700; T_int++){
		T = 0.8 - T_int*0.001;
		beta = 1.0/T;
		double e_c = -5000;
		// if(T>0.3){
		// 	e_c = -0.5*(e_max+e_min)*beta*N;
		// }else{
		// 	e_c = -e_min*beta*N;
		// }
		for(int j=50; j<bin_E; j++){
			if(e_c < -beta*e_in[j]*N + g_E[j]){
				e_c = -beta*e_in[j]*N + g_E[j];
			}
		}
		for(int j=0; j<num_phys; j++){
			physical_q_ave[j] = 0.0;
		}
		denominator = 0.0;
		for(int j=50; j<bin_E; j++){
			denominator += count_E[j]*exp(-beta*e_in[j]*N - e_c + g_E[j]);
			count_E_re[j] = count_E[j]*exp((-beta)*e_in[j]*N - e_c + g_E[j]);
			for(int i=0; i<num_phys; i++){
				physical_q_ave[i] += physical_q[i][j]*count_E[j]*exp((-beta)*e_in[j]*N - e_c + g_E[j]);
			}
		}
		// cout << T << " " << denominator << endl;
		for(int i=0; i<num_phys; i++){
			physical_q_ave[i] /= denominator;
		}
		for(int j=0; j<bin_E; j++){
			count_E_re[j] /= denominator;
		}
		fstream f4;
	    stringstream name4;
	    name4 << "eq_physical/hist/" << L << "_" << mu_in
	    << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_" << (int)((T+0.00001)*10000) << "_" << step1 << "_" << seed << "_reweight.txt";
	    f4.open(name4.str().c_str(), ios::out);
	    for(int j=50; j<bin_E; j++){
	    	f4 << j << " " << e_in[j] << " " << count_E_re[j]/bin_wid << endl;
	    }
	    f4.close();
	    E_ave = physical_q_ave[0];
	    E_sq_ave = physical_q_ave[1];
	    C=specific_heat(E_sq_ave,E_ave,beta);
	    Nads_ave = physical_q_ave[2];
	    Nads_sq_ave = physical_q_ave[3];
	    Nads_ave_sq = Nads_ave*Nads_ave;
	    chi = susceptibility(Nads_sq_ave,Nads_ave_sq,beta);
	    volume_ave = physical_q_ave[4];
	    volume_sq_ave = physical_q_ave[5];
	    volume_ave_sq = volume_ave*volume_ave;
	    kappa = susceptibility(volume_sq_ave,volume_ave_sq,beta)/volume_ave;
	    f1 << T//1
         << " " << E_ave/(double)N//2
         << " " << C/(double)N//3
         << " " << Nads_ave/(double)N//4
         << " " << Nads_sq_ave/((double)N*N)//5
         << " " << chi/(double)N//6
         << " " << volume_ave//7
         << " " << volume_sq_ave//8
         << " " << kappa//9
         << endl;
	}
	f1.close();
}
