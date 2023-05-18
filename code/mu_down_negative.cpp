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

const double mu_start = 1.5;
const double mu_finish = 0.4;

int n_b[N][2*edge];
int sigma[N];
double position[N][2];
double a1;
double physical_q[num_phys];
int T_in;
int Pre_in;
int alpha_in;
int k_in;
int seed;
double p_bond[N][4];
double p_plaquette[N];
double f_bond[N][4];

void initial_conf(int sigma[N], double position[N][2], double alpha, double k, double &a1){
  int sum_sig=0;
	for(int i=0; i<N; i++){
    sigma[i] = 1;
    sum_sig += 1;
  }
  a1 = 1 + alpha*k*sum_sig/((1+k)*N);
  for(int i=0; i<N; i++){
    position[i][0] = (i%L)*a1;
    position[i][1] = (i/L)*a1;
  }
}

int main(int argc, char **argv){
	if(argc > 1){
		T_in = atoi(argv[1]);
    Pre_in = atoi(argv[2]);
		alpha_in = atoi(argv[3]);
    k_in = atoi(argv[4]);
		seed = atoi(argv[5]);
	}else{
		return 1;
	}
	double T = T_in/10000.0;
	double alpha = -alpha_in/100.0;
  double Pre = Pre_in/100.0;
  double k = k_in/100.0;
	make_network(n_b);
	init_genrand(seed);
  initial_conf(sigma, position, alpha,k,a1);
  double mu = mu_start;
  double E, E_ave, E_sq_ave, C;
  int Nads;
  double Nads_ave, Nads_sq_ave, Nads_ave_sq, chi;
  double volu;
  double volume_ave, volume_sq_ave, volume_ave_sq, kappa;

  fstream f1;
  stringstream name1;
  name1 << "eq_physical" << "/" << L 
  << "_" << T_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_" << step1 << "_" << seed << "_mu_down_negative.txt";
  f1.open(name1.str().c_str(), ios::out);
  double beta = 1.0/T;
  while(mu > mu_finish+0.00001){
    for(int i=0; i<num_phys; i++){
      physical_q[i] = 0.0;
    }
    for(int t=0; t<step1; t++){
    	metropolis_ads_next(sigma, position, n_b, alpha, k, mu, a1, beta);
    	for(int l=0; l<L; l++){
    		metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
        metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
    	}
    }
    for(int t=0; t<step2; t++){
      metropolis_ads_next(sigma, position, n_b, alpha, k, mu, a1, beta);
      for(int l=0; l<L; l++){
        metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
        metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
      }
      Nads = num_ads(sigma);
      E = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);
      volu = a1*a1*L*L;
    	physical_q[0] += E;
    	physical_q[1] += E*E;
    	physical_q[2] += Nads;
    	physical_q[3] += Nads*Nads;
      physical_q[4] += volu;
      physical_q[5] += volu*volu;
    }
    E_ave = physical_q[0]/(double)(step2);
    E_sq_ave = physical_q[1]/(double)(step2);
    C=specific_heat(E_sq_ave,E_ave,beta);
    Nads_ave = physical_q[2]/(double)(step2);
    Nads_sq_ave = physical_q[3]/(double)(step2);
    Nads_ave_sq = Nads_ave*Nads_ave;
    chi = susceptibility(Nads_sq_ave,Nads_ave_sq,beta);
    volume_ave = physical_q[4]/(double)(step2);
    volume_sq_ave = physical_q[5]/(double)(step2);
    volume_ave_sq = volume_ave*volume_ave;
    kappa = susceptibility(volume_sq_ave,volume_ave_sq,beta)/volume_ave;
   	f1 << mu//1
         << " " << E_ave/(double)N//2
         << " " << C/(double)N//3
         << " " << Nads_ave/(double)N//4
         << " " << Nads_sq_ave/((double)N*N)//5
         << " " << chi/(double)N//6
         << " " << volume_ave//7
         << " " << volume_sq_ave//8
         << " " << kappa//9
         << endl;
    fstream f4;
    stringstream name4;
    name4 << "conf" << "/" << L << "_" << (int)((mu+0.0001)*100) << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_"
    << T_in << "_" << step1 << "_" << seed << "_mu_down_negative.txt";
    f4.open(name4.str().c_str(), ios::out);
    f4 << a1 << " " << 0 << " " << 0 << endl;
    for(int i=0; i<N; i++){
      f4 << position[i][0] << " " << position[i][1] << " " << sigma[i] << endl;
    }
    f4.close();
    mu -= 0.01;
  }
}