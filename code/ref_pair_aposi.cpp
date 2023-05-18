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

int n_b[N][2*edge];
int sigma[N];
double position[N][2];
double a1;
int alpha_in = 60;
int seed;
double p_bond[N][4];
double p_plaquette[N];
double p_plaquette2[N];
double p_plaquette3[N];
double f_bond[N][4];

void initial_conf(int sigma[N], double position[N][2], double alpha, double k, double &a1, int p_dis){
  int sum_sig=0;
	for(int i=0; i<N; i++){
    if(i == L*L/2+L/4){
      sigma[i] = 1;
      sum_sig += 1;
    }else if(i == L*L/2+L/4+p_dis){
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
		seed = atoi(argv[1]);
	}else{
		return 1;
	}
	double mu = 0.0;
  double Pre = 0.0;
  double alpha = alpha_in/100.0;
	make_network(n_b);
	init_genrand(seed);
  double T = 0.0000000001;
  double E1, E2;

  fstream f1;
  stringstream name1;
  name1 << "eq_physical" << "/" << L << "_" << alpha_in << "_" << step1 << "_" << seed << "_pair_ref_aposi.txt";
  f1.open(name1.str().c_str(), ios::out);
  double beta = 1.0/T;
  for(int k_in = -50; k_in<10; k_in++){
    double k = k_in/100.0;
    int p_d = 1;
    initial_conf(sigma, position, alpha,k,a1,p_d);
    for(int t=0; t<step1; t++){
      metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
      metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
    }
    E1 = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);

    p_d = L/2-1;
    initial_conf(sigma, position, alpha,k,a1,p_d);
    for(int t=0; t<step1; t++){
      metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
      metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
    }
    E2 = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);
    f1 << k << " " << E1 << " " << E2 << " " << E1-E2 << endl;
  }
  for(int k_in = 1; k_in<71; k_in++){
    double k = k_in/10.0;
    int p_d = 1;
    initial_conf(sigma, position, alpha,k,a1,p_d);
    for(int t=0; t<step1; t++){
      metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
      metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
    }
    E1 = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);

    p_d = L/2-1;
    initial_conf(sigma, position, alpha,k,a1,p_d);
    for(int t=0; t<step1; t++){
      metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
      metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
    }
    E2 = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);
    f1 << k << " " << E1 << " " << E2 << " " << E1-E2 << endl;
  }
  f1.close();
}