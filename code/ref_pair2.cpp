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
int mu_in;
int Pre_in;
int alpha_in;
int k_in;
int seed;
double p_bond[N][4];
double p_plaquette[N];
double f_bond[N][4];

void initial_conf(int sigma[N], double position[N][2], double alpha, double k, double &a1, int p_dis){
  int sum_sig=0;
	for(int i=0; i<N; i++){
    if(i == L*L/2+L/4){
      sigma[i] = 1;
      sum_sig += 1;
    }else if(i == L*L/2+L/4+p_dis*(L+1)){
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
	double alpha = alpha_in/100.0;
  double Pre = Pre_in/100.0;
  double k = k_in/100.0;
	make_network(n_b);
	init_genrand(seed);
  double T = 0.0000000001;
  double E;

  fstream f1;
  stringstream name1;
  name1 << "eq_physical" << "/" << L 
  << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_" << step1 << "_" << seed << "_pair_ref2.txt";
  f1.open(name1.str().c_str(), ios::out);
  double beta= 1.0/T;
  for(int p_d = 1; p_d < L/2; p_d++){
    initial_conf(sigma, position, alpha,k,a1,p_d);
    for(int t=0; t<step1; t++){
  		metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
      metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
    }
    E = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);
   	f1 << p_d << " " << E/(double)N << endl;
  }
  f1.close();
}
