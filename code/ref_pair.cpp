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
  << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_" << step1 << "_" << seed << "_pair_ref.txt";
  f1.open(name1.str().c_str(), ios::out);
  double beta;
  for(int p_d = 1; p_d < L/2; p_d++){
    beta = 1.0/T;
    initial_conf(sigma, position, alpha,k,a1,p_d);
    for(int t=0; t<step1; t++){
  		metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
      metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
    }
    E = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);
   	local_potential(sigma, position, n_b, alpha, k, mu, a1, p_plaquette, p_plaquette2, p_plaquette3,f_bond);
    double sum_p1 = 0;
    double sum_p2 = 0;
    double sum_p3 = 0;
    for(int i=0; i<N; i++){
      sum_p1 += p_plaquette2[i];
      if(sigma[i] == 0){
        sum_p2 += p_plaquette2[i];
      }else{
        sum_p3 += p_plaquette2[i];
      }
    }
    f1 << p_d << " " << E << " " << sum_p1 << " " << sum_p2 << " " << sum_p3 << endl;
    fstream f4;
    stringstream name4;
    name4 << "conf" << "/" << L << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_"
    << step1 << "_" << seed << "_" << p_d << "_ref_pair.txt";
    f4.open(name4.str().c_str(), ios::out);
    for(int i=0; i<N; i++){
      f4 << a1 << " " //1
        << position[i][0] << " " //2
        << position[i][1] << " "  //3
        << sigma[i] << " " //4
        << p_plaquette[i] << " " //5
        << p_plaquette2[i] << " " //6
        << p_plaquette3[i] << endl; //7
    }
    f4.close();
  }
  f1.close();
}