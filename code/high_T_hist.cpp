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

const double T_start = 30.0;

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
    double alpha = alpha_in/100.0;
    double Pre = Pre_in/100.0;
    double k = k_in/100.0;
    make_network(n_b);
    init_genrand(seed);
    initial_conf(sigma, position, alpha,k,a1);
    double T = T_start;
    
    double E;
    fstream f1;
    stringstream name1;
    name1 << "eq_physical" << "/" << L << "_" << mu_in << "_" << Pre_in 
    << "_" << alpha_in << "_" << k_in << "_" << (int)((T+0.00001)*10000) 
    << "_" << step1 << "_" << seed << "_h_hist.txt";
    f1.open(name1.str().c_str(), ios::out);
    double beta;
    beta = 1.0/T;
    for(int k=0; k<step1; k++){
       metropolis_ads_next(sigma, position, n_b, alpha, k, mu, a1, beta);
        for(int l=0; l<L; l++){
            metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
            metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
        }
        E = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);
        f1 << k << " " << E/(double)N << endl;
    }
    f1.close();
}