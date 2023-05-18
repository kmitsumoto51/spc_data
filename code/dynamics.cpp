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

int T_in = 1900;

int n_b[N][2*edge];
int sigma[N];
double position[N][2];
double position_eq[N][2];
double a1;
double a1_eq;
int mu_in;
int alpha_in;
int Pre_in;
int k_in;
int seed;
double p_plaquette[N];
double p_plaquette2[N];
double p_plaquette3[N];
double f_bond[N][4];
double p_vol[N];
double p_vol_eq[N];
double p_plaquette_eq[N];
double p_plaquette2_eq[N];
double p_plaquette3_eq[N];
double f_bond_eq[N][4];
double debye_waller[N];
int contact_f[N][edge];
int contact_v[N][edge];
int max_cluster_size_f;
int max_cluster_size_v;
int cluster_num_f[N];
int cluster_num_v[N];
int num_cluster_f[N];
int num_cluster_v[N];
double c_ratio_f[2];
double c_ratio_v[2];

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
  double k = k_in/100.0;
  double Pre = Pre_in/100.0;
	make_network(n_b);
	init_genrand(seed);
  double T = T_in/10000.0;;
  double beta;
  double volu;
  double E;
  int Nads;
  fstream f1;
  stringstream name1;
  name1 << "conf" << "/" << L << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_"
  << T_in+100 << "_" << step1 << "_" << seed << ".txt";
  f1.open(name1.str().c_str(), ios::in);
  double dami1, dami2;
  f1 >> a1 >> dami1 >> dami2;
  for(int i=0; i<N; i++){
    f1 >> position[i][0] >> position[i][1] >> sigma[i];
  }
  f1.close();
  beta = 1.0/T;
  fstream f2;
  stringstream name2;
  name2 << "eq_physical" << "/" << L 
  << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in 
  << "_" << T_in << "_" << step1 << "_" << seed << "_dynamics.txt";
  f2.open(name2.str().c_str(), ios::out);
  for(int t=0; t<500; t++){
		metropolis_ads_next(sigma, position, n_b, alpha, k, mu, a1, beta);
    for(int l=0; l<L; l++){
      metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
      metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
    }
    Nads = num_ads(sigma);
    E = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);
    volu = a1*a1*L*L;
    contact_list(sigma, n_b, contact_f, contact_v);
    cluster_ratio_f(sigma, n_b, cluster_num_f, c_ratio_f, contact_f);
    cluster_ratio_v(sigma, n_b, cluster_num_v, c_ratio_v, contact_v);
    for(int i=0; i<N; i++){
      num_cluster_f[i] = 0;
      num_cluster_v[i] = 0;
    }
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        if(cluster_num_f[i] == j){
          num_cluster_f[j] += 1;
        }
        if(cluster_num_v[i] == j){
          num_cluster_v[j] += 1;
        }
      }
    }
    int max_size_f = 0;
    int max_size_v = 0;
    int max_c_f;
    int max_c_v;
    for(int i=0; i<N; i++){
      if(max_size_f < num_cluster_f[i]){
        max_size_f = num_cluster_f[i];
        max_c_f = i;
      }
      if(max_size_v < num_cluster_v[i]){
        max_size_v = num_cluster_v[i];
        max_c_v = i;
      }
    }
    double c_ratio_max_f = 4*max_size_f;
    double c_ratio_max_v = 4*max_size_v;
    for(int i=0; i<N; i++){
      if(max_c_f == cluster_num_f[i]){
        for(int j=0; j<4; j++){
          c_ratio_max_f -= contact_f[i][j];
        }
      }
      if(max_c_v == cluster_num_v[i]){
        for(int j=0; j<4; j++){
          c_ratio_max_v -= contact_v[i][j];
        }
      }
    }
    c_ratio_max_f /= max_size_f;
    c_ratio_max_v /= max_size_v;
    f2 << t << " " << Nads << " " << E << " " << volu << " " << c_ratio_f[0] << " " << c_ratio_f[1]
    << " " << c_ratio_v[0] << " " << c_ratio_v[1] << " " << max_size_f << " " << c_ratio_max_f
    << " " << max_size_v << " " << c_ratio_max_v << endl;
    fstream f4;
    stringstream name4;
    name4 << "conf" << "/" << L << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_"
    << T_in << "_" << step1 << "_" << t << "_" << seed << "_dynamics.txt";
    f4.open(name4.str().c_str(), ios::out);
    f4 << a1 << " " << 0 << " " << 0 << endl;
    for(int i=0; i<N; i++){
      f4 << position[i][0] << " " << position[i][1] << " " << sigma[i] << endl;
    }
    f4.close();
  }
  f2.close();
}