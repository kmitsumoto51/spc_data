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

const double T_start = 0.02;
const double T_finish = 0.0;
const int ave_t = 1000;

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
int N_per;
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
double c_ratio_f[2];
double c_ratio_v[2];


int main(int argc, char **argv){
	if(argc > 1){
		N_per = atoi(argv[1]);
    Pre_in = atoi(argv[2]);
    alpha_in = atoi(argv[3]);
    k_in = atoi(argv[4]);
    seed = atoi(argv[5]);
	}else{
		return 1;
	}
  double mu = 0.0;
	double alpha = -alpha_in/100.0;
  double k = -k_in/100.0;
  double Pre = Pre_in/100.0;
	make_network(n_b);
	init_genrand(seed);
  double T = T_start;
  double beta;
  while(T > T_finish+0.00001){
    a1_eq = 0;
    for(int i=0; i<N; i++){
      position_eq[i][0] = 0;
      position_eq[i][1] = 0;
      for(int j=0; j<4; j++){
        f_bond_eq[i][j] = 0;
      }
      p_plaquette_eq[i] = 0;
      p_plaquette2_eq[i] = 0;
      p_plaquette3_eq[i] = 0;
      p_vol_eq[i] = 0;
      debye_waller[i] = 0;
    }
    fstream f1;
    stringstream name1;
    name1 << "conf" << "/" << L << "_" << N_per << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_"
    << (int)((T+0.00001)*10000) << "_" << step1 << "_" << seed << "_both_negative_Nfix.txt";
    f1.open(name1.str().c_str(), ios::in);
    double dami1, dami2;
    f1 >> a1 >> dami1 >> dami2;
    for(int i=0; i<N; i++){
      f1 >> position[i][0] >> position[i][1] >> sigma[i];
    }
    f1.close();
    beta = 1.0/T;
    for(int t=0; t<ave_t; t++){
  		metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
      metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
      local_potential(sigma, position, n_b, alpha, k, mu, a1, p_plaquette, p_plaquette2, p_plaquette3,f_bond);
      plaquette_volume(position, n_b, a1, p_vol);
      a1_eq += a1;
      for(int i=0; i<N; i++){
        position_eq[i][0] += position[i][0];
        position_eq[i][1] += position[i][1];
        for(int j=0; j<4; j++){
          f_bond_eq[i][j] += f_bond[i][j];
        }
        p_plaquette_eq[i] += p_plaquette[i];
        p_plaquette2_eq[i] += p_plaquette2[i];
        p_plaquette3_eq[i] += p_plaquette3[i];
        p_vol_eq[i] += p_vol[i];
      }
    }
    for(int t=0; t<ave_t; t++){
      metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
      metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
      for(int i=0; i<N; i++){
        debye_waller[i] += (position_eq[i][0]/ave_t - position[i][0])*(position_eq[i][0]/ave_t - position[i][0])
                         + (position_eq[i][1]/ave_t - position[i][1])*(position_eq[i][1]/ave_t - position[i][1]);
      }
    }
    contact_list(sigma, n_b, contact_f, contact_v);
    cluster_ratio_f(sigma, n_b, cluster_num_f, c_ratio_f, contact_f);
    cluster_ratio_v(sigma, n_b, cluster_num_v, c_ratio_v, contact_v);
    fstream f4;
    stringstream name4;
    name4 << "conf" << "/" << L << "_" << N_per << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_"
    << (int)((T+0.00001)*10000) << "_" << step1 << "_" << seed << "_d_noise_both_negative_Nfix.txt";
    f4.open(name4.str().c_str(), ios::out);
    for(int i=0; i<N; i++){
      f4 << a1_eq/ave_t << " " //1
      << position_eq[i][0]/ave_t << " " //2
      << position_eq[i][1]/ave_t << " "  //3
      << sigma[i] << " " //4
      << f_bond_eq[i][0]/ave_t << " "  //5
      << f_bond_eq[i][1]/ave_t << " "  //6
      << f_bond_eq[i][2]/ave_t << " "  //7
      << f_bond_eq[i][3]/ave_t << " " //8
      << p_plaquette_eq[i]/ave_t << " " //9
      << p_plaquette2_eq[i]/ave_t << " " //10
      << p_plaquette3_eq[i]/ave_t << " " //11
      << p_vol_eq[i]/ave_t << " " //12
      << debye_waller[i]/ave_t << " "//13
      << cluster_num_f[i] << " "//14
      << cluster_num_v[i] << " "//15
      << c_ratio_f[0] << " "//16
      << c_ratio_f[1] << " "//17
      << c_ratio_v[0] << " "//18 
      << c_ratio_v[1] << endl; //19
    }
    f4.close();
    T -= 0.01;
  }
}