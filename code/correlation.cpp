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

const double T_start = 1.0;
const double T_finish = 0.0;
const int step_c = 1000;
const int num_seed = 5;

int n_b[N][2*edge];
int sigma[N];
double position[N][2];
int mu_in;
int Pre_in;
int alpha_in;
int k_in;
int seed1;
double a1;
double p_bond[N][4];
double p_plaquette[N];
double p_plaquette2[N];
double p_plaquette3[N];
double f_bond[N][4];
double f_correlation[L/2][4];
double ads_correlation[L/2];
double f_correlation_eq[L/2][4];
double ads_correlation_eq[L/2];
double f_ave[2];
double f_ave_eq[2];

int main(int argc, char **argv){
	if(argc > 1){
		mu_in = atoi(argv[1]);
    Pre_in = atoi(argv[2]);
    alpha_in = atoi(argv[3]);
    k_in = atoi(argv[4]);
    seed1 = atoi(argv[5]);
	}else{
		return 1;
	}
	double mu = mu_in/100.0;
  double alpha = alpha_in/100.0;
  double Pre = Pre_in/100.0;
  double k = k_in/100.0;
	make_network(n_b);
  double T = T_start;
  int Nads;
  double Nads_eq;
  double beta;
  while(T > T_finish+0.00001){
    for(int seed2=1; seed2<num_seed+1; seed2++){
      Nads_eq = 0.0;
      f_ave_eq[0] = 0.0;
      f_ave_eq[1] = 0.0;
      for(int j=0; j<L/2; j++){
        f_correlation_eq[j][0] = 0.0;
        f_correlation_eq[j][1] = 0.0;
        f_correlation_eq[j][2] = 0.0;
        f_correlation_eq[j][3] = 0.0;
        ads_correlation_eq[j] = 0.0;
      }
      init_genrand(seed2);
      fstream f1;
      stringstream name1;
      name1 << "conf" << "/" << L << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_"
      << (int)((T+0.00001)*10000) << "_" << step1 << "_" << seed1 << ".txt";
      f1.open(name1.str().c_str(), ios::in);
      double dami1, dami2;
      f1 >> a1 >> dami1 >> dami2;
      for(int i=0; i<N; i++){
        f1 >> position[i][0] >> position[i][1] >> sigma[i];
      }
      f1.close();
      Nads = num_ads(sigma);
      beta = 1.0/T;
      for(int t=0; t<step_c; t++){
        metropolis_ads_next(sigma, position, n_b, alpha, k, mu, a1, beta);
        for(int l=0; l<L; l++){
          metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
          metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
        }
        Nads = num_ads(sigma);
        local_potential(sigma, position, n_b, alpha, k, mu, a1, p_plaquette, p_plaquette2, p_plaquette3,f_bond);
        force_correlation(f_bond, f_correlation);
        particle_correlation(sigma, ads_correlation);
        Nads_eq += Nads/(double)N;
        f_ave_eq[0] += f_ave[0];
        f_ave_eq[1] += f_ave[1];
        for(int j=0; j<L/2; j++){
          f_correlation_eq[j][0] += f_correlation[j][0];
          f_correlation_eq[j][1] += f_correlation[j][1];
          f_correlation_eq[j][2] += f_correlation[j][2];
          f_correlation_eq[j][3] += f_correlation[j][3];
          ads_correlation_eq[j] += ads_correlation[j];
        }
      }
      fstream f2;
      stringstream name2;
      name2 << "correlation/" << L << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_" 
      << (int)((T+0.00001)*10000) << "_" << num_seed*(seed1-1)+seed2 << ".txt";
      f2.open(name2.str().c_str(), ios::out);
      for(int j=0; j<L/2; j++){
        f2 << f_ave_eq[0]/(double)step_c << " " //1
           << f_ave_eq[1]/(double)step_c << " " //2
           << Nads_eq/(double)step_c << " " //3
           << f_correlation_eq[j][0]/(double)step_c << " " //4
           << f_correlation_eq[j][1]/(double)step_c << " " //5
           << f_correlation_eq[j][2]/(double)step_c << " " //6
           << f_correlation_eq[j][3]/(double)step_c << " " //7
           << ads_correlation_eq[j]/(double)step_c << endl; //8
      }
    }
    T -= 0.01;
  }
}