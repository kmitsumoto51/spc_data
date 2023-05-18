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

const double T_start = 0.6;
const double T_finish = 0.3;
const int step_c = 10000;

int n_b[N][2*edge];
int sigma[N];
double position[N][2];
double a1;
double p_vol[N];
int mu_in;
int Pre_in;
int alpha_in;
int k_in;
int seed;
int vol_hist[num_bin];
int size_hist_v[N-1]; 
double ratio_hist_v[N-1];
int size_hist_f[N-1]; 
double ratio_hist_f[N-1];
int contact_f[N][edge];
int contact_v[N][edge];
int cluster_num_f[N];
int cluster_num_v[N];

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
  double T = T_start;
  double beta;
  int total_f;
  int total_v;
  while(T > T_finish+0.00001){
    for(int j=0; j<num_bin; j++){
      vol_hist[j] = 0.0;
    }
    for(int i=0; i<N-1; i++){
      size_hist_f[i] = 0.0;
      size_hist_v[i] = 0.0;
      ratio_hist_f[i] = 0.0;
      ratio_hist_v[i] = 0.0;
    }
    total_f = 0;
    total_v = 0;
    fstream f1;
    stringstream name1;
    name1 << "conf" << "/" << L << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_"
    << (int)((T+0.00001)*10000) << "_" << step1 << "_" << seed << ".txt";
    f1.open(name1.str().c_str(), ios::in);
    double dami1, dami2;
    f1 >> a1 >> dami1 >> dami2;
    for(int i=0; i<N; i++){
      f1 >> position[i][0] >> position[i][1] >> sigma[i];
    }
    f1.close();
    beta = 1.0/T;
    for(int t=0; t<step_c; t++){
      metropolis_ads_next(sigma, position, n_b, alpha, k, mu, a1, beta);
      for(int l=0; l<L; l++){
        metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
        metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
      }
      plaquette_volume(position, n_b, a1, p_vol);
      contact_list(sigma, n_b, contact_f, contact_v);
      cluster_ratio_hist_f(sigma, n_b, cluster_num_f, contact_f, size_hist_f, ratio_hist_f);
      cluster_ratio_hist_v(sigma, n_b, cluster_num_v, contact_v, size_hist_v, ratio_hist_v);
      for(int j=0; j<num_bin; j++){
        for(int i=0; i<N; i++){
          if(p_vol[i]<=(j+1)*3.0/num_bin && p_vol[i]>j*3.0/num_bin){
            vol_hist[j] += 1;
          }
        }
      }
    }
    for(int i=0; i<N-1; i++){
      total_f += size_hist_f[i];
      total_v += size_hist_v[i];
      if(size_hist_f[i] == 0){
        ratio_hist_f[i] = 0;
      }else{
        ratio_hist_f[i] /= size_hist_f[i];
      }
      if(size_hist_v[i] == 0){
        ratio_hist_v[i] = 0;
      }else{
        ratio_hist_v[i] /= size_hist_v[i];
      }
    }

    fstream f2;
    stringstream name2;
    name2 << "hist/" << L << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_" 
    << (int)((T+0.00001)*10000) << "_" << seed << ".txt";
    f2.open(name2.str().c_str(), ios::out);
    for(int j=0; j<num_bin; j++){
      f2 << (j+0.5)/(num_bin/3.0) << " " << vol_hist[j]/((double)N*step_c/num_bin) << endl;
    }
    f2.close();
    fstream f3;
    stringstream name3;
    name3 << "hist/" << L << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in << "_" 
    << (int)((T+0.00001)*10000) << "_" << seed << "_cluster.txt";
    f3.open(name3.str().c_str(), ios::out);
    for(int j=0; j<N-1; j++){
      f3 << j+1 << " " << size_hist_f[j]/(double)total_f << " " << size_hist_v[j]/(double)total_v << " " << ratio_hist_f[j] << " " << ratio_hist_v[j] << endl;
    }
    f3.close();
    T -= 0.01;
  }
}