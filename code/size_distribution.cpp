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

int T_in = 2000;

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
int contact_f[N][edge];
int contact_v[N][edge];
int cluster_num_f[N];
double c_ratio_f[2];
int num_cluster_f[N];
double bin_size[N+1];
double bin_max_size[N];
double prob_s_cluster[N];


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
  double size_sum, max_size_sum;
  int size_count;
  beta = 1.0/T;
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
  size_sum = 0;
  size_count = 0;
  max_size_sum = 0;
  for(int i=0; i<N; i++){
    bin_size[i] = 0;
    bin_max_size[i] = 0;
  }
  for(int t=0; t<20000; t++){
		metropolis_ads_next(sigma, position, n_b, alpha, k, mu, a1, beta);
    for(int l=0; l<L; l++){
      metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
      metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
    }
    contact_list(sigma, n_b, contact_f, contact_v);
    cluster_ratio_f(sigma, n_b, cluster_num_f, c_ratio_f, contact_f);
    for(int i=0; i<N; i++){
      num_cluster_f[i] = 0;
    }
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        if(cluster_num_f[i] == j){
          num_cluster_f[j] += 1;
        }
      }
    }
    int max_size_f = 0;
    for(int i=0; i<N; i++){
      if(max_size_f < num_cluster_f[i]){
        max_size_f = num_cluster_f[i];
      }
    }
    max_size_sum += max_size_f;
    for(int i=0; i<N; i++){
      if(i==max_size_f){
        bin_max_size[i] += 1.0;
      }
    }
    for(int i=0; i<N; i++){
      if(num_cluster_f[i] != 0){
        bin_size[num_cluster_f[i]] += 1;
        prob_s_cluster[num_cluster_f[i]] += num_cluster_f[i];
        size_sum += num_cluster_f[i];
        size_count += 1;
      }
    }
  }
  for(int i=0; i<N; i++){
    prob_s_cluster[i] /= 20000.0*N;
  }
  double c_size_sum = 0.0;
  double denom = 0.0;
  for(int i=0; i<N; i++){
    c_size_sum += i*prob_s_cluster[i];
    denom += prob_s_cluster[i];
  }
  c_size_sum /= denom;
  fstream f2;
  stringstream name2;
  name2 << "eq_physical" << "/" << L 
  << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in 
  << "_" << T_in << "_" << step1 << "_" << seed << "_dist_f.txt";
  f2.open(name2.str().c_str(), ios::out);
  for(int i=0; i<N; i++){
    f2 << i << " " << size_sum/(double)size_count << " " << bin_size[i]/(double)size_count 
    << " " << max_size_sum/20000.0 << " " << bin_max_size[i]/20000.0
    << " " << c_size_sum << " " << prob_s_cluster[i] << endl;
  }
  f2.close();
}