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

int count_E[bin_E];
double g_E[bin_E];
double delta_g_E;
int n_b[N][2*edge];
int sigma[N];
double position[N][2];
double a1;
int mu_in;
int Pre_in;
int alpha_in;
int k_in;
int seed;

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
	double alpha;
  if(alpha_in == 67){
    alpha = 2.0/3.0;
  }else{
    alpha = alpha_in/100.0;
  }
  double mu = mu_in/100.0;
  	double Pre = Pre_in/100.0;
  	double k = k_in/100.0;
	make_network(n_b);
	init_genrand(seed);
  	initial_conf(sigma, position, alpha,k,a1);
  	double e_max = 1.0;
  	double e_min = -mu + 3*k*alpha*alpha/(1+k);
  	double min_count;
  	double mean_count;
  	double bin_wid = (e_max-e_min)/(double)bin_E;
    int loop_c;
    for(int i=0; i<bin_E; i++){
      g_E[i] = 0;
    }
    fstream f2;
    stringstream name2;
    name2 << "eq_physical" << "/" << L << "_" << mu_in << "_" << Pre_in 
    << "_" << alpha_in << "_" << k_in << "_" << step1 << "_" << seed << "_w_l.txt";
    f2.open(name2.str().c_str(), ios::in);
    for(int i=50; i<bin_E; i++){
      int a,c;
      double b;
      f2 >> a >> b >> delta_g_E >> loop_c >> c >> g_E[i];
    }
    f2.close();
  	double max_g;
  	double E = enthalpy_next(sigma, position, n_b, alpha, k, mu, a1, Pre);
  	while(delta_g_E > 0.000001){
  		min_count = 100;
  		mean_count = 200;
  		for(int i=0; i<bin_E; i++){
	  		count_E[i] = 0;
	  	}
	  	loop_c = 0;
  		while(min_count < 0.8*mean_count){
  			metropolis_entropic_ads_next(sigma, position, n_b, alpha, k, mu, a1, E, e_max, e_min, bin_wid, count_E, g_E, delta_g_E);
  			for(int l=0; l<L; l++){
  				metropolis_entropic_position_next(sigma, position, n_b, alpha, k, a1, E, e_max, e_min, bin_wid, count_E, g_E, delta_g_E);
	  			metropolis_entropic_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, E, e_max, e_min, bin_wid, count_E, g_E, delta_g_E);
  			}
  			min_count = count_E[50];
  			mean_count = 0.0;
  			mean_count += count_E[50];
  			for(int i=51; i<bin_E; i++){
  				mean_count += count_E[i];
  				if(min_count > count_E[i]){
  					min_count = count_E[i];
  				}
  			}
  			mean_count /= (double)(bin_E-50);
  			if(loop_c%1000 == 0){
  				fstream f1;
				stringstream name1;
				name1 << "eq_physical" << "/" << L << "_" << mu_in << "_" << Pre_in 
				<< "_" << alpha_in << "_" << k_in << "_" << step1 << "_" << seed << "_w_l.txt";
				f1.open(name1.str().c_str(), ios::out);
				for(int i=50; i<bin_E; i++){
					f1 << i << " " << e_min + (i+0.5)*bin_wid << " " << delta_g_E << " " << loop_c << " " << count_E[i] << " " << g_E[i] << endl;
				}
				f1.close();
  			}
  			loop_c += 1;
  		}
  		max_g = 0;
  		for(int i=0; i<bin_E; i++){
  			if(max_g < g_E[i]){
  				max_g = g_E[i];
  			}
  		}
  		for(int i=0; i<bin_E; i++){
	  		g_E[i] -= max_g;
	  	}
  		delta_g_E /= 2.0;
  	}
  	fstream f1;
	stringstream name1;
	name1 << "eq_physical" << "/" << L << "_" << mu_in << "_" << Pre_in 
	<< "_" << alpha_in << "_" << k_in << "_" << step1 << "_" << seed << "_w_l.txt";
	f1.open(name1.str().c_str(), ios::out);
	for(int i=50; i<bin_E; i++){
		f1 << i << " " << e_min + (i+0.5)*bin_wid << " " << delta_g_E << " " << loop_c << " " << count_E[i] << " " << g_E[i] << endl;
	}
	f1.close();
}