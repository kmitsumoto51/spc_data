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

const int ave_t = 1000;

int n_b[N][2*edge];
int sigma[N];
double position[N][2];
double a1;
int mu_in;
int alpha_in;
int Pre_in;
int k_in;
int seed;
double posix1, posiy1, posix2, posiy2, posix3, posiy3, posix4, posiy4;
double p_vol[N];


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
  double T = 0.2;
  double beta;
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
  double sum_length_a = 0.0;
  double sum_length_a_sq = 0.0;
  double sum_length_d = 0.0;
  double sum_length_d_sq = 0.0;
  double length_a_ave, length_d_ave;
  double length_a_sq_ave, length_d_sq_ave;
  double length_sq1, length_sq2, length_sq3, length_sq4;
  double sum_vol_a = 0.0;
  double sum_vol_a_sq = 0.0;
  double sum_vol_d = 0.0;
  double sum_vol_d_sq = 0.0;
  double vol_a_ave, vol_d_ave;
  double vol_a_sq_ave, vol_d_sq_ave;
  int sum_sig = 0;
  for(int t=0; t<ave_t; t++){
		metropolis_ads_next(sigma, position, n_b, alpha, k, mu, a1, beta);
    for(int l=0; l<L; l++){
      metropolis_position_next(sigma, position, n_b, alpha, k, a1, beta);
      metropolis_volume_next(sigma, position, n_b, alpha, k, mu, a1, Pre, beta);
    }
    plaquette_volume(position, n_b, a1, p_vol);
    for(int i=0; i<N; i++){
      posix1 = position[i][0];
      posiy1 = position[i][1];
      posix2 = position[n_b[i][0]][0];
      posiy2 = position[n_b[i][0]][1];
      posix3 = position[n_b[i][4]][0];
      posiy3 = position[n_b[i][4]][1];
      posix4 = position[n_b[i][1]][0];
      posiy4 = position[n_b[i][1]][1];
      if(i%L==L-1){
        posix2 = position[n_b[i][0]][0] + L*a1;
        posix3 = position[n_b[i][4]][0] + L*a1;
      }
      if(i/L==L-1){
        posiy3 = position[n_b[i][4]][1] + L*a1;
        posiy4 = position[n_b[i][1]][1] + L*a1;
      }
      length_sq1 = (posix1-posix2)*(posix1-posix2)+(posiy1-posiy2)*(posiy1-posiy2);
      length_sq2 = (posix2-posix3)*(posix2-posix3)+(posiy2-posiy3)*(posiy2-posiy3);
      length_sq3 = (posix3-posix4)*(posix3-posix4)+(posiy3-posiy4)*(posiy3-posiy4);
      length_sq4 = (posix4-posix1)*(posix4-posix1)+(posiy4-posiy1)*(posiy4-posiy1);
      if(sigma[i]==1){
        sum_length_a += sqrt(length_sq1)+sqrt(length_sq2)+sqrt(length_sq3)+sqrt(length_sq4);
        sum_length_a_sq += length_sq1+length_sq2+length_sq3+length_sq4;
        sum_vol_a += p_vol[i];
        sum_vol_a_sq += p_vol[i]*p_vol[i];
        sum_sig += 1;
      }else{
        sum_length_d += sqrt(length_sq1)+sqrt(length_sq2)+sqrt(length_sq3)+sqrt(length_sq4);
        sum_length_d_sq += length_sq1+length_sq2+length_sq3+length_sq4;
        sum_vol_d += p_vol[i];
        sum_vol_d_sq += p_vol[i]*p_vol[i];
      }
    }
  }
  length_a_ave = sum_length_a/(4.0*sum_sig);
  length_d_ave = sum_length_d/(4.0*(N*ave_t-sum_sig));
  length_a_sq_ave = sum_length_a_sq/(4.0*sum_sig);
  length_d_sq_ave = sum_length_d_sq/(4.0*(N*ave_t-sum_sig));

  vol_a_ave = sum_vol_a/(sum_sig);
  vol_d_ave = sum_vol_d/((N*ave_t-sum_sig));
  vol_a_sq_ave = sum_vol_a_sq/(sum_sig);
  vol_d_sq_ave = sum_vol_d_sq/((N*ave_t-sum_sig));

  cout << length_a_ave << " " << length_d_ave << " " 
  << (length_a_sq_ave - length_a_ave*length_a_ave) << " " 
  << (length_d_sq_ave - length_d_ave*length_d_ave) << " "
  << vol_a_ave << " " << vol_d_ave << " "
  << (vol_a_sq_ave - vol_a_ave*vol_a_ave) << " " 
  << (vol_d_sq_ave - vol_d_ave*vol_d_ave) << endl;
}