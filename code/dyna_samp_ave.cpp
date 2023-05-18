#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

const int T_in = 1900;
const int L = 48;
const int mu_in = 170;
const int alpha_in = 60;
const int Pre_in = 0;
const int k_in = 500;
const int step1 = 20000;
const int sample = 150;
const int num_bin = 100;
double phys[12];
double cmax_f;
double cmax_v;
int count_max_f[num_bin];
int count_max_v[num_bin];
double max_f_sum[num_bin];
double max_v_sum[num_bin];
double max_f_sq_sum[num_bin];
double max_v_sq_sum[num_bin];
double st_div_f[num_bin];
double st_div_v[num_bin];
double nads_ave[2000];
double nads_sq_ave[2000];

int main(){
	for(int i=0; i<num_bin; i++){
		max_f_sum[i] = 0;
		max_v_sum[i] = 0;
		max_f_sq_sum[i] = 0;
		max_v_sq_sum[i] = 0;
	}
	for(int i=0; i<2000; i++){
		nads_ave[i] = 0.0;
	}
	for(int samp = 1; samp<sample+1; samp++){
		fstream f2;
	    stringstream name2;
	    name2 << "eq_physical" << "/" << L 
	    << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in 
	    << "_" << T_in << "_" << step1 << "_" << samp << "_dynamics_sample.txt";
	    f2.open(name2.str().c_str(), ios::in);
	    for(int t=0; t<2000; t++){
	    	f2 >> phys[0] >> phys[1] >> phys[2] >> phys[3] >> phys[4] >> phys[5] 
	    	   >> phys[6] >> phys[7] >> phys[8] >> phys[9] >> phys[10] >> phys[11];
	    	cmax_f = phys[8]/((double)L*L);
	    	cmax_v = phys[10]/((double)L*L);
	    	nads_ave[t] += phys[1]/(double)(L*L);
	    	nads_sq_ave[t] += (phys[1]/(double)(L*L))*(phys[1]/(double)(L*L));
	    	if(cmax_f <= 1/(double)num_bin){
	    		count_max_f[0] += 1;
	    		max_f_sum[0] += phys[9];
	    		max_f_sq_sum[0] += phys[9]*phys[9];
	    	}
	    	if(cmax_v <= 1/(double)num_bin){
	    		count_max_v[0] += 1;
	    		max_v_sum[0] += phys[11];
	    		max_v_sq_sum[0] += phys[11]*phys[11];
	    	}
	    	for(int i=1; i<num_bin; i++){
	    		if(cmax_f <= (i+1)/(double)num_bin && cmax_f > i/(double)num_bin){
		    		count_max_f[i] += 1;
		    		max_f_sum[i] += phys[9];
		    		max_f_sq_sum[i] += phys[9]*phys[9];
		    	}
		    	if(cmax_v <= (i+1)/(double)num_bin && cmax_v > i/(double)num_bin){
		    		count_max_v[i] += 1;
		    		max_v_sum[i] += phys[11];
		    		max_v_sq_sum[i] += phys[11]*phys[11];
		    	}
	    	}
	    }
	    f2.close();
	}
	for(int i=0; i<2000; i++){
    	nads_ave[i] /= (double)sample;
    	nads_sq_ave[i] /= (double)sample;
    }
	fstream f3;
    stringstream name3;
    name3 << "eq_physical" << "/" << L 
    << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in 
    << "_" << T_in << "_" << step1 << "_dynamics_nads_ave.txt";
    f3.open(name3.str().c_str(), ios::out);
    for(int i=0; i<2000; i++){
    	f3 << i << " " << nads_ave[i] << " " << sqrt((nads_sq_ave[i] - nads_ave[i]*nads_ave[i])/(double)sample) << endl;
    }
    f3.close();
	for(int i=0; i<num_bin; i++){
		if(count_max_f[i] != 0){
			max_f_sum[i] /= count_max_f[i];
			max_f_sq_sum[i] /= count_max_f[i];
			st_div_f[i] = sqrt((max_f_sq_sum[i] - max_f_sum[i]*max_f_sum[i])/count_max_f[i]);
		}
		if(count_max_v[i] != 0){
			max_v_sum[i] /= count_max_v[i];
			max_v_sq_sum[i] /= count_max_v[i];
			st_div_v[i] = sqrt((max_v_sq_sum[i] - max_v_sum[i]*max_v_sum[i])/count_max_v[i]);
		}
	}
	fstream f1;
    stringstream name1;
    name1 << "eq_physical" << "/" << L 
    << "_" << mu_in << "_" << Pre_in << "_" << alpha_in << "_" << k_in 
    << "_" << T_in << "_" << step1 << "_dynamics_sample_ave.txt";
    f1.open(name1.str().c_str(), ios::out);
    for(int i=0; i<num_bin; i++){
    	f1 << ((i+0.5)/(double)num_bin) << " " << max_f_sum[i] << " " << st_div_f[i] << " " << max_v_sum[i] << " " << st_div_f[i] << endl;
    }
}