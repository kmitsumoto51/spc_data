#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include "mt19937ar.h"
#include "const.h"
#include "utils.h"

using namespace std;

void make_network(int n_b[N][2*edge]){
	for(int i=0; i<N; i++){
		n_b[i][0]=(i/L)*L+(i+1)%L;
		n_b[i][1]=(i+L)%N;
		n_b[i][2]=(i/L)*L+(L+i-1)%L;
		n_b[i][3]=(N+i-L)%N;
	}
	for(int i=0; i<N; i++){
		n_b[i][4]=n_b[n_b[i][1]][0];
		n_b[i][5]=n_b[n_b[i][1]][2];
		n_b[i][6]=n_b[n_b[i][3]][2];
		n_b[i][7]=n_b[n_b[i][3]][0];
	}
}

int num_ads(int sigma[N]){
	int ad_sum = 0;
	for(int i=0; i<N; i++){
		ad_sum += sigma[i];
	}
	return ad_sum;
}

double lattice_distance(double x1, double y1, double x2, double y2){
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

double susceptibility(double m_sq_ave,double m_ave_sq,double beta){
	return beta*(m_sq_ave - m_ave_sq);
}

double specific_heat(double E_sq_ave,double E_ave,double beta){
	return beta*beta*(E_sq_ave - E_ave*E_ave);
}

void local_potential(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu,
	double a1, double p_plaquette[N], double p_plaquette2[N], double p_plaquette3[N], double f_bond[N][4]){
	double posix1,posiy1,posix2,posiy2,posix3,posiy3,posix4,posiy4;
	double l_d1,l_d2,l_d3,l_d4,l_d5,l_d6;
	double size_l = a1*L;
	double l_1 = 1+alpha;
	double l_2 = sqrt(2)*(1+alpha);
	for(int i=0; i<N; i++){
		posix1 = position[i][0];
		posiy1 = position[i][1];
		posix2 = position[n_b[i][0]][0];
		posiy2 = position[n_b[i][0]][1];
		posix3 = position[n_b[i][4]][0];
		posiy3 = position[n_b[i][4]][1];
		posix4 = position[n_b[i][1]][0];
		posiy4 = position[n_b[i][1]][1];
		if(i%L == L-1){
			posix2 = position[n_b[i][0]][0] + size_l;
			posix3 = position[n_b[i][4]][0] + size_l;
		}
		if(i/L == L-1){
			posiy3 = position[n_b[i][4]][1] + size_l;
			posiy4 = position[n_b[i][1]][1] + size_l;
		}
		l_d1 = lattice_distance(posix1, posiy1, posix2, posiy2);
		l_d2 = lattice_distance(posix2, posiy2, posix3, posiy3);
		l_d3 = lattice_distance(posix3, posiy3, posix4, posiy4);
		l_d4 = lattice_distance(posix4, posiy4, posix1, posiy1);
		l_d5 = lattice_distance(posix1, posiy1, posix3, posiy3);
		l_d6 = lattice_distance(posix2, posiy2, posix4, posiy4);

		p_plaquette[i] = 0.25*((1 - l_d1)*(1 - l_d1) + (1 - l_d2)*(1 - l_d2) + (1 - l_d3)*(1 - l_d3) + (1 - l_d4)*(1 - l_d4))
			  +0.5*((sqrt(2) - l_d5)*(sqrt(2) - l_d5) + (sqrt(2) - l_d6)*(sqrt(2) - l_d6))
			  +0.25*k*sigma[i]*((l_1 - l_d1)*(l_1 - l_d1) + (l_1 - l_d2)*(l_1 - l_d2) 
			  				  + (l_1 - l_d3)*(l_1 - l_d3) + (l_1 - l_d4)*(l_1 - l_d4))
			  +0.5*k*sigma[i]*((l_2 - l_d5)*(l_2 - l_d5) + (l_2 - l_d6)*(l_2 - l_d6));
		p_plaquette2[i] = p_plaquette[i] - (mu - 3*k*alpha*alpha/(1.0+k))*sigma[i];
		p_plaquette3[i] = 0.25*((1 - l_d1)*(1 - l_d1) + (1 - l_d2)*(1 - l_d2) + (1 - l_d3)*(1 - l_d3) + (1 - l_d4)*(1 - l_d4))
			  +0.5*((sqrt(2) - l_d5)*(sqrt(2) - l_d5) + (sqrt(2) - l_d6)*(sqrt(2) - l_d6));

		f_bond[i][0] = -(1 - l_d1) - 0.5*k*(sigma[i]+sigma[n_b[i][3]])*(l_1 - l_d1);
		f_bond[i][1] = -(1 - l_d4) - 0.5*k*(sigma[i]+sigma[n_b[i][2]])*(l_1 - l_d4);
		f_bond[i][2] = -(sqrt(2) - l_d5) - k*sigma[i]*(l_2 - l_d5);
		f_bond[i][3] = -(sqrt(2) - l_d6) - k*sigma[i]*(l_2 - l_d6);
	}
}

double local_potential2(double posix[4], double posiy[4], int n_b[N][2*edge], double alpha, double k){
	double l_d1 = lattice_distance(posix[0], posiy[0], posix[1], posiy[1]);
	double l_d2 = lattice_distance(posix[1], posiy[1], posix[2], posiy[2]);
	double l_d3 = lattice_distance(posix[2], posiy[2], posix[3], posiy[3]);
	double l_d4 = lattice_distance(posix[3], posiy[3], posix[0], posiy[0]);
	double l_1 = 1+alpha;
	double p_plaquette2 = 0.25*k*((l_1 - l_d1)*(l_1 - l_d1) + (l_1 - l_d2)*(l_1 - l_d2) 
			  				 	+ (l_1 - l_d3)*(l_1 - l_d3) + (l_1 - l_d4)*(l_1 - l_d4));
	return p_plaquette2;
}

double local_potential2_next(double posix[4], double posiy[4], int n_b[N][2*edge], double alpha, double k){
	double l_d1 = lattice_distance(posix[0], posiy[0], posix[1], posiy[1]);
	double l_d2 = lattice_distance(posix[1], posiy[1], posix[2], posiy[2]);
	double l_d3 = lattice_distance(posix[2], posiy[2], posix[3], posiy[3]);
	double l_d4 = lattice_distance(posix[3], posiy[3], posix[0], posiy[0]);
	double l_d5 = lattice_distance(posix[0], posiy[0], posix[2], posiy[2]);
	double l_d6 = lattice_distance(posix[1], posiy[1], posix[3], posiy[3]);
	double l_1 = 1+alpha;
	double l_2 = sqrt(2)*(1+alpha);
	double p_plaquette2 = 0.25*k*((l_1 - l_d1)*(l_1 - l_d1) + (l_1 - l_d2)*(l_1 - l_d2) 
			  				 	+ (l_1 - l_d3)*(l_1 - l_d3) + (l_1 - l_d4)*(l_1 - l_d4))
			  		     + 0.5*k*((l_2 - l_d5)*(l_2 - l_d5) + (l_2 - l_d6)*(l_2 - l_d6));
	return p_plaquette2;
}

void force_average(double f_bond[N][4], double f_ave[2]){
	double f_sum1 = 0;
	double f_sum2 = 0;
	for(int i=0; i<N; i++){
		f_sum1 += f_bond[i][0] + f_bond[i][1];
		f_sum2 += f_bond[i][2] + f_bond[i][3];
	}
	f_ave[0] = f_sum1/(double)(2*N);
	f_ave[1] = f_sum2/(double)(2*N);
}

void force_correlation(double f_bond[N][4], double f_correlation[L/2][4]){
	for(int j=0; j<L/2; j++){
		f_correlation[j][0] = 0.0;
		f_correlation[j][1] = 0.0;
		f_correlation[j][2] = 0.0;
		f_correlation[j][3] = 0.0;
	}
	for(int i=0; i<N; i++){
		for(int j=0; j<L/2; j++){
			f_correlation[j][0] += f_bond[i][0]*f_bond[(i/L)*L+(i+j)%L][0] + f_bond[i][1]*f_bond[(i+j*L)%N][1];
			f_correlation[j][1] += f_bond[i][1]*f_bond[(i/L)*L+(i+j)%L][1] + f_bond[i][0]*f_bond[(i+j*L)%N][0];
			f_correlation[j][2] += f_bond[i][2]*f_bond[(i/L)*L+(i+j)%L][2] + f_bond[i][2]*f_bond[(i+j*L)%N][2];
			f_correlation[j][3] += f_bond[i][0]*f_bond[(i/L)*L+(i+j)%L][1] + f_bond[i][1]*f_bond[(i+j*L)%N][0];
		}
	}
	for(int j=0; j<L/2; j++){
		f_correlation[j][0] /= (double)(2*N);
		f_correlation[j][1] /= (double)(2*N);
		f_correlation[j][2] /= (double)(2*N);
		f_correlation[j][3] /= (double)(2*N);
	}
}

void particle_correlation(int sigma[N], double ads_correlation[L/2]){
	for(int j=0; j<L/2; j++){
		ads_correlation[j] = 0.0;
	}
	for(int i=0; i<N; i++){
		for(int j=0; j<L/2; j++){
			ads_correlation[j] += sigma[i]*(sigma[(i/L)*L+(i+j)%L] + sigma[(i+j*L)%N]);
		}
	}
	for(int j=0; j<L/2; j++){
		ads_correlation[j] /= (double)(2*N);
	}
}

void plaquette_volume(double position[N][2], int n_b[N][2*edge], double a1, double p_vol[N]){
	double posix1,posiy1,posix2,posiy2,posix3,posiy3,posix4,posiy4;
	double size_l = a1*L;
	for(int i=0; i<N; i++){
		posix1 = position[i][0];
		posiy1 = position[i][1];
		posix2 = position[n_b[i][0]][0];
		posiy2 = position[n_b[i][0]][1];
		posix3 = position[n_b[i][4]][0];
		posiy3 = position[n_b[i][4]][1];
		posix4 = position[n_b[i][1]][0];
		posiy4 = position[n_b[i][1]][1];
		if(i%L == L-1){
			posix2 = position[n_b[i][0]][0] + size_l;
			posix3 = position[n_b[i][4]][0] + size_l;
		}
		if(i/L == L-1){
			posiy3 = position[n_b[i][4]][1] + size_l;
			posiy4 = position[n_b[i][1]][1] + size_l;
		}
		p_vol[i] = abs(0.5*(posix1*posiy2 - posix2*posiy1 
						   +posix2*posiy3 - posix3*posiy2
						   +posix3*posiy4 - posix4*posiy3
						   +posix4*posiy1 - posix1*posiy4));
	}
}

void contact_list(int sigma[N], int n_b[N][2*edge], int contact_f[N][edge], int contact_v[N][edge]){
	for(int i=0; i<N; i++){
		for(int j=0; j<edge; j++){
			contact_f[i][j] = 0;
			contact_v[i][j] = 0;
		}
	}
	for(int i=0; i<N; i++){
		if(sigma[i]==1){
			for(int j=0; j<edge; j++){
				if(sigma[n_b[i][j]]==1){
					contact_f[i][j]=1;
				}
			}
			
		}
		if(sigma[i]==0){
			for(int j=0; j<edge; j++){
				if(sigma[n_b[i][j]]==0){
					contact_v[i][j]=1;
				}
			}
		}
	}
}

int site_percolation_f(int sigma[N], int n_b[N][2*edge], int cluster_num_f[N], int &max_cluster_size_f, int contact_f[N][edge]){
	int *bcn;
	bcn = new int[3];
	int *bl;
	bl = new int[3];
	int *bl_low;
	bl_low = new int[2];
	int *label_num;
	label_num = new int[N];
	int perc = 0;
	int new_num = 0;
	for(int i=0; i<N; i++){
		cluster_num_f[i] = 0;
		label_num[i] = 0;
	}
	if(sigma[0]==0){
		cluster_num_f[0] = N;
	}else{
		cluster_num_f[0]=new_num;
		label_num[new_num]=new_num;
	}
	for(int i=0; i<N; i++){
		if(sigma[i]==0){
			cluster_num_f[i] = N;
		}else{
			for(int j=0; j<3; j++){
				bcn[j] = N;
				bl[j] = N;
			}
			if(i%L != 0){
				if(contact_f[i][2]==1){
					bcn[0] = cluster_num_f[n_b[i][2]];
					bl[0]=label_num[bcn[0]];
				}
			}
			if(i > L-1){
				if(contact_f[i][3]==1){
					bcn[1] = cluster_num_f[n_b[i][3]];
					bl[1]=label_num[bcn[1]];
				}
			}
			if(i%L == L-1){
				if(contact_f[i][0]==1){
					bcn[2] = cluster_num_f[n_b[i][0]];
					bl[2]=label_num[bcn[2]];
				}
			}
			int flg = 0;
			for(int j=0; j<3; j++){
				if(bcn[j] != N){
					flg = 1;
				}
			}
			if(flg == 0){
				new_num += 1;
				cluster_num_f[i]=new_num;
				label_num[new_num]=new_num;
			}else{
				int bl_min=N;
	        	for(int j=0; j<3; j++){
					if(bl[j] < bl_min){
						bl_min = bl[j];
					}
				}
				int j2 = 0;
				bl_low[0] = N;
				bl_low[1] = N;
				for(int j=0; j<3; j++){
					if(bl[j] != bl_min && bl[j] != N){
						bl_low[j2] = bl[j];
						j2+=1;
					}
				}
				cluster_num_f[i]=bl_min;
				for(int j=0; j<2; j++){
					if(bl_low[j] != N){
						label_num[bl_low[j]]=label_num[bl_min];
						for(int k=0; k<new_num+1; ++k){
		                    if(label_num[k]==bl_low[j]){
		                        label_num[k]=bl_min;
		                    }
		                }
					}
				}
			}
		}
	}
	int *cluster_size;
	cluster_size = new int[new_num+1];
    for(int k=0; k<new_num+1; k++){
        cluster_size[k] =0;
    }
    int max_size =-1;
    for(int k=0; k<new_num+1; k++){
        for(int i=0; i<N; i++){
            if(cluster_num_f[i]==k){
                cluster_num_f[i]=label_num[k];
                cluster_size[label_num[k]] +=1;
            }
        }
    }
    for(int i=0; i<L; i++){
        for(int j=0; j<L; j++){
            if(cluster_num_f[i]==cluster_num_f[N-L+j] && cluster_num_f[i] != N){
                perc=1;
            }
        }
    }
    for(int k=0; k<new_num+1; k++){
        if(cluster_size[k]>max_size){
            max_size=cluster_size[k];
        }
    }
    if(max_size==-1){
        max_cluster_size_f = 0;
    }else{
        max_cluster_size_f = max_size;
    }
	delete[] label_num;
	delete[] cluster_size;
	delete[] bcn;
	delete[] bl;
	delete[] bl_low;
	return perc;
}

void cluster_ratio_f(int sigma[N], int n_b[N][2*edge], int cluster_num_f[N], double c_ratio_f[2], int contact_f[N][edge]){
	int *bcn;
	bcn = new int[4];
	int *bl;
	bl = new int[4];
	int *bl_low;
	bl_low = new int[3];
	int *label_num;
	label_num = new int[N];
	int new_num = 0;
	for(int i=0; i<N; i++){
		cluster_num_f[i] = 0;
		label_num[i] = 0;
	}
	if(sigma[0]==0){
		cluster_num_f[0] = N;
	}else{
		cluster_num_f[0]=new_num;
		label_num[new_num]=new_num;
	}
	for(int i=0; i<N; i++){
		if(sigma[i]==0){
			cluster_num_f[i] = N;
		}else{
			for(int j=0; j<4; j++){
				bcn[j] = N;
				bl[j] = N;
			}
			if(i%L != 0){
				if(contact_f[i][2]==1){
					bcn[0] = cluster_num_f[n_b[i][2]];
					bl[0]=label_num[bcn[0]];
				}
			}
			if(i > L-1){
				if(contact_f[i][3]==1){
					bcn[1] = cluster_num_f[n_b[i][3]];
					bl[1]=label_num[bcn[1]];
				}
			}
			if(i%L == L-1){
				if(contact_f[i][0]==1){
					bcn[2] = cluster_num_f[n_b[i][0]];
					bl[2]=label_num[bcn[2]];
				}
			}
			if(i/L == L-1){
				if(contact_f[i][1]==1){
					bcn[3] = cluster_num_f[n_b[i][1]];
					bl[3]=label_num[bcn[3]];
				}
			}
			int flg = 0;
			for(int j=0; j<4; j++){
				if(bcn[j] != N){
					flg = 1;
				}
			}
			if(flg == 0){
				new_num += 1;
				cluster_num_f[i]=new_num;
				label_num[new_num]=new_num;
			}else{
				int bl_min=N;
	        	for(int j=0; j<4; j++){
					if(bl[j] < bl_min){
						bl_min = bl[j];
					}
				}
				int j2 = 0;
				bl_low[0] = N;
				bl_low[1] = N;
				bl_low[2] = N;
				for(int j=0; j<4; j++){
					if(bl[j] != bl_min && bl[j] != N){
						bl_low[j2] = bl[j];
						j2+=1;
					}
				}
				cluster_num_f[i]=bl_min;
				for(int j=0; j<3; j++){
					if(bl_low[j] != N){
						label_num[bl_low[j]]=label_num[bl_min];
						for(int k=0; k<new_num+1; ++k){
		                    if(label_num[k]==bl_low[j]){
		                        label_num[k]=bl_min;
		                    }
		                }
					}
				}
			}
		}
	}
	int *cluster_size;
	cluster_size = new int[new_num+1];
	int *cluster_circum;
	cluster_circum = new int[new_num+1];
    for(int k=0; k<new_num+1; k++){
        cluster_size[k] =0;
    }
    for(int k=0; k<new_num+1; k++){
        for(int i=0; i<N; i++){
            if(cluster_num_f[i]==k){
                cluster_num_f[i]=label_num[k];
                cluster_size[label_num[k]] +=1;
            }
        }
    }
    int count_f1 = 0;
    int count_f2 = 0;
    for(int k=0; k<new_num+1; k++){
    	if(cluster_size[k] != 0){
    		count_f1 += 1;
    		if(cluster_size[k] > 3){
    			count_f2 += 1;
    		}
    		cluster_circum[k] = 4*cluster_size[k];
    		for(int i=0; i<N; i++){
    			if(cluster_num_f[i]==k){
    				for(int j=0; j<4; j++){
    					cluster_circum[k] -= contact_f[i][j];
    				}
    			}
    		}
    	}
    }
    double sum_ratio1 = 0.0;
    double sum_ratio2 = 0.0;
    for(int k=0; k<new_num+1; k++){
    	if(cluster_size[k] != 0){
    		sum_ratio1 += cluster_circum[k]/(double)cluster_size[k];
    	}
    	if(cluster_size[k] > 3){
    		sum_ratio2 += cluster_circum[k]/(double)cluster_size[k];
    	}
    }
    if(count_f1 != 0){
    	c_ratio_f[0] = sum_ratio1/count_f1;
    }else{
    	c_ratio_f[0] = 0;
    }
    if(count_f2 != 0){
    	c_ratio_f[1] = sum_ratio2/count_f2;
    }else{
    	c_ratio_f[1] = 0;
    }

	delete[] label_num;
	delete[] cluster_size;
	delete[] bcn;
	delete[] bl;
	delete[] bl_low;
	delete[] cluster_circum;
}


int site_percolation_v(int sigma[N], int n_b[N][2*edge], int cluster_num_v[N], int &max_cluster_size_v, int contact_v[N][edge]){
	int *bcn;
	bcn = new int[3];
	int *bl;
	bl = new int[3];
	int *bl_low;
	bl_low = new int[2];
	int *label_num;
	label_num = new int[N];
	int perc = 0;
	int new_num = 0;
	for(int i=0; i<N; i++){
		cluster_num_v[i] = 0;
		label_num[i] = 0;
	}
	if(sigma[0]==1){
		cluster_num_v[0] = N;
	}else{
		cluster_num_v[0]=new_num;
		label_num[new_num]=new_num;
	}
	for(int i=0; i<N; i++){
		if(sigma[i]==1){
			cluster_num_v[i] = N;
		}else{
			for(int j=0; j<3; j++){
				bcn[j] = N;
				bl[j] = N;
			}
			if(i%L != 0){
				if(contact_v[i][2]==1){
					bcn[0] = cluster_num_v[n_b[i][2]];
					bl[0]=label_num[bcn[0]];
				}
			}
			if(i > L-1){
				if(contact_v[i][3]==1){
					bcn[1] = cluster_num_v[n_b[i][3]];
					bl[1]=label_num[bcn[1]];
				}
			}
			if(i%L == L-1){
				if(contact_v[i][0]==1){
					bcn[2] = cluster_num_v[n_b[i][0]];
					bl[2]=label_num[bcn[2]];
				}
			}
			int flg = 0;
			for(int j=0; j<3; j++){
				if(bcn[j] != N){
					flg = 1;
				}
			}
			if(flg == 0){
				new_num += 1;
				cluster_num_v[i]=new_num;
				label_num[new_num]=new_num;
			}else{
				int bl_min=N;
	        	for(int j=0; j<3; j++){
					if(bl[j] < bl_min){
						bl_min = bl[j];
					}
				}
				int j2 = 0;
				bl_low[0] = N;
				bl_low[1] = N;
				for(int j=0; j<3; j++){
					if(bl[j] != bl_min && bl[j] != N){
						bl_low[j2] = bl[j];
						j2+=1;
					}
				}
				cluster_num_v[i]=bl_min;
				for(int j=0; j<2; j++){
					if(bl_low[j] != N){
						label_num[bl_low[j]]=label_num[bl_min];
						for(int k=0; k<new_num+1; ++k){
		                    if(label_num[k]==bl_low[j]){
		                        label_num[k]=bl_min;
		                    }
		                }
					}
				}
			}
		}
	}
	int *cluster_size;
	cluster_size = new int[new_num+1];
    for(int k=0; k<new_num+1; k++){
        cluster_size[k] =0;
    }
    int max_size =-1;
    for(int k=0; k<new_num+1; k++){
        for(int i=0; i<N; i++){
            if(cluster_num_v[i]==k){
                cluster_num_v[i]=label_num[k];
                cluster_size[label_num[k]] +=1;
            }
        }
    }
    for(int i=0; i<L; i++){
        for(int j=0; j<L; j++){
            if(cluster_num_v[i]==cluster_num_v[N-L+j] && cluster_num_v[i] != N){
                perc=1;
            }
        }
    }
    for(int k=0; k<new_num+1; k++){
        if(cluster_size[k]>max_size){
            max_size=cluster_size[k];
        }
    }
    if(max_size==-1){
        max_cluster_size_v = 0;
    }else{
        max_cluster_size_v = max_size;
    }
	delete[] label_num;
	delete[] cluster_size;
	delete[] bcn;
	delete[] bl;
	delete[] bl_low;
	return perc;
}

void cluster_ratio_v(int sigma[N], int n_b[N][2*edge], int cluster_num_v[N], double c_ratio_v[2], int contact_v[N][edge]){
	int *bcn;
	bcn = new int[4];
	int *bl;
	bl = new int[4];
	int *bl_low;
	bl_low = new int[3];
	int *label_num;
	label_num = new int[N];
	int new_num = 0;
	for(int i=0; i<N; i++){
		cluster_num_v[i] = 0;
		label_num[i] = 0;
	}
	if(sigma[0]==1){
		cluster_num_v[0] = N;
	}else{
		cluster_num_v[0]=new_num;
		label_num[new_num]=new_num;
	}
	for(int i=0; i<N; i++){
		if(sigma[i]==1){
			cluster_num_v[i] = N;
		}else{
			for(int j=0; j<4; j++){
				bcn[j] = N;
				bl[j] = N;
			}
			if(i%L != 0){
				if(contact_v[i][2]==1){
					bcn[0] = cluster_num_v[n_b[i][2]];
					bl[0]=label_num[bcn[0]];
				}
			}
			if(i > L-1){
				if(contact_v[i][3]==1){
					bcn[1] = cluster_num_v[n_b[i][3]];
					bl[1]=label_num[bcn[1]];
				}
			}
			if(i%L == L-1){
				if(contact_v[i][0]==1){
					bcn[2] = cluster_num_v[n_b[i][0]];
					bl[2]=label_num[bcn[2]];
				}
			}
			if(i/L == L-1){
				if(contact_v[i][1]==1){
					bcn[3] = cluster_num_v[n_b[i][1]];
					bl[3]=label_num[bcn[3]];
				}
			}
			int flg = 0;
			for(int j=0; j<4; j++){
				if(bcn[j] != N){
					flg = 1;
				}
			}
			if(flg == 0){
				new_num += 1;
				cluster_num_v[i]=new_num;
				label_num[new_num]=new_num;
			}else{
				int bl_min=N;
	        	for(int j=0; j<4; j++){
					if(bl[j] < bl_min){
						bl_min = bl[j];
					}
				}
				int j2 = 0;
				bl_low[0] = N;
				bl_low[1] = N;
				bl_low[2] = N;
				for(int j=0; j<4; j++){
					if(bl[j] != bl_min && bl[j] != N){
						bl_low[j2] = bl[j];
						j2+=1;
					}
				}
				cluster_num_v[i]=bl_min;
				for(int j=0; j<3; j++){
					if(bl_low[j] != N){
						label_num[bl_low[j]]=label_num[bl_min];
						for(int k=0; k<new_num+1; ++k){
		                    if(label_num[k]==bl_low[j]){
		                        label_num[k]=bl_min;
		                    }
		                }
					}
				}
			}
		}
	}
	int *cluster_size;
	cluster_size = new int[new_num+1];
	int *cluster_circum;
	cluster_circum = new int[new_num+1];
    for(int k=0; k<new_num+1; k++){
        cluster_size[k] =0;
    }
    for(int k=0; k<new_num+1; k++){
        for(int i=0; i<N; i++){
            if(cluster_num_v[i]==k){
                cluster_num_v[i]=label_num[k];
                cluster_size[label_num[k]] +=1;
            }
        }
    }
    int count_f1 = 0;
    int count_f2 = 0;
    for(int k=0; k<new_num+1; k++){
    	if(cluster_size[k] != 0){
    		count_f1 += 1;
    		if(cluster_size[k] > 3){
    			count_f2 += 1;
    		}
    		cluster_circum[k] = 4*cluster_size[k];
    		for(int i=0; i<N; i++){
    			if(cluster_num_v[i]==k){
    				for(int j=0; j<4; j++){
    					cluster_circum[k] -= contact_v[i][j];
    				}
    			}
    		}
    	}
    }
    double sum_ratio1 = 0.0;
    double sum_ratio2 = 0.0;
    for(int k=0; k<new_num+1; k++){
    	if(cluster_size[k] != 0){
    		sum_ratio1 += cluster_circum[k]/(double)cluster_size[k];
    	}
    	if(cluster_size[k] > 3){
    		sum_ratio2 += cluster_circum[k]/(double)cluster_size[k];
    	}
    }
    if(count_f1 != 0){
    	c_ratio_v[0] = sum_ratio1/count_f1;
    }else{
    	c_ratio_v[0] = 0;
    }
    if(count_f2 != 0){
    	c_ratio_v[1] = sum_ratio2/count_f2;
    }else{
    	c_ratio_v[1] = 0;
    }

	delete[] label_num;
	delete[] cluster_size;
	delete[] bcn;
	delete[] bl;
	delete[] bl_low;
	delete[] cluster_circum;
}

void cluster_ratio_hist_f(int sigma[N], int n_b[N][2*edge], int cluster_num_f[N], int contact_f[N][edge], int size_hist_f[N-1], double ratio_hist_f[N-1]){
	int *bcn;
	bcn = new int[4];
	int *bl;
	bl = new int[4];
	int *bl_low;
	bl_low = new int[3];
	int *label_num;
	label_num = new int[N];
	int new_num = 0;
	for(int i=0; i<N; i++){
		cluster_num_f[i] = 0;
		label_num[i] = 0;
	}
	if(sigma[0]==0){
		cluster_num_f[0] = N;
	}else{
		cluster_num_f[0]=new_num;
		label_num[new_num]=new_num;
	}
	for(int i=0; i<N; i++){
		if(sigma[i]==0){
			cluster_num_f[i] = N;
		}else{
			for(int j=0; j<4; j++){
				bcn[j] = N;
				bl[j] = N;
			}
			if(i%L != 0){
				if(contact_f[i][2]==1){
					bcn[0] = cluster_num_f[n_b[i][2]];
					bl[0]=label_num[bcn[0]];
				}
			}
			if(i > L-1){
				if(contact_f[i][3]==1){
					bcn[1] = cluster_num_f[n_b[i][3]];
					bl[1]=label_num[bcn[1]];
				}
			}
			if(i%L == L-1){
				if(contact_f[i][0]==1){
					bcn[2] = cluster_num_f[n_b[i][0]];
					bl[2]=label_num[bcn[2]];
				}
			}
			if(i/L == L-1){
				if(contact_f[i][1]==1){
					bcn[3] = cluster_num_f[n_b[i][1]];
					bl[3]=label_num[bcn[3]];
				}
			}
			int flg = 0;
			for(int j=0; j<4; j++){
				if(bcn[j] != N){
					flg = 1;
				}
			}
			if(flg == 0){
				new_num += 1;
				cluster_num_f[i]=new_num;
				label_num[new_num]=new_num;
			}else{
				int bl_min=N;
	        	for(int j=0; j<4; j++){
					if(bl[j] < bl_min){
						bl_min = bl[j];
					}
				}
				int j2 = 0;
				bl_low[0] = N;
				bl_low[1] = N;
				bl_low[2] = N;
				for(int j=0; j<4; j++){
					if(bl[j] != bl_min && bl[j] != N){
						bl_low[j2] = bl[j];
						j2+=1;
					}
				}
				cluster_num_f[i]=bl_min;
				for(int j=0; j<3; j++){
					if(bl_low[j] != N){
						label_num[bl_low[j]]=label_num[bl_min];
						for(int k=0; k<new_num+1; ++k){
		                    if(label_num[k]==bl_low[j]){
		                        label_num[k]=bl_min;
		                    }
		                }
					}
				}
			}
		}
	}
	int *cluster_size;
	cluster_size = new int[new_num+1];
	int *cluster_circum;
	cluster_circum = new int[new_num+1];
    for(int k=0; k<new_num+1; k++){
        cluster_size[k] =0;
    }
    for(int k=0; k<new_num+1; k++){
        for(int i=0; i<N; i++){
            if(cluster_num_f[i]==k){
                cluster_num_f[i]=label_num[k];
                cluster_size[label_num[k]] +=1;
            }
        }
    }
    for(int k=0; k<new_num+1; k++){
    	for(int j=0; j<N; j++){
			if(cluster_size[k] == j+1){
				size_hist_f[j] += 1;
			}
		}
    	if(cluster_size[k] != 0){
    		cluster_circum[k] = 4*cluster_size[k];
    		for(int i=0; i<N; i++){
    			if(cluster_num_f[i]==k){
    				for(int j=0; j<4; j++){
    					cluster_circum[k] -= contact_f[i][j];
    				}
    			}
    		}
    	}
    }
    for(int k=0; k<new_num+1; k++){
    	for(int j=0; j<N; j++){
			if(cluster_size[k] == j+1){
				ratio_hist_f[j] += cluster_circum[k]/(double)cluster_size[k];
			}
		}
    }

	delete[] label_num;
	delete[] cluster_size;
	delete[] bcn;
	delete[] bl;
	delete[] bl_low;
	delete[] cluster_circum;
}

void cluster_ratio_hist_v(int sigma[N], int n_b[N][2*edge], int cluster_num_v[N], int contact_v[N][edge], int size_hist_v[N-1], double ratio_hist_v[N-1]){
	int *bcn;
	bcn = new int[4];
	int *bl;
	bl = new int[4];
	int *bl_low;
	bl_low = new int[3];
	int *label_num;
	label_num = new int[N];
	int new_num = 0;
	for(int i=0; i<N; i++){
		cluster_num_v[i] = 0;
		label_num[i] = 0;
	}
	if(sigma[0]==1){
		cluster_num_v[0] = N;
	}else{
		cluster_num_v[0]=new_num;
		label_num[new_num]=new_num;
	}
	for(int i=0; i<N; i++){
		if(sigma[i]==1){
			cluster_num_v[i] = N;
		}else{
			for(int j=0; j<4; j++){
				bcn[j] = N;
				bl[j] = N;
			}
			if(i%L != 0){
				if(contact_v[i][2]==1){
					bcn[0] = cluster_num_v[n_b[i][2]];
					bl[0]=label_num[bcn[0]];
				}
			}
			if(i > L-1){
				if(contact_v[i][3]==1){
					bcn[1] = cluster_num_v[n_b[i][3]];
					bl[1]=label_num[bcn[1]];
				}
			}
			if(i%L == L-1){
				if(contact_v[i][0]==1){
					bcn[2] = cluster_num_v[n_b[i][0]];
					bl[2]=label_num[bcn[2]];
				}
			}
			if(i/L == L-1){
				if(contact_v[i][1]==1){
					bcn[3] = cluster_num_v[n_b[i][1]];
					bl[3]=label_num[bcn[3]];
				}
			}
			int flg = 0;
			for(int j=0; j<4; j++){
				if(bcn[j] != N){
					flg = 1;
				}
			}
			if(flg == 0){
				new_num += 1;
				cluster_num_v[i]=new_num;
				label_num[new_num]=new_num;
			}else{
				int bl_min=N;
	        	for(int j=0; j<4; j++){
					if(bl[j] < bl_min){
						bl_min = bl[j];
					}
				}
				int j2 = 0;
				bl_low[0] = N;
				bl_low[1] = N;
				bl_low[2] = N;
				for(int j=0; j<4; j++){
					if(bl[j] != bl_min && bl[j] != N){
						bl_low[j2] = bl[j];
						j2+=1;
					}
				}
				cluster_num_v[i]=bl_min;
				for(int j=0; j<3; j++){
					if(bl_low[j] != N){
						label_num[bl_low[j]]=label_num[bl_min];
						for(int k=0; k<new_num+1; ++k){
		                    if(label_num[k]==bl_low[j]){
		                        label_num[k]=bl_min;
		                    }
		                }
					}
				}
			}
		}
	}
	int *cluster_size;
	cluster_size = new int[new_num+1];
	int *cluster_circum;
	cluster_circum = new int[new_num+1];
    for(int k=0; k<new_num+1; k++){
        cluster_size[k] =0;
    }
    for(int k=0; k<new_num+1; k++){
        for(int i=0; i<N; i++){
            if(cluster_num_v[i]==k){
                cluster_num_v[i]=label_num[k];
                cluster_size[label_num[k]] +=1;
            }
        }
    }
    for(int k=0; k<new_num+1; k++){
    	for(int j=0; j<N; j++){
			if(cluster_size[k] == j+1){
				size_hist_v[j] += 1;
			}
		}
    	if(cluster_size[k] != 0){
    		cluster_circum[k] = 4*cluster_size[k];
    		for(int i=0; i<N; i++){
    			if(cluster_num_v[i]==k){
    				for(int j=0; j<4; j++){
    					cluster_circum[k] -= contact_v[i][j];
    				}
    			}
    		}
    	}
    }
    for(int k=0; k<new_num+1; k++){
    	for(int j=0; j<N; j++){
			if(cluster_size[k] == j+1){
				ratio_hist_v[j] += cluster_circum[k]/(double)cluster_size[k];
			}
		}
    }

	delete[] label_num;
	delete[] cluster_size;
	delete[] bcn;
	delete[] bl;
	delete[] bl_low;
	delete[] cluster_circum;
}