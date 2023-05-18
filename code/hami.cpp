#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include "const.h"
#include "utils.h"

using namespace std;

double enthalpy(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double a1, double Pre){
	double ene = 0;
	double posix1,posiy1,posix2,posiy2,posix3,posiy3,posix4,posiy4;
	double l_d1,l_d2,l_d3,l_d4;
	double size_l = L*a1;
	double l_1 = 1+alpha;
	for(int i=0; i<N; i++){
		posix1 = position[i][0];
		posiy1 = position[i][1];
		posiy2 = position[n_b[i][0]][1];
		posix4 = position[n_b[i][1]][0];
		if(i%L == L-1){
			posix2 = position[n_b[i][0]][0] + size_l;
			posix3 = position[n_b[i][4]][0] + size_l;
			if(i/L == L-1){
				posiy3 = position[n_b[i][4]][1] + size_l;
				posiy4 = position[n_b[i][1]][1] + size_l;
			}else{
				posiy3 = position[n_b[i][4]][1];
				posiy4 = position[n_b[i][1]][1];
			}
		}else{
			posix2 = position[n_b[i][0]][0];
			posix3 = position[n_b[i][4]][0];
			if(i/L == L-1){
				posiy3 = position[n_b[i][4]][1] + size_l;
				posiy4 = position[n_b[i][1]][1] + size_l;
			}else{
				posiy3 = position[n_b[i][4]][1];
				posiy4 = position[n_b[i][1]][1];
			}
		}
		l_d1 = lattice_distance(posix1, posiy1, posix2, posiy2);
		l_d2 = lattice_distance(posix2, posiy2, posix3, posiy3);
		l_d3 = lattice_distance(posix3, posiy3, posix4, posiy4);
		l_d4 = lattice_distance(posix4, posiy4, posix1, posiy1);
		
		ene += 0.25*((1 - l_d1)*(1 - l_d1) + (1 - l_d2)*(1 - l_d2) + (1 - l_d3)*(1 - l_d3) + (1 - l_d4)*(1 - l_d4))
			  +0.25*k*sigma[i]*((l_1 - l_d1)*(l_1 - l_d1) + (l_1 - l_d2)*(l_1 - l_d2) 
			  				  + (l_1 - l_d3)*(l_1 - l_d3) + (l_1 - l_d4)*(l_1 - l_d4));
		ene -= mu*sigma[i];
	}
	ene += Pre*size_l*size_l;
	return ene;
}

double enthalpy_next(int sigma[N], double position[N][2], int n_b[N][2*edge], double alpha, double k, double mu, double a1, double Pre){
	double ene = 0;
	double posix1,posiy1,posix2,posiy2,posix3,posiy3,posix4,posiy4;
	double l_d1,l_d2,l_d3,l_d4,l_d5,l_d6;
	double size_l = L*a1;
	double l_1 = 1+alpha;
	double l_2 = sqrt(2)*(1+alpha);
	for(int i=0; i<N; i++){
		posix1 = position[i][0];
		posiy1 = position[i][1];
		posiy2 = position[n_b[i][0]][1];
		posix4 = position[n_b[i][1]][0];
		if(i%L == L-1){
			posix2 = position[n_b[i][0]][0] + size_l;
			posix3 = position[n_b[i][4]][0] + size_l;
			if(i/L == L-1){
				posiy3 = position[n_b[i][4]][1] + size_l;
				posiy4 = position[n_b[i][1]][1] + size_l;
			}else{
				posiy3 = position[n_b[i][4]][1];
				posiy4 = position[n_b[i][1]][1];
			}
		}else{
			posix2 = position[n_b[i][0]][0];
			posix3 = position[n_b[i][4]][0];
			if(i/L == L-1){
				posiy3 = position[n_b[i][4]][1] + size_l;
				posiy4 = position[n_b[i][1]][1] + size_l;
			}else{
				posiy3 = position[n_b[i][4]][1];
				posiy4 = position[n_b[i][1]][1];
			}
		}
		l_d1 = lattice_distance(posix1, posiy1, posix2, posiy2);
		l_d2 = lattice_distance(posix2, posiy2, posix3, posiy3);
		l_d3 = lattice_distance(posix3, posiy3, posix4, posiy4);
		l_d4 = lattice_distance(posix4, posiy4, posix1, posiy1);
		l_d5 = lattice_distance(posix1, posiy1, posix3, posiy3);
		l_d6 = lattice_distance(posix2, posiy2, posix4, posiy4);
		
		
		ene += 0.25*((1 - l_d1)*(1 - l_d1) + (1 - l_d2)*(1 - l_d2) + (1 - l_d3)*(1 - l_d3) + (1 - l_d4)*(1 - l_d4))
			  +0.5*((sqrt(2) - l_d5)*(sqrt(2) - l_d5) + (sqrt(2) - l_d6)*(sqrt(2) - l_d6))
			  +0.25*k*sigma[i]*((l_1 - l_d1)*(l_1 - l_d1) + (l_1 - l_d2)*(l_1 - l_d2) 
			  				  + (l_1 - l_d3)*(l_1 - l_d3) + (l_1 - l_d4)*(l_1 - l_d4))
			  +0.5*k*sigma[i]*((l_2 - l_d5)*(l_2 - l_d5) + (l_2 - l_d6)*(l_2 - l_d6));
		ene -= mu*sigma[i];
	}
	ene += Pre*size_l*size_l;
	return ene;
}