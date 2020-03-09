

#include <cstdlib>
#include <iomanip>
#include <cmath>
#include "Poisson.h"
Poisson::Poisson(){
	Nx = 0;
	Ny = 0; 
	Pv = {};
	Ps = {};
}
Poisson::~Poisson(){
	
}

void Poisson::SetGridSize(int nx, int ny)
{
	Nx = nx;
	Ny = ny;
}

void Poisson::SetDomainSize(double xlen, double ylen){
	Lx = xlen;
	Ly = ylen;
}

void Poisson::ComputeStreamFunction(double* fi, double* omag,double* A2){
	Pv = omag;
	Ps = fi;
	PoissonMA = A2;
	int N = (Nx-2)*(Ny-2);
	double det_y = Ly/double (Ny-1);
	double det_x = Lx/double (Nx-1);
	//const int N = (Nx-2)*(Ny-2);    // 	Matrix dimension
	const int Ldh = 2*Nx+1;       // leading diamention
	int ld = Nx-2;
	A2[ld   ] = 2*(det_y*det_y+det_x*det_x);
	A2[ld +1] = det_x*det_x;
	A2[Ldh-1] = det_y*det_y;
	for (int i = 1; i < N; ++i) {
		if (i < Nx-2){	
			A2[i*Ldh + Ldh-1] = det_y*det_y;
			A2[i*Ldh + ld +1] = det_x*det_x;
			A2[i*Ldh + ld   ] = 2*(det_y*det_y+det_x*det_x);
			A2[i*Ldh + ld -1] = det_x*det_x;
		}
		else {
			if (i%(Ldh-1) == 0){
				A2[i*Ldh + ld  ] = 2*(det_y*det_y+det_x*det_x);
				A2[i*Ldh       ] = det_y*det_y;	
				A2[i*Ldh + Ldh-1] = det_y*det_y;				
			}
			else {
				A2[i*Ldh + ld   ] = 2*(det_y*det_y+det_x*det_x);	
				A2[i*Ldh + ld -1] = det_x*det_x;
				A2[i*Ldh        ] = det_y*det_y;
				A2[i*Ldh + ld +1] = det_x*det_x;		
				A2[i*Ldh + Ldh-1] = det_y*det_y;				
			}
		}				
    }	
	
}