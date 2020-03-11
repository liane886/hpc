#include <cstdlib>
#include <cmath>
#include "OPoisson.h"
#include <fstream>
#include <iomanip>
#include "/home/li/Desktop/header/cblas/CBLAS/include/cblas.h"
#define F77NAME(x) x##_

// This tells the C++ compiler to use C-style naming of symbols for these
// functions, as opposed to the "mangled" symbols used in C++. The need to
// mangle symbol names is because C++ supports overloading the same function
// name with different parameters and data types so needs to be able to
// distinguish them by name in the object code.
extern "C" {
    // LAPACK routine for solving systems of linear equations8void 
	void F77NAME(dgbsv)(const int& n, const int& KL,
						const int& KU,const int& NRHS, const double* AB,
						const int& ldab, int* ipiv, double* B,
						const int& ldb, int& info);
}
//double cblas_dnrm2(const int N, const double *X, const int incX);

Poisson::Poisson(){
	Nx = 0;
	Ny = 0; 
	Pv = {};
	Ps = {};
	PoissonMA = {};
}
Poisson::~Poisson(){
	cout<<"object is being deleted"<<endl;
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
	cout<<det_x*det_x<<endl;
	//const int N = (Nx-2)*(Ny-2);    // 	Matrix dimension
	const int Ldh = 3*Nx+1;       // leading diamention
	int ld = 2*(Nx-2);
	A2[ld   ] = 2*(det_y*det_y+det_x*det_x);
	A2[ld +1] = det_x*det_x;
	A2[Ldh-1] = det_y*det_y;
	for (int i = 1; i < N; ++i) {
		if (i < Nx-2){	
			//A2[i*Ldh + Ldh-1] = 1;
			A2[i*Ldh + Ldh-1] = - det_y*det_y;
			A2[i*Ldh + ld +1] = - det_x*det_x;
			A2[i*Ldh + ld   ] = - 2*(det_y*det_y+det_x*det_x);
			A2[i*Ldh + ld -1] = - det_x*det_x;
		}
		else  {
			if (i%(Nx-2) == 0){
				A2[i*Ldh + ld  ] = - 2*(det_y*det_y+det_x*det_x);
				A2[i*Ldh + Nx-2] = - det_y*det_y;	
				A2[i*Ldh + Ldh-1] = - det_y*det_y;				
			}
			else {
				A2[i*Ldh + ld   ] = - 2*(det_y*det_y+det_x*det_x);	
				A2[i*Ldh + ld -1] = - det_x*det_x;
				A2[i*Ldh + Nx-2 ] = - det_y*det_y;
				A2[i*Ldh + ld +1] = - det_x*det_x;		
				A2[i*Ldh + Ldh-1] = - det_y*det_y;				
			}
		}				
    }
	
//    cout.precision(4);
//    for (int i = 0; i < Ldh; ++i) {
//        for (int j = 0; j < N; ++j) {
//            cout << setw(6) << A2[j*Ldh+i] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;
	
	int KUL = Nx-2;
	int* ipiv = new int[N];
	int info;
	double* r = new double[N];

	cblas_dcopy(N, omag, 1, r, 1);   // r_0 = b (i.e. omag)
	

	F77NAME(dgbsv)(N,KUL,KUL,1,A2,Ldh,ipiv,r,N,info);	
	cblas_dcopy(N, r, 1, fi, 1);   // fi = x (i.e. r)
	
	if (info != 0){
		cout <<"ERROR: An error occurred in DGBSV:"<<endl;
		cout<<info<<endl;
	}
	

	
//	ofstream stream;
//	stream.open("/home/li/Desktop/hpc/hpc-cw/data.txt",ios::trunc);
//	for(int i =0;i<N;++i){
//		stream<<r[i];
//		stream<<" ";
//	}
//	stream.close();
	
	
	
	
}