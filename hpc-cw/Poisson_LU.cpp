#include <cstdlib>
#include <cmath>
#include "Poisson.h"
#include <fstream>
#include <iomanip>
#include <cstring>

#include "/home/li/Desktop/header/cblas/CBLAS/include/cblas.h"


// This tells the C++ compiler to use C-style naming of symbols for these
// functions, as opposed to the "mangled" symbols used in C++. The need to
// mangle symbol names is because C++ supports overloading the same function
// name with different parameters and data types so needs to be able to
// distinguish them by name in the object code.

#define F77NAME(x) x##_
extern "C" { 
	
	double F77NAME(dcopy) (const int& n,
                          const double *x, const int& incx,
                          const double *y, const int& incy);
						  
	double F77NAME(daxpy) (const int& n, const double& alpha,
                          const double *x, const int& incx,
                          const double *y, const int& incy);
                                                             
	double F77NAME(dpbtrf) (const char& UPLO, const int& n,
                            const int& kd, const double* AB,
                            const int& ldab, int& info);

    double F77NAME(dpbtrs) (const char& UPLO, const int& n,
                            const int& kd, const int& nRHS,
                            const double* AP, const int& ldab,
                            const double* b, const int& ldb,
                            int& info);


    
}

Poisson::Poisson(){

}
Poisson::~Poisson(){
	cout<<"object is being deleted"<<endl;
	delete[]A;
	
}

void Poisson::PInitialise(int& Nx, int& Ny)
{
	this->Nx = Nx-2;
	this->Ny = Ny-2;
	this->N = (this->Nx)*(this->Ny);
	det_y = Ly/double (Ny+1);
	det_x = Lx/double (Nx+1);
	this->beta[0] = 2*(det_x*det_x+det_y*det_y);
	this->beta[1] = -det_y*det_y;
	this->beta[2] = -det_x*det_x;
	A = new double[N*(this->Ny+1)]{};
	
	this ->buildMA();
	this ->PPTRF();
}

void Poisson::SetDomainSize(double xlen, double ylen){
	Lx = xlen;
	Ly = ylen;
}

void Poisson::buildMA(){

	//generate matrix A;

	
int N = (Nx-2)*(Ny-2);
	double det_y = Ly/double (Ny-1);
	double det_x = Lx/double (Nx-1);
	//cout<<det_x*det_x<<endl;
	//const int N = (Nx-2)*(Ny-2);    // 	Matrix dimension
	const int Ldh = 3*Nx+1;       // leading diamention
	int ld = 2*(Nx-2);
	A[ld   ] = 2*(det_y*det_y+det_x*det_x);
	A[ld +1] = det_x*det_x;
	A[Ldh-1] = det_y*det_y;
	for (int i = 1; i < N; ++i) {
		if (i < Nx-2){	
			//A2[i*Ldh + Ldh-1] = 1;
			A[i*Ldh + Ldh-1] = - det_y*det_y;
			A[i*Ldh + ld +1] = - det_x*det_x;
			A[i*Ldh + ld   ] = - 2*(det_y*det_y+det_x*det_x);
			A[i*Ldh + ld -1] = - det_x*det_x;
		}
		else  {
			if (i%(Nx-2) == 0){
				A[i*Ldh + ld  ] = - 2*(det_y*det_y+det_x*det_x);
				A[i*Ldh + Nx-2] = - det_y*det_y;	
				A[i*Ldh + Ldh-1] = - det_y*det_y;				
			}
			else {
				A[i*Ldh + ld   ] = - 2*(det_y*det_y+det_x*det_x);	
				A[i*Ldh + ld -1] = - det_x*det_x;
				A[i*Ldh + Nx-2 ] = - det_y*det_y;
				A[i*Ldh + ld +1] = - det_x*det_x;		
				A[i*Ldh + Ldh-1] = - det_y*det_y;				
			}
		}				
    }
   
///////////////////////////////////////////////////IO 	
//    cout.precision(4);
//    for (int i = 0; i < Ldh; ++i) {
//        for (int j = 0; j < N; ++j) {
//            cout << setw(6) << A2[j*Ldh+i] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;

	
}

void Poisson::PPTRF(){
	int info;
	
	F77NAME(dpbtrf)('U', N,Ny,A,Ny+1,info);
	if (info != 0){
		cout <<"ERROR: An error occurred in dpptrf:";
		cout<<info<<endl;
	}

}

void Poisson::SolvePoisson(double* s,double*v){

	
	double* vCopy = new double[N]{};  // update vorticity vector	     
	fill_n(vCopy,N,0.0);
	int counter = Ny + 2;
	int KUL = Nx-2;
	int* ipiv = new int[N];
	int info;
	
	cblas_dcopy(N, v, 1, vCopy, 1);   // r_0 = b (i.e. omag)
	

	F77NAME(dgbsv)(N,KUL,KUL,1,A2,Ldh,ipiv,r,N,info);	
	cblas_dcopy(N, r, 1, fi, 1);   // fi = x (i.e. r)
	
	if (info != 0){
		cout <<"ERROR: An error occurred in DGBSV:"<<endl;
		cout<<info<<endl;
	}



}
