#include <cstdlib>
#include <cmath>
#include "Poisson.h"
#include <fstream>
#include <iomanip>
#include "/home/li/Desktop/header/cblas/CBLAS/include/cblas.h"


// This tells the C++ compiler to use C-style naming of symbols for these
// functions, as opposed to the "mangled" symbols used in C++. The need to
// mangle symbol names is because C++ supports overloading the same function
// name with different parameters and data types so needs to be able to
// distinguish them by name in the object code.
#define F77NAME(x) x##_
extern "C" { 
                                                             //Chos
	double F77NAME(dpptrf) (const char& UPLO, const int& n,
                            const double* AP, int& info);


//mutiply 
    double F77NAME(daxpy) (const int& n, const double& alpha,
                          const double *x, const int& incx,
                          const double *y, const int& incy);
	
// copy the array x to y					  
	double F77NAME(dcopy) (const int& n,
                          const double *x, const int& incx, 
                          const double *y, const int& incy);
						  
	double F77NAME(dpptrs) (const char& UPLO, const int& n,
                            const int& nRHS, const double* AP, 
                            const double* b, const int& ldb,
                            int& info);
}

Poisson::Poisson(){

}
Poisson::~Poisson(){
	cout<<"object is being deleted"<<endl;
}

void Poisson::SetGridSize(int& Nx, int& Ny)
{
	this->Nx = Nx-2;
	this->Ny = Ny-2;
	this->N = (this->Nx)*(this->Ny);
	det_y = Ly/double (Ny-1);
	det_x = Lx/double (Nx-1);
	beta[0] = 2*(det_x*det_x+det_y*det_y);
	beta[1] = -det_y*det_y;
	beta[2] = -det_x*det_x;
	A = new double[N*(N+1)/2];
	this ->buildMA();
	cout<<"123"<<endl;
	this ->PPTRF();
	cout<<"1455523"<<endl;
}

void Poisson::SetDomainSize(double xlen, double ylen){
	Lx = xlen;
	Ly = ylen;
}

void Poisson::buildMA(){
	
	//generate matrix A;
	A[0] = 	2*(det_x*det_x+det_y*det_y);	
	for (int i = 1; i < N; ++i) {
		int counter = (i+1)*(i+2)/2 -1;
		
		A[counter] = beta[0]; // diagonal
		if (i >= Ny) {
			A[counter - Ny] = beta[1];  //super-diagonal
		}	

		if (i % Ny != 0){
			A[counter - 1] = beta[2];  //sub-diagonal
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
	F77NAME(dpptrf)('U', N, A, info);
	if (info != 0){
		cout <<"ERROR: An error occurred in dpptrf:"<<endl;
		cout<<info<<endl;
	}
}

void Poisson::SolvePoisson(double* s,double*v){

	
	double* vCopy = new double[N];  // update vorticity vector	     
	fill_n(vCopy,N,0.0);	
// add  streamfunction BC 

	vCopy[1] = -beta[2]*s[1];
	vCopy[Ny] = -beta[1]*s[Ny+2];
	int counter = Ny + 2;	
	for (int i = 0; i<Nx; i++){
		
		 vCopy[i*Ny] = v[counter*(i+1) + 1]; 	
 //F77NAME(dcopy) (Ny, &v[counter*(i+1) + 1], 1, &vCopy[i*Ny], 1);	 
	}

        F77NAME(daxpy) (Ny, -beta[2], &s[counter*(Nx+1) + 1], 1, &vCopy[Ny*(Nx-1)], 1);


        F77NAME(daxpy) (Nx, -beta[0], &s[2*(Ny +2) - 1], (Ny + 2), &vCopy[Ny-1], Ny);
		int info;
//        // Solve linear system with LAPACK:
//        // Packed storage, Symetric Positive definite matrix
        F77NAME(dpptrs) ('U', N, 1, A, vCopy, N, info);
//
//        // Position solution back in streamfunction array
        //int offset = Ny + 2;
        for(int i =0; i<Nx; i++){
            F77NAME(dcopy) (Ny, &vCopy[i*Ny], 1, &s[counter*(i+1) + 1], 1);
        }
//	
//	
//	
//	int KUL = Nx-2;
//	int* ipiv = new int[N];
//	int info;
//	double* r = new double[N];
//
//	cblas_dcopy(N, v, 1, r, 1);   // r_0 = b (i.e. omag)
//	
//
//	F77NAME(dgbsv)(N,KUL,KUL,1,A2,Ldh,ipiv,r,N,info);	
//	cblas_dcopy(N, r, 1, s, 1);   // fi = x (i.e. r)
//	
//	if (info != 0){
//		cout <<"ERROR: An error occurred in DGBSV:"<<endl;
//		cout<<info<<endl;
//	}


}
