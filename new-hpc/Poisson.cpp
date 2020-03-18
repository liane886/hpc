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
	det_y = Ly/double (Ny-1);
	det_x = Lx/double (Nx-1);
	this->beta[2] = -1.0/(det_y*det_y);
	this->beta[1] = -1.0/(det_x*det_x);
	this->beta[0] = -2.0*(beta[1] + beta[2]);

	A = new double[N*(this->Ny+1)]{};
	
	this ->buildMA();
	this ->CholeskyFactorzation();
}

void Poisson::SetDomainSize(double xlen, double ylen){
	Lx = xlen;
	Ly = ylen;
}

void Poisson::buildMA(){

	//generate matrix A;
	int Ud = Ny +1;
	
	A[Ny] = beta[0];
	for (int i = 1; i < N; ++i) {		
		A[Ud*(i+1)-1] = beta[0]; // diagonal
		
		if (i >= Ny-1) {
			A[(i+1)*Ud] = beta[1];  //super-diagonal
		}	

		if (i % Ny != 0){
			A[(i+1)*Ud -2] = beta[2];  //sub-diagonal
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

void Poisson::CholeskyFactorzation(){
	int info;
	//conpute Cholesky factoration
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

	// set the vorticity vector for interial vorticity (length:(Nx-2)*(Ny-2))
//        for (int i = 1; i<Nx+1; i++){
//			for (int j=1;j<Ny+1;j++){
//				int counter1 = 0;
//				vCopy[counter1] = v[i*(Ny+2)+j];
//				counter1+=1; 
//			}
//        }

      for (int i = 0; i<Nx; i++){
            F77NAME(dcopy) (Ny, &v[counter*(i+1) + 1], 1, &vCopy[i*Ny], 1);
        }
		
//        // Populate streamfunction BC values x-direction
//        F77NAME(daxpy) (Ny, -beta[1], &s[1], 1, vCopy, 1);
//        F77NAME(daxpy) (Ny, -beta[1], &s[counter*(Nx+1) + 1], 1, &vCopy[Ny*(Nx-1)], 1);
//
//        // Populate streamfunction BC values y-direction
//        F77NAME(daxpy) (Nx, -beta[2], &s[Ny +2], (Ny + 2), vCopy, Ny);
//        F77NAME(daxpy) (Nx, -beta[2], &s[2*(Ny +2) - 1], (Ny + 2), &vCopy[Ny-1], Ny);


		int info;
// Solve linear system with LAPACK:
     // solve the linear equaiton A*x=B with a symmertric positive define band matrix A 
	 //using Cholesky factorization A computed by DPBTRF.
		F77NAME(dpbtrs) ('U', N, Ny,1, A,Ny+1,vCopy, N, info);
		
	// Position solution back in streamfunction array	
      for (int i = 0; i<Nx; i++){
            F77NAME(dcopy) (Ny,&vCopy[i*Ny] , 1, &s[counter*(i+1) + 1], 1);
        }
		
delete[]vCopy;

	if (info != 0){
		cout <<"ERROR: An error occurred in DPBTRS:"<<endl;
		cout<<info<<endl;
	}


}
