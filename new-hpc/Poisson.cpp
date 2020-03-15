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
// add  streamfunction BC 
	//vCopy[1] = -beta[2]*s[1];
	//vCopy[Ny] = -beta[1]*s[Ny+2];
	//	cout<<N<<endl;
//	for (int i = 0; i<Nx; i++){
//		for(int j = 0; j< Ny;j++){
//			 vCopy[i] = v[(i+1)*(Ny+2)+1+j]; 	
//		}
//		//F77NAME(dcopy) (Ny, &v[counter*(i+1) + 1], 1, &vCopy[i*Ny], 1);
//		
//	}
//	
	        // Populate vorticity forcing values
        for (int i = 0; i<Nx; i++){
            F77NAME(dcopy) (Ny, &v[counter*(i+1) + 1], 1, &vCopy[i*Ny], 1);
        }

        // Populate streamfunction BC values x-direction
        F77NAME(daxpy) (Ny, -beta[2], &s[1], 1, vCopy, 1);
        F77NAME(daxpy) (Ny, -beta[2], &s[counter*(Nx+1) + 1], 1, &vCopy[Ny*(Nx-1)], 1);

        // Populate streamfunction BC values y-direction
        F77NAME(daxpy) (Nx, -beta[1], &s[Ny +2], (Ny + 2), vCopy, Ny);
        F77NAME(daxpy) (Nx, -beta[1], &s[2*(Ny +2) - 1], (Ny + 2), &vCopy[Ny-1], Ny);

//	for (int i = 0; i<Ny; i++){	
//		vCopy[i] = -beta[2]*s[1];
//	    
//	}	
//	for (int i =Ny*(Nx-1); i<Ny*Nx+1;i++){
//		vCopy[Ny*(Nx-1)+i] = -beta[2]*s[counter*(Nx+1) + 1+i];
//	}
//	
//	for (int i = 0; i<Nx;){
//			i += Ny;
//		vCopy[i] = -beta[1]*s[(Ny+2)+(i+2)];
//		vCopy[Ny-1+i] = -beta[1]*s[2*(Ny+2)-1+(i+2)];
//	}
for (int j = 0;j<Ny;j++){
for (int i = 0;i<Nx;i++){
	cout<<setw(8)<<setprecision(3)<<v[j*Ny+i]<<" ";
}	
cout<<endl;
}

		int info;
// Solve linear system with LAPACK:
     // Packed storage, Banded  definite matrix
		F77NAME(dpbtrs) ('U', N, Ny,1, A,Ny+1,vCopy, N, info);
	// Position solution back in streamfunction array
	
		for(int i =0; i<Nx; i++){
		F77NAME(dcopy) (Ny, &vCopy[i*Ny], 1, &s[counter*(i+1) + 1], 1);
		}
		
		
delete[]vCopy;

//	if (info != 0){
//		cout <<"ERROR: An error occurred :"<<endl;
//		cout<<info<<endl;
//	}


}
