#include <iostream>
#include <fstream>
using namespace std;

//#include "Poisson.h"
#include "LidDrivenCavity.h"
#include "Poisson.h"

/*
 * • What the program does (name, author, date and description).
	• A brief description of each routine and what the parameters are.
//	• A description of what the routine returns.
// */


//#include "cblas.h"
//#define F77NAME(x) x##_
////
////// This tells the C++ compiler to use C-style naming of symbols for these
////// functions, as opposed to the "mangled" symbols used in C++. The need to
////// mangle symbol names is because C++ supports overloading the same function
////// name with different parameters and data types so needs to be able to
////// distinguish them by name in the object code.
//extern "C" {
//   double F77NAME(ddot) (const int& n, 
//                          const double *x, const int& incx,
//                         const double *y, const int& incy);
//}
//double cblas_dnrm2(const int N, const double *X, const int incX);
// 

int main(int argc, char **argv)
{	
	
	int Nx = 161;
	int Ny = 161;
	int NumberofPoints = (Nx-2)*(Ny-2);
	double xlen = 1.0;
	double ylen = 1.0;
	double Re = 100.0;
	double T = 1;
	double dt = 0.0005;
	//int N = (Nx-2)*(Ny-2);
    // Create a new instance of the LidDrivenCavity class
	
    LidDrivenCavity* solver = new LidDrivenCavity (dt, T,  Nx, Ny, xlen, ylen, Re);
	

	//double *VortiInter = new double [(Ny-2)*(Nx-2)];
	
    // Configure the solver here...
	solver->SetDomainSize(xlen,ylen);
	
	solver->SetGridSize(Nx,Ny,NumberofPoints);
	
	solver->SetReynoldsNumber(Re);
	solver->SetTimeStep(dt);
	solver->SetFinalTime(T);			
		
	solver->Initialise();
 // Run the solver
    solver->Integrate();
	

//////////////////////////////////////////IO
//	int Ldh = 3*Nx+1;
//	for(int i =0;i<2*Nx;++i){
//		for (int j =0;j<2*Nx; ++j){
//			cout<<A[j*Ldh+i]<<" ";
//		}
//			cout<<endl;
//	
//	}

	
//	ofstream stream;
//	stream.open("/home/li/Desktop/hpc/hpc-cw/data.txt");
//	int Ldh = 3*Nx+1;
//	for(int i =0;i<(Nx-1);++i){
//		for (int j =0;j<(Nx-2)*(Ny-2); ++j){
//		stream<<A[j*(Nx-1)+i];
//		stream<<" ";
//		}
//		stream<<"\n";
//	}
//	stream.close();
 
	delete solver;
	return 0;
}