#include <iostream>
#include <fstream>
using namespace std;

#include "Poisson.h"
#include "LidDrivenCavity.h"

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
    // Create a new instance of the LidDrivenCavity class
	
    LidDrivenCavity* solver = new LidDrivenCavity();
	Poisson* Psolver = new Poisson();
	int Nx = 161;
	int Ny = 161;
	int NumberofPoints = (Nx-2)*(Ny-2);
	double xlen = 1.0;
	double ylen = 1.0;
	double Re = 200.0;
	double T = Re*Nx*Ny/4.0;
	double dt = 0.0001;
	double *omag = new double [(Ny-2)*(Nx-2)];
	double *fi = new double [(Ny-2)*(Nx-2)];
	//fi = {};
	double *A = new double[(Nx-1)*(Ny-2)*(Nx-2)];
	double *A2 = new double[(3*Nx+1)*(Ny-2)*(Nx-2)];
	//double *VortiInter = new double [(Ny-2)*(Nx-2)];
	//double *streamInter = new double[(Ny-2)*(Nx-2)];
	//double *S_i = new double[(Ny-2)*(Nx-2)];
	//double *S_j = new double[(Ny-2)*(Nx-2)];
	double *V_i = new double[(Ny-2)*(Nx-2)];
	double *V_j = new double[(Ny-2)*(Nx-2)];
	
    // Configure the solver here...
	solver->SetDomainSize(xlen,ylen);
	Psolver->SetDomainSize(xlen,ylen);
	solver->SetGridSize(Nx,Ny,NumberofPoints);
	Psolver->SetGridSize(Nx,Ny);
	solver->SetReynoldsNumber(Re);
	solver->SetTimeStep(dt);
	solver->SetFinalTime(T);
	int counter = 0; //iteration counter
	
	do{		
		
	solver->Initialise(omag,fi,V_i,V_j);
	
	solver->CalVorticityT(A,omag,fi);

	solver->CalVorticityTplus(A,omag,fi,V_i,V_j);

	Psolver->ComputeStreamFunction(fi,omag,A2);
	
	cout<<counter<<endl;
	counter++;
	}while(counter<21);
	
	
//	ofstream stream;
//	stream.open("/home/li/Desktop/hpc/hpc-cw/data.txt");
//	for(int i =0;i<(Ny-2)*(Nx-2);++i){
//		stream<<fi[i];
//		stream<<" ";
//	}
	//stream.close();

	
	delete [] A;
	delete [] omag;
	delete [] fi;
	delete [] A2;
	delete [] V_j;
	delete [] V_i;
	

    // Run the solver
    solver->Integrate();

	return 0;
}