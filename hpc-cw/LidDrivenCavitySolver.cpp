#include <iostream>
using namespace std;

#include "LidDrivenCavity.h"
/*
 * • What the program does (name, author, date and description).
	• A brief description of each routine and what the parameters are.
//	• A description of what the routine returns.
// */
#include "cblas.h"
#define F77NAME(x) x##_
//
//// This tells the C++ compiler to use C-style naming of symbols for these
//// functions, as opposed to the "mangled" symbols used in C++. The need to
//// mangle symbol names is because C++ supports overloading the same function
//// name with different parameters and data types so needs to be able to
//// distinguish them by name in the object code.
extern "C" {
   double F77NAME(ddot) (const int& n, 
                          const double *x, const int& incx,
                         const double *y, const int& incy);
}
double cblas_dnrm2(const int N, const double *X, const int incX);
 
int main(int argc, char **argv)
{
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();
	int Nx = 161;
	int Ny = 161;
	double xlen = 1.0;
	double ylen = 1.0;
	double Re = 200.0;
	double T = Re*Nx*Ny/4.0;
	double dt = 0.0001;
	double *omag = new double [Ny*Nx];
	double *fi = new double [Ny*Nx];
	double *A = new double[(Nx-1)*(Ny-2)*(Nx-2)];
	double *VortiInter = new double [(Ny-2)*(Nx-2)];
	double *streamInter = new double[(Ny-2)*(Nx-2)];
	double *S_i = new double[(Ny-2)*(Nx-2)];
	double *S_j = new double[(Ny-2)*(Nx-2)];
	double *V_i = new double[(Ny-2)*(Nx-2)];
	double *V_j = new double[(Ny-2)*(Nx-2)];
	
    // Configure the solver here...
	solver->SetDomainSize(1.0,1.0);
	solver->SetGridSize(Nx,Ny);
	solver->SetReynoldsNumber(Re);
	solver->SetTimeStep(dt);
	solver->SetFinalTime(T);
    solver->Initialise(omag, fi,S_i,S_j,V_i,V_j);
	//cout<<omag[162];
	//cout<<"....";
	//cout<<omag[323]<<endl;
//	cout<<omag[161];
//	cout<<"....";
//	cout<<omag[322];
//	cout<<"....";
//	cout<<omag[25761]<<endl;
	solver->CalVorticityT(A,VortiInter,streamInter);
	solver->CalVorticityTplus(A,VortiInter,streamInter,S_i,S_j,V_i,V_j);

	cout <<endl;
	delete A;
	//delete b;

    // Run the solver
    solver->Integrate();

	return 0;
}