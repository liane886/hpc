#include <iostream>
using namespace std;

#include "LidDrivenCavity.h"
/*
 * • What the program does (name, author, date and description).
	• A brief description of each routine and what the parameters are.
	• A description of what the routine returns.
 */
 
 
int main(int argc, char **argv)
{
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();
	double Nx = 161.0;
	double Ny = 161.0;
	double Re = 200.0;
	double T = Re*Nx*Ny/4.0;
	double dt = 0.00001;
	double *omag = new double [int(Ny)*int(Nx)];
	double *fi = new double [int(Ny)*int(Nx)];
	double *A = new double[160*159*159];
	double *b = new double [int(Nx-2)*int(Ny-2)];
	 //CalV1_b[159*159] = {};
//	double *A = new double [(int(Ny)-2)*(int(Nx)-2)*(int(Ny)-2)*(int(Nx)-2)];
//	double *b = new double [(int(Ny)-2)*(int(Nx)-2)];
    // Configure the solver here...
	solver->SetDomainSize(1.0,1.0);
	solver->SetGridSize(Nx,Ny);
	solver->SetReynoldsNumber(Re);
	solver->SetTimeStep(dt);
	solver->SetFinalTime(T);
    solver->Initialise(omag, fi);
	//cout<<omag[162];
	//cout<<"....";
	//cout<<omag[323]<<endl;
//	cout<<omag[161];
//	cout<<"....";
//	cout<<omag[322];
//	cout<<"....";
//	cout<<omag[25761]<<endl;
	solver->CalVorticityT(A,b);
	cout<<A[160];
	cout<<"....";
	cout<<A[2]<<endl;
//	cout<<A[160];
//	cout<<A[161];
//	cout<<A[162]<<endl;
	delete A;
	//delete b;

    // Run the solver
    solver->Integrate();

	return 0;
}