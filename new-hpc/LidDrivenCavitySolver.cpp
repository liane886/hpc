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



int main(int argc, char **argv)
{	
	
	int Nx = 161;
	int Ny = 161;
	int NumberofPoints = (Nx-2)*(Ny-2);
	double xlen = 1.0;
	double ylen = 1.0;
	double Re = 200.0;
	double T = 1;
	double dt = 0.0005;


    // Initialise MPI
    MPI_Init(&argc, &argv);

    // MPI related variables
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

// SPLITE DOMIAN TO SUBDOMAIN

    // Initialize cartesian grid communicator
    MPI_Comm cartGrid;
    int periods[2] = {0, 0}, coords[2];
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, partitionSize, periods, reorder, &cartGrid);
    MPI_Cart_coords(cartGrid, rank, 2, coords);

    int rankShift[4] = {rank,rank,rank,rank};
    MPI_Cart_shift(cartGrid, 0, 1, &rankShift[0], &rankShift[1]);
    MPI_Cart_shift(cartGrid, 1, 1, &rankShift[2], &rankShift[3]);




    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity (dt, T,  Nx, Ny, xlen, ylen, Re);
	

	
    // Configure the solver here...
	solver->SetDomainSize(xlen,ylen);
	
	solver->SetGridSize(Nx,Ny,NumberofPoints);
	
	solver->SetReynoldsNumber(Re);
	solver->SetTimeStep(dt);
	solver->SetFinalTime(T);			
		
	solver->Initialise();
 // Run the solver
    solver->Integrate();
	




 
	delete solver;
	return 0;
}