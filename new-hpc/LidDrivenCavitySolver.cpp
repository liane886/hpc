#include <iostream>
#include <fstream>
using namespace std;

#include "LidDrivenCavity.h"
#include "Poisson.h"
//#include <mpi.h>
#include <boost/program_options.hpp>
// An alias to reduce typing
namespace po = boost::program_options;

/*
 * •Parallel code to solve lid-cavity problem; Author: Liping Wang;  Date: 23/3/2020
	This prpgram is using to solve the lid-cavity problem, 
	 * enter the variables dirctly in the terminal to run the program. 
// */

//void splitDomian4MPI(int grids,int partitionNum){                   <------------------------  for mpi
//	
////Split the grid points according to the number of partitions
//	int subX;subX_Remaind;subY;subY_Remaind;
//	//in x-direction
//	subX = grids[0] / partitionNum[0];
//	subX_Remaind = grids[0] % partitionNum[0];
//	//in y -direction
//	subY = grids[1] / partitionNum[1];
//	subY_Remaind = grids[1] % partitionNum[1];
//
//	
//	if (coords[0] <  subX_Remaind ){
//		subGrids [0] = subX + 1;
//	}
//	else{
//		subGrids [0] = subX;
//	}
//
//	if (coords[1] <  subY_Remaind ){
//		subGrids [1] = subY + 1;
//	}
//	else{
//		subGrids [1] = subY;
//	}
//	
//}


int main(int argc, char **argv)
{	

// set the program option using boost library
	po::options_description opts(
		"Enter the variable for solving lid-cavity problem");
    opts.add_options()
        ("help" , "Produce help message.")
        ("Lx" , po::value<double>() -> default_value(1.0)    , "Length of the domain in the x-direction.")
        ("Ly" , po::value<double>() -> default_value(1.0)    , "Length of the domain in the y-direction.")
        ("Nx" , po::value<int>() -> default_value(161)       , "Number of grid points in x-direction.")
        ("Ny" , po::value<int>() -> default_value(161)       , "Number of grid points in y-direction.")
        ("Px" , po::value<int>() -> default_value(1)         , "Number of partitions in x-direction.")
        ("Py" , po::value<int>() -> default_value(1)         , "Number of partitions in y-direction.")
        ("dt" , po::value<double>() -> default_value(0.01)   , "Time step size.")
        ("T"  , po::value<double>() -> default_value(10.0)   , "Final time.")
        ("Re" , po::value<double>() -> default_value(1000.0) , "Reynolds number.");
		
    // Tell Boost to parse the command-line arguments using the list of
    // possible options and generate a map (vm) containing the options and
    // values actually specified by the user.
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);
	
	//Save the inputs in variables 
	double Lx = vm["Lx"].as<double>();
	double Ly = vm["Ly"].as<double>();
	int Nx = vm["Nx"].as<int>();
	int Ny = vm["Ny"].as<int>();
	int partitionNum[2];
	partitionNum[0] = vm["Py"].as<int>();
	partitionNum[1] = vm["Px"].as<int>();
	double dt = vm["dt"].as<double>();
	double T = vm["T"].as<double>();
	double Re = vm["Re"].as<double>();
	int NumberofPoints = (Nx-2)*(Ny-2);

	double dx = Lx/Nx;
	double dy = Ly/Ny;
	double coef = (Re* dx*dy)/4;
	
	if (dt>= coef){
		cout<<"Please reselecte dt. It shouble be less than (Re*dx*dy)/4"<<endl;
		return 0;
	}
	
	//   -------------------------------------------------------------Set MPI----------------------------------------------//
//	int grids[2];
//	splitDomian4MPI(grids,partitionNum);
//	int np;
//	np = partitionNum[0]*partitionNum[1];
//
//
//    // Initialise MPI & test MPI return values
//	int retval;
//    retval = MPI_Init(&argc, &argv);
//	if (retval != MPI_SUCCESS){
//		cout<<"An error occurred initialising MPI"<<endl;
//	}
//    // get the rank of the processes and the number of the processes  
//    int rank, size;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//	
//	if(size!= np){
//		cout<<"MPI size error"<<endl;
//		break;
//	}
//
//// SPLITE DOMIAN TO SUBDOMAIN
//
//    // Initialize cartesian grid communicator
//    MPI_Comm subdomain;
//    int periods[2] = {0, 0}, coords[2];
//    int reorder = 0;
//    MPI_Cart_create(MPI_COMM_WORLD, 2, partitionDims, periods, reorder, &subdomain);
//    MPI_Cart_coords(subdomain, rank, 2, coords);
//
//    int ShiftRank[4] = {rank,rank,rank,rank}; // bottom;top;left;right;
//    MPI_Cart_shift(subdomain, 0, 1, &ShiftRank[0], &ShiftRank[1]);
//    MPI_Cart_shift(subdomain, 1, 1, &ShiftRank[2], &ShiftRank[3]);
// -----------------------------------------------------------------------------------------------------------------------------------//

    // Create a new instance of the LidDrivenCavity class
   // LidDrivenCavity* solver = new LidDrivenCavity (MPI_COMM_WORLD,ShiftRank, dt, T,  Nx, Ny, Lx, Ly, Re);  <---used for mpi
	LidDrivenCavity* solver = new LidDrivenCavity ( dt, T, Nx, Ny, Lx, Ly, Re);
	
	
    // Configure the solver here...
	solver->SetDomainSize(Lx,Ly);
	solver->SetGridSize(Nx,Ny,NumberofPoints);
	solver->SetReynoldsNumber(Re);
	solver->SetTimeStep(dt);
	solver->SetFinalTime(T);			
	solver->Initialise();
	// Run the solver
    solver->Integrate();
	
// output the result 
//	string filepath_stream;
//	string filepath_vorticity,filepath_vh,filepath_vv;
//	cout<<"PLease enter the file path and name to get the result of stream:"<<endl;
//	getline(cin,filepath_stream);
//	cout<<"PLease enter the file path and name to get the result of vorticity:"<<endl;
//	getline(cin,filepath_vorticity);
//	cout<<"PLease enter the file path and name to get the result of horizontal velocity:"<<endl;
//	getline(cin,filepath_vh);
//	cout<<"PLease enter the file path and name to get the result of vertical velocity:"<<endl;
//	getline(cin,filepath_vv);

	solver->PrintResult2file();
	
	
	
	delete solver; //clear the memory
  //-------------------------------------------------------------  MPI_Finalize();
	return 0;
}