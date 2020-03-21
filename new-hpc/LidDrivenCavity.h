//#pragma once
#include <iostream>
#include <string>
#include "Poisson.h"
using namespace std;

/** : Using Doxygen to document a class.
* @class LidDrivenCavity
 * @brief Describes a position on the 2D plane.
 */

class LidDrivenCavity
{
public:
    LidDrivenCavity(double dt,double T, int Nx,int Ny,double Lx, double Ly,double Re);
    ~LidDrivenCavity();
	
    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny,int n);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void Initialise();
	void CalVorticityT(double alpha,double* vorticity_inter,double* s);
	void CalVorticityTplus( double* v,double* s,double* vorticity_inter);
	void BoundaryCondition();
    void Integrate();
	void PrintResult2file();
	void Velocity(double* s,double* velocity_dx,double* velocity_dy);
    void splitDomian4MPI();
	void MPIrecv(double* recv)
	void MPIsend(double* send);
private:
	Poisson* Psolver;
    double* v = nullptr;
    double* s = nullptr;
    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
	int N;  //the size of matrix A and array of inter vorticity & stream  
	double* vorticity_inter = nullptr;
	double* velocity_dx = nullptr;
	double* velocity_dy = nullptr;
};
