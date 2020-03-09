//#pragma once
#include <iostream>
#include <string>
using namespace std;

/** : Using Doxygen to document a class.
* @class LidDrivenCavity
 * @brief Describes a position on the 2D plane.
 */

class LidDrivenCavity
{
public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void Initialise(double* omag,double* fi,double*V_i,double* V_j);
	void CalVorticityT(double* A,double* VortiInter,double* streamInter);
	void CalVorticityTplus(double* A, double* VortiInter,double* streamInter,double*V_i,double* V_j);
	
    void Integrate();

    // Add any other public functions

private:
    double* v = nullptr;
    double* s = nullptr;
    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
	double* CalVI_A;
	double* CalVI_y;
	double* CalVI_x;
	const int N = (Nx-2)*(Ny-2);  //the size of matrix A and array of inter vorticity & stream  
	double *LhsVortiAns_i = nullptr;
	
	double *LhsVortiAns_j = nullptr;
	
	
//	double* vorticity = new double [Ny*Nx];
//	double* stream = new double [Ny*Nx];
//	double* CalV1_A = new double [(Ny-2)*(Nx-2)*(Ny-2)*(Nx-2)];
//	double* CalV1_b = new double [(Ny-2)*(Nx-2)];
};
