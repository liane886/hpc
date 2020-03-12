#pragma once
#include <iostream>
#include <string>
using namespace std;

class Poisson
{
	
private:
    int    Nx;
    int    Ny;
	double Lx;
    double Ly;
	double* A = nullptr; // PoissonMatrixA
	double beta[3];
	int N;
	double det_y;
	double det_x;
	
public: 
	Poisson();
	~Poisson();
	//void Initialize(int& Nx,int& Ny);
	void buildMA();
	void PPTRF();
	void SolvePoisson(double* s,double*v);
	void SetGridSize (int& Nx,int& Ny);
	void SetDomainSize (double xlen,double ylen);
};