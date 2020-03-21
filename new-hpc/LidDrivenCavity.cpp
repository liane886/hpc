#include "LidDrivenCavity.h"
#include "Poisson.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <time.h>
//#include "cblas.h"
//#include "/home/li/Desktop/header/cblas/CBLAS/include/cblas.h"



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
}
LidDrivenCavity::LidDrivenCavity(double dt,double T, int Nx,int Ny,double Lx, double Ly,double Re)
{	 
     this ->T = T;
     this ->Nx = Nx;
	 this ->Ny = Ny;
     this ->Lx = Lx;
     this ->Ly = Ly;
     this ->Re = Re;
	 this ->dt = dt;
	 //this ->N = N;
	
	// this ->vorticity_inter = {};
	 cout << "Object is being created" << endl;
}

// test

LidDrivenCavity::~LidDrivenCavity()
{
	cout << "Object is being deleted" << endl;	// print information
//	delete[]v;
//	delete[]vorticity_inter;
//	delete[]s;
//	delete[]velocity_dx;
//	delete[]velocity_dy;
	
}


/**  Using Doxygen to document a function
 * @brief Perform a single insert step in the insertion sort.
 * @param xlen Length of hte domain in the x-dirction
 * @param ylen Length of hte domain in the y-dirction
 *
 */
void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
	this ->Lx = xlen;
	this ->Ly = ylen;

}

void LidDrivenCavity::SetGridSize(int nx, int ny,int NumberofPoints)
{
	this ->Nx = nx;
	this ->Ny = ny;
	this ->N = NumberofPoints;
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
	this ->dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
	this ->T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
	this ->Re = re;
}


void LidDrivenCavity::Initialise()
{
	//count the size of subdomain
	if (ShiftRank[0] != MPI_PROC_NULL){
		Ny += 1;
	}
	if (ShiftRank[1] != MPI_PROC_NULL){
		Ny += 1;
	}
	if (ShiftRank[2] != MPI_PROC_NULL){
		Nx += 1;
	}
	if (ShiftRank[3] != MPI_PROC_NULL){
		Nx += 1;
	}

	
	this ->v = new double[Nx*Ny]{};
	this ->s = new double[Nx*Ny]{};
	this ->vorticity_inter = new double[Nx*Ny]{};
	this ->	velocity_dx =new double [Nx*Ny]{};
	this -> velocity_dy =new double [Nx*Ny]{};

	//Initialise boundary condition
	this->BoundaryCondition();
        

	Psolver = new Poisson();	
	Psolver -> SetDomainSize(Lx,Ly);
	
	Psolver -> PInitialise(Nx,Ny);

}

void LidDrivenCavity:: BoundaryCondition(){
	double U = 1.0;
	double det_y = Ly/double (Ny-1);
	double det_x = Lx/double (Nx-1);

//  column major with [+ve right & down]

//BC for streamfunction & vorticity
	//Bottom BC
	if (ShiftRank[0] == MPI_PROC_NULL){
		for (int i = 0; i < this -> Nx; i++){
			
			v[i*Ny] = (s[i*(Ny)]-s[i*(Ny)+1])*(2.0/(det_y*det_y)); 
			s[i*Ny] = 0;
		
	}
	
	//top BC
	if (ShiftRank[1] == MPI_PROC_NULL){
		
		for (int i = 0; i < this -> Nx ; i++){ 
			v[(i+1)*Ny - 1] = (s[(i+1)*Ny - 1]-s[(i+1)*Ny - 2])*(2.0/(det_y*det_y)) - (2.0*U/det_y);  
			s[(i+1)*Ny-1] = 0;	
	}
	//left BC		
	if (ShiftRank[2] == MPI_PROC_NULL){
		for (int i = 0; i < this -> Ny; i++){
			v[i] = (s[i]-s[i+Ny])*(2.0/(det_x*det_x));
			s[i] =0;
		}
	}	
	//right BC	
	if (ShiftRank[2] == MPI_PROC_NULL){
		for (int i = 0; i < this -> Ny; i++){
			s[j+(Nx-1)*Ny] = 0;
			v[i + (Nx - 1)*Ny] = (s[i + (Nx - 1)*Ny]-s[i + (Nx - 2)*Ny])*(2.0/(det_x*det_x)); 
		}
	}
			

       

        
}

/**  Solve inter vortisity 
 * @brief Solve inter vortisity by y = A*x --> y    && using cblas-dsbmv (band& symmetric matrix)
 * @param 		A 			Matrix  A in the equation y = A*x 
 * @param 	VortiInter 		inter vorticity (y in the equation y = A*x)
 *
 */
 
 void LidDrivenCavity::Integrate()
{
	double t = 0.0;
	int counter =0;

	do{
//clock_t start,finish;
		

	CalVorticityT(-1.0,this->v,this->s);
	
	// solve interior vorticity at t+dt
	F77NAME(dcopy)(Ny*Nx, v, 1, vorticity_inter, 1); 
			
	CalVorticityT(dt/Re,vorticity_inter,v);

	CalVorticityTplus(v,s,vorticity_inter);
	
//	 // Add the value of v to v_new
//     for (int i=1; i < Nx-1; i++ ){
//          for (int j=1; j< Ny -1; j++){
//			vorticity_inter[ j + Ny*i ] += v[ j + Ny*i ];
//		}
//	}
		
	F77NAME(dcopy)(Ny*Nx, vorticity_inter, 1, v, 1);
//
//double  totaltime;
//start = clock();
	//for (int k = 0; k < 6; k++ ){
	Psolver -> SolvePoisson(this->s,this->v);
//finish = clock();
//totaltime=(double)(finish - start)/CLOCKS_PER_SEC;
//cout<<"TIME takes using 'dpbtrs' for each iteration     "<<totaltime<<endl;

	counter +=1;
	cout<<counter<<endl;
	
	//update bc
	BoundaryCondition();

	t += dt;
	}while(t<T);
	
	// calculate horizontal and vertical velocity 
	Velocity(s,velocity_dx,velocity_dy);
	PrintResult2file();
	
}


void LidDrivenCavity::CalVorticityT(double alpha, double* y,double* x)
{
		
//generate matrix A 
	double det_y = Ly/double (Ny-1);
	double det_x = Lx/double (Nx-1);

	for (int i=1; i < Nx-1; i++){
		for (int j=1; j< Ny-1; j++){
			y[ j + Ny*i ] = alpha*(( x[(j+1) + Ny*i] - 2.0*x[j + Ny*i] + x[(j-1) + Ny*i] )/(det_y*det_y) +
                                ( x[j + Ny*(i+1)] - 2.0*x[j + Ny*i] + x[j + Ny*(i-1)] )/(det_x*det_x));
		}
	}
	
}

/**  Solve the equation 11
 * @brief Solve inter vortisity by y = A*x --> y    && using cblas-dsbmv (band& symmetric matrix)
 * @param 		A 			Pointer to vector of length N*(Nx-1). Matrix A in the equation y = A*x is stored in this way.
 * @param 	VortiInter 		Pointer to vector of length N & inter vorticity (y in the equation y = A*x)
 *
 */
 
void LidDrivenCavity::CalVorticityTplus(double* v,double* s,double*vorticity_inter)
{	
	double det_y = Ly/double (Ny-1);
	double det_x = Lx/double (Nx-1);
	double beta_i = 0.5/det_x;
	double beta_j = 0.5/det_y;

	for (int i=1; i < Nx-1; i++ ){
            for (int j=1; j< Ny-1; j++){
			vorticity_inter [j + Ny*i] += dt*(beta_i*(s[j + Ny*(i+1)]-s[j + Ny*(i-1)])*beta_j*(v[(j+1) + Ny*i] - v[(j-1) + Ny*i]) -
							beta_j*(s[(j+1) + Ny*i] - s[(j-1) + Ny*i])*beta_i*(v[j + Ny*(i+1)] - v[j + Ny*(i-1)])) + v[ j + Ny*i ];
            }
	}

}


void LidDrivenCavity::Velocity(double* str, double* velocoty_dx, double* velocoty_dy){
	
	double det_y = Ly/double (Ny-1);
	double det_x = Lx/double (Nx-1);

	for (int i=1; i < Nx-1; i++ ){
		for (int j=1; j< Ny-1; j++){
			velocity_dx[j + Ny*i] = (str[Ny*(i+1)+j] - str[ Ny*(i-1) +j])/(2.0*det_x);
			velocity_dy[j + Ny*i] = (str[ Ny*i + (j+1)] - str[ Ny*i + (j-1)])/(2.0*det_y);
		}
	}
}


void LidDrivenCavity::PrintResult2file(filepath_stream,filepath_vorticity,filepath_vh,filepath_vv){
	
// Get the number of processes
	int size;
	MPI_Comm_size(MPIcomm, &size);

	// Counte to get the interface point
	int yStart = 0, yEnd = 0, xStart = 0, xEnd = 0;

	if (ShiftRank[0] != -2) yStart++;
	if (ShiftRank[1] != -2) yEnd++;
	if (ShiftRank[2] != -2) xStart++;
	if (ShiftRank[3] != -2) xEnd++;
		
// Open file and set as output only and overwrite
	ofstream stream(filepath_stream, ios::out | ios::trunc);
	ofstream vorticity(filepath_vorticity, ios::out | ios::trunc);
	ofstream V_h(filepath_vh, ios::out | ios::trunc);  //file stored horizontal velocity
	ofstream V_v(filepath_vv, ios::out | ios::trunc);  //file stored vertical velocity

	for (int k = 0; k < size; k++){
		if (k == 0 && k ==rank  ){

			for (int i = xStart; i < (Nx - xEnd); i++ ){
				for (int j = yStart; j < (Ny - yEnd); j++ ){
					stream << s[j + i*Ny] << endl;
					vorticity << v[j + i*Ny] << endl;
					V_h << velU[j + i*Ny] << endl;
					V_v << velV[j + i*Ny] << endl;
				}
			}
		}
		else if (k == rank){
				
			for (int i = xShift_Start; i < (Nx - xShift_End); i++ ){
				for (int j = yShift_Start; j < (Ny - yShift_End); j++ ){
					stream << s[j + i*Ny] << endl;
					vorticity << v[j + i*Ny] << endl;
					V_h << velU[j + i*Ny] << endl;
					V_v << velV[j + i*Ny] << endl;
				}
			}

		}
			
            // Ensure the data is writed to the file by a single process at a time
		MPI_Barrier(MPIcomm);
	}
	stream.close();
	vorticity.close();				
	V_h.close();
	V_v.close();
}

          

void LidDrivenCavity::MPIsend(double* send){
        
		bufNx   = new double [Nx]{};   
		bufNy = new double [Ny]{};
		
//send values to neighbors in all directions
        
        // Neighbor below
        if (ShiftRank[0] != -2){
			
			F77NAME(dcopy) (Nx, &send[1], Ny, bufNy, 1);
			MPI_Send(bufNy, Nx, MPI_DOUBLE, ShiftRank[0], rank, MPIcomm)
           
        }   

        // Nbor above
        if (ShiftRank[1] != -2){
			
			F77NAME(dcopy) (Nx, &send[Ny-2], Ny, bufNx, 1);
			MPI_Send(bufNx, Nx, MPI_DOUBLE, ShiftRank[1], rank, MPIcomm)
			
        }

        // Nbor leftward
        if (ShiftRank[2] != -2){
		
			F77NAME(dcopy) (Ny, &send[Ny], Ny, bufNy, 1);
			MPI_Send(bufNy, Nx, MPI_DOUBLE, ShiftRank[2], rank, MPIcomm)
			
            
        }

        // Nbor rightward
        if (ShiftRank[3] != -2){
			
			F77NAME(dcopy) (Ny, &send[Ny*(Ny-2)], Ny, bufNy, 1);
			MPI_Send(bufNy, Ny, MPI_DOUBLE, ShiftRank[3], rank, MPIcomm)
			
           
        }

}



void LidDrivenCavity::MPIrecv(double* recv){
        
		bufNx   = new double [Nx]{};   
		bufNy = new double [Ny]{};
    
 //recv the values from neighbor
 
        //Neighbor below
        if (ShiftRank[0] != -2){
		MPI_Recv(bufNx, count, MPI_DOUBLE, ShiftRank[0], ShiftRank[0], MPIcomm, MPI_STATUS_IGNORE);
        F77NAME(dcopy) (Nx, bufNx, 1, recv, Ny);
        }
		
        //Neighbor above
        if (ShiftRank[1] != -2){
			
			
		MPI_Recv(bufNx, Nx, MPI_DOUBLE, ShiftRank[1], ShiftRank[1], MPIcomm, MPI_STATUS_IGNORE);
        F77NAME(dcopy) (Nx, bufNx, 1, &recv[Ny-1], Ny);
			
		
        }
		
       //Neighbor leftward
        if (ShiftRank[2] != -2){
		MPI_Recv(bufNy, count, MPI_DOUBLE, ShiftRank[2], ShiftRank[2], MPIcomm, MPI_STATUS_IGNORE);
        F77NAME(dcopy) (Ny, bufNy, 1, recv, Nx);
		
        }
		

        //Neighbor rightward
        if (ShiftRank[3] != -2){
		MPI_Recv(bufNy, count, MPI_DOUBLE, ShiftRank[3], ShiftRank[3], MPIcomm, MPI_STATUS_IGNORE);
        F77NAME(dcopy) (Ny, bufNy, 1, &recv[Ny*(Nx-1)], Nx);
		
        }

 
        
}


