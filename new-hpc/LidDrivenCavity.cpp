#include "LidDrivenCavity.h"
#include "Poisson.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>

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

	 //  row major with [+ve y-direction -> down , +ve x-direction -> right]
	 
//BC for vorticity
//	for (int i = 0; i < this -> Nx ; i++){ 
//		v[i*Ny] = (s[i*Ny]-s[i*Ny +1 ])*(2.0/(det_y*det_y)) - (2.0*U/det_y);  	//top BC
//		v[(i+1)*Ny - 1] = (s[(i+1)*Ny - 1]-s[(i+1)*Ny - 2])*(2.0/(det_y*det_y));  //Bottom BC
//	} 
//
//	for (int i = 0; i < this -> Ny; i++){
//		v[i] = (s[i]-s[i+Ny])*(2.0/(det_x*det_x));  //left BC
//		v[i + (Nx - 1)*Ny] = (s[i + (Nx - 1)*Ny]-s[i + (Nx - 2)*Ny])*(2.0/(det_x*det_x)); //right BC
//	}	
//	
       
//BC for streamfunction

	for (int i =0;i<Nx;i++){
		s[i*Ny] = 0;
		s[(i+1)*Ny-1] = 0;
		
	}

	for (int j=0;j<Ny;j++){
		s[j] =0;
		s[j+(Nx-1)*Ny] = 0;
	}
//Bottom BC
	for (int i = 0; i < this -> Nx; i++){
		v[i*Ny] = (s[i*(Ny)]-s[i*(Ny)+1])*(2.0/(det_y*det_y)); 
	}
//top BC			
	for (int i = 0; i < this -> Nx ; i++){ 
		v[(i+1)*Ny - 1] = (s[(i+1)*Ny - 1]-s[(i+1)*Ny - 2])*(2.0/(det_y*det_y)) - (2.0*U/det_y);  	
	}
			
			
//left BC
	for (int i = 0; i < this -> Ny; i++){
		v[i] = (s[i]-s[i+Ny])*(2.0/(det_x*det_x));
	}
			
//right BC
       
	for (int i = 0; i < this -> Ny; i++){
		v[i + (Nx - 1)*Ny] = (s[i + (Nx - 1)*Ny]-s[i + (Nx - 2)*Ny])*(2.0/(det_x*det_x));
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
		
	
	CalVorticityT(-1.0,this->v,this->s);
	//cout<<v[180]<<endl;	//Initialise();	
	
	// solve interior vorticity at t+dt
	F77NAME(dcopy)(Ny*Nx, v, 1, vorticity_inter, 1); 
			
	CalVorticityT(dt/Re,vorticity_inter,v);

	CalVorticityTplus(v,s,vorticity_inter);
	        // Add the value of v to v_new
     for (int i=1; i < Nx-1; i++ ){
          for (int j=1; j< Ny -1; j++){
			vorticity_inter[ j + Ny*i ] += v[ j + Ny*i ];
		}
	}
		
	F77NAME(dcopy)(Ny*Nx, vorticity_inter, 1, v, 1);
	
	
	//for (int k = 0; k < 6; k++ ){
	Psolver -> SolvePoisson(this->s,this->v);
	//}
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
							beta_j*(s[(j+1) + Ny*i] - s[(j-1) + Ny*i])*beta_i*(v[j + Ny*(i+1)] - v[j + Ny*(i-1)]));
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


void LidDrivenCavity::PrintResult2file(){
	
	ofstream stream;
	ofstream vorticity;
	ofstream V_h;
	ofstream V_v;
	stream.open("/home/li/Desktop/hpc/new-hpc/s.txt");
	vorticity.open("/home/li/Desktop/hpc/new-hpc/v.txt");
	V_h.open("/home/li/Desktop/hpc/new-hpc/horizontalV.txt");
	V_v.open("/home/li/Desktop/hpc/new-hpc/verticalV.txt");
	//stream<<"998"<<endl;
	for(int i =0;i<Nx*Ny;++i){
			
			stream<<s[i]<<endl;
			vorticity<<v[i]<<endl;
			V_h << velocity_dx[i]<<endl;
			V_v << velocity_dy[i]<<endl;
		
	}
	stream.close();
	vorticity.close();
	V_h.close();
	V_v.close();
}



// void LidDrivenCavity::InterfaceBroadcast(double* field){
//        
//        //Sequentially send interface values to neighbors in all directions
//        MPI_Barrier(MPIcomm);
//        // Neighbor below
//        if (rankShift[0] != -2){
//            LidDrivenCavity::InterfaceSend(Nx, &field[1], bufNx, Ny, rankShift[0],rank,MPIcomm);
//        }
//
//        // Neighbor above
//        if (rankShift[1] != -2){
//            LidDrivenCavity::InterfaceSend(Nx, &field[Ny-2], bufNx, Ny, rankShift[1],rank,MPIcomm);
//        }
//
//        // Neighbor leftward
//        if (rankShift[2] != -2){
//            LidDrivenCavity::InterfaceSend(Ny, &field[Ny], bufNy, 1, rankShift[2],rank,MPIcomm);
//        }
//
//        // Neighbor rightward
//        if (rankShift[3] != -2){
//            LidDrivenCavity::InterfaceSend(Ny, &field[Ny*(Nx-2)], bufNy, 1, rankShift[3],rank,MPIcomm);
//        }
//
//    }
//
//    void LidDrivenCavity::InterfaceGather(double* field){
//        
//        //Sequentially send interface values to neighbors in all directions
//
//        //Neighbor above
//        if (rankShift[1] != -2){
//            LidDrivenCavity::InterfaceRecv(Nx, &field[Ny-1], bufNx, Ny, rankShift[1], rankShift[1], MPIcomm);
//        }
//
//        //Neighbor below
//        if (rankShift[0] != -2){
//            LidDrivenCavity::InterfaceRecv(Nx, field, bufNx, Ny, rankShift[0], rankShift[0], MPIcomm);
//        }
//
//        //Neighbor rightward
//        if (rankShift[3] != -2){
//            LidDrivenCavity::InterfaceRecv(Ny, &field[Ny*(Nx-1)], bufNy, 1, rankShift[3], rankShift[3], MPIcomm);
//        }
//
//        //Neighbor leftward
//        if (rankShift[2] != -2){
//            LidDrivenCavity::InterfaceRecv(Ny, field, bufNy, 1, rankShift[2], rankShift[2], MPIcomm);
//        }
//        MPI_Barrier(MPIcomm);
//    }
//
//
//    void LidDrivenCavity::InterfaceSend(int& count, double* field, double* buff, int disp, int& dest, int& tag, MPI_Comm MPIcomm){
//        F77NAME(dcopy) (count, field, disp, buff, 1);
//        MPI_Send(buff, count, MPI_DOUBLE, dest, tag, MPIcomm);
//    }
//
//    void LidDrivenCavity::InterfaceRecv(int& count, double* field, double* buff, int disp, int& source, int& tag, MPI_Comm MPIcomm){
//        MPI_Recv(buff, count, MPI_DOUBLE, source, tag, MPIcomm, MPI_STATUS_IGNORE);
//        F77NAME(dcopy) (count, buff, 1, field, disp);
//    }