#include "LidDrivenCavity.h"

LidDrivenCavity::LidDrivenCavity()
{	 
     T = 0.0;
     Nx = 0;
	 Ny = 0;
     Lx = 0.0;
     Ly = 0.0;
     Re = 0.0;
	 dt = 0.0;
	 vorticity= new double [int(Nx)*int(Nx)];
	 stream = new double [int(Nx)*int(Nx)];
	 CalV1_A = new double;
	 CalV1_b = new double [int(Nx-2)*int(Ny-2)];
	 //-------------------------------- CalV1_A[(Ny-2)*(Nx-2)*(Ny-2)*(Nx-2)] = {};
	 //-------------------------------- CalV1_b[(Ny-2)*(Nx-2)] = {};
	 cout << "Object is being created" << endl;
}

// test

LidDrivenCavity::~LidDrivenCavity()
{
	cout << "Object is being deleted" << endl;
}


/**  Using Doxygen to document a function
 * @brief Perform a single insert step in the insertion sort.
 * @param xlen Length of hte domain in the x-dirction
 * @param ylen Length of hte domain in the y-dirction
 *
 */
void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
	Lx = xlen;
	Ly = ylen;

}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
	Nx = nx;
	Ny = ny;
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
	dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
	T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
	Re = re;
}

void LidDrivenCavity::Initialise(double* omag, double* fi)
{
	vorticity = omag;
	stream = fi;
	double U = 1.0;
	double det_y = 1.0/(Ny-1.0);
	double det_x = 1.0/(Nx-1.0);
	//Initialise boundary condition
	for(int i=0; i<int(Nx*Ny); i++){
		if (i < int(Nx)){			
			omag[i] =(fi[i]-fi[int(Ny+i)])*(2.0/(det_y*det_y)); //Bottom			
			}			
		else if(i>int(Nx*(Ny-1))){
			omag[i] = fi[i-int (Ny)]-fi[i-2*int(Nx)]*(2.0/(det_y*det_y)) - (2.0*U/det_y); //Top		
		}
		else if (i%Nx ==1){		
			omag[i] = (fi[i]-fi[i+1]); //left
		}
		else if (i%Nx == 0){			
			omag[i] = (fi[i]-fi[i-1])*(2/(det_x*det_x)); // right
		}	
 	
	}
	
}

void LidDrivenCavity::CalVorticityT(double* A,double* b)
{
	
	CalV1_A = A;
	CalV1_b = b;
	//double det_y = 1.0/(Ny-1.0);
	//double det_x = 1.0/(Nx-1.0);
	int Ayn = Nx-1;
	int Axn = (Nx-2)*(Ny-2);
	cout<<Ny<<endl;
	int Id = Nx-1; //leading dimention  		
	for (int j = 0; j < Ayn+1; j++){
		for (int i = 0; i < Axn; i++){
		if (j== 0 && i > Nx-2){
			A[j*Id + i] = 1;
		}
		else if (j == Ayn){
			A[j*Id +i] = 2;
		}
		else if (j == Ayn-1){
			A[j*Id +i] = 3;
		}
		//else if (j == Ayn-1 && i%)
		
//			cout<<j%sub_An<<endl;
////			if (j%sub_An == 0 ){
////				CalV1_A[j*An+(j-1)] = 1;
////				CalV1_A[j*An+j] = 2;
////			}
////			else if (j%sub_An == sub_An-1){
////				CalV1_A[j*An+(j-1)] = 1;
////				
////				CalV1_A[j*An+(j-2)] = 3;
////			}	
////			else {
////				CalV1_A[j*An+(j-1)] = 1; //i,j
////				CalV1_A[j*An+j] = 2; //(i+1),j
////				CalV1_A[j*An+(j-2)] = 3; //(i+1),j
////			}
		}		
	
	}	
}

void LidDrivenCavity::CalVorticityTplus()
{		
}
void LidDrivenCavity::Integrate()
{
}
