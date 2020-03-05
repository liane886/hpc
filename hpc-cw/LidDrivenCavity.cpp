#include "LidDrivenCavity.h"
#include <cstdlib>
#include <iomanip>
#include <cmath>
//#include "cblas.h"
#include "/home/li/Desktop/header/cblas/CBLAS/include/cblas.h"
//#include "cblas_f77.h"
#define F77NAME(x) x##_

// This tells the C++ compiler to use C-style naming of symbols for these
// functions, as opposed to the "mangled" symbols used in C++. The need to
// mangle symbol names is because C++ supports overloading the same function
// name with different parameters and data types so needs to be able to
// distinguish them by name in the object code.
extern "C" {
    double F77NAME(ddot) (const int& n, 
                          const double *x, const int& incx,
                          const double *y, const int& incy);
}
double cblas_dnrm2(const int N, const double *X, const int incX);

LidDrivenCavity::LidDrivenCavity()
{	 
     T = 0.0;
     Nx = 0;
	 Ny = 0;
     Lx = 0.0;
     Ly = 0.0;
     Re = 0.0;
	 dt = 0.0;
	 vorticity= new double [Nx*Nx];
	 stream = new double [Nx*Nx];
	 CalVI_A = new double;
	 CalVI_x = new double [(Nx-2)*(Ny-2)];
	 CalVI_y = new double [(Nx-2)*(Ny-2)];
	 //-------------------------------- CalV1_A[(Ny-2)*(Nx-2)*(Ny-2)*(Nx-2)] = {};
	 //-------------------------------- CalV1_b[(Ny-2)*(Nx-2)] = {};
	 cout << "Object is being created" << endl;
}

// test

LidDrivenCavity::~LidDrivenCavity()
{
	cout << "Object is being deleted" << endl;	// print information
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
	double det_y = 1/(Ny-1);
	double det_x = 1/(Nx-1);
	//Initialise boundary condition
	for(int i=0; i<Nx*Ny; i++){
		if (i < Nx){			
			omag[i] =(fi[i]-fi[Ny+i])*(2.0/(det_y*det_y)); //Bottom			
			}			
		else if(i>Nx*(Ny-1)){
			omag[i] = fi[i-Ny]-fi[i-2*Nx]*(2.0/(det_y*det_y)) - (2.0*U/det_y); //Top		
		}
		else if (i%Nx ==1){		
			omag[i] = (fi[i]-fi[i+1]); //left
		}
		else if (i%Nx == 0){			
			omag[i] = (fi[i]-fi[i-1])*(2/(det_x*det_x)); // right
		}	
 	
	}
	
}

/**  Solve inter vortisity 
 * @brief Solve inter vortisity by y = A*x --> y    && using cblas-dsbmv (band& symmetric matrix)
 * @param 		A 			Matrix  A in the equation y = A*x 
 * @param 	VortiInter 		inter vorticity (y in the equation y = A*x)
 *
 */
 
void LidDrivenCavity::CalVorticityT(double* A, double* VortiInter,double* streamInter)
{
	
	CalVI_A = A;	
	CalVI_x = VortiInter;
	CalVI_y = streamInter;
//generate matrix A 
	double det_y = 1.0/double (Ny-1);
	double det_x = 1.0/double (Nx-1);
	//const int N = (Nx-2)*(Ny-2);    // 	Matrix dimension
	const int Ldh = Nx-1;       // leading diamention
	const int Axn = (Nx-2)*(Ny-2);
	 		
//	for (int j = 0; j < Ayn ; j++){
//		for (int i = 0; i < Axn ; i++){ //
//			if (j == 0 && i > Nx-3){
//				A[j*Axn +i] = det_y*det_y;  // for superdiagonal
//				//cout<<j*Id +i<<endl;
//			}
//			else if (j == Ayn-1){
//				A[j*Axn +i] = 2*(det_y*det_y+det_x*det_x);  // diagonal
//				cout<<j*Axn +i<<endl;
//			}
//			else if (j == Ayn-2 && i%(Nx-2) != 0){
//				A[j*Axn +i] = det_x*det_x;   // for first superdiagonal
//			}
//		}		
//	}

	A[159] = 2*(det_y*det_y+det_x*det_x);
	for (int i = 1; i < N; ++i) {
		if (i < Nx-2){			
			A[i*Ldh + Ldh-2] = det_x*det_x;
			A[i*Ldh + Ldh-1] = 2*(det_y*det_y+det_x*det_x);
		//cout<<i*Ayn + Axn-1<<endl;
	
		}
		else if (i < N){
			if (i%(Ldh-1) == 0){
				A[i*Ldh + Ldh-1] = 2*(det_y*det_y+det_x*det_x);
				A[i*Ldh        ] = det_y*det_y;					
			}
			else {
				A[i*Ldh + Ldh-1] = 2*(det_y*det_y+det_x*det_x);
				//cout<<i<<endl;	
				A[i*Ldh + Ldh-2] = det_x*det_x;
				A[i*Ldh        ] = det_y*det_y;		
				//cout<<i*Ayn + Ayn-1<<endl;
			}
		}				
    }

	
//	for (int i = 0; i < Ldh; ++i) {                                               //////print matrix A
//        for (int j = 200; j > 70; --j) {
//            cout <<setprecision(3)<< A[j*Ldh+i] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;

// solve the equation using blas (matrix*array)
	const double beta = 0.0;
	const double alpha = -1.0;
	const int Inc = 1; 
	const int k = Ldh -2;	
	cblas_dsbmv(CblasColMajor,CblasUpper,N, k, alpha,A,Ldh,streamInter,Inc,beta,VortiInter,Inc);

}

/**  Solve the equation 11
 * @brief Solve inter vortisity by y = A*x --> y    && using cblas-dsbmv (band& symmetric matrix)
 * @param 		A 			Pointer to vector of length N*(Nx-1). Matrix A in the equation y = A*x is stored in this way.
 * @param 	VortiInter 		Pointer to vector of length N & inter vorticity (y in the equation y = A*x)
 *
 */
 
void LidDrivenCavity::CalVorticityTplus(double* A,double* VortiInter,double* streamInter)
{	
	
	const double beta = 0.0;
	const int Inc = 1; 	
	const int Ldh = Nx-1;
	const int k = Ldh -2;	
	
	
	// build the martix A to solve W_(i+1,j) - W_(i+1,j) using y = A* W
	double *MatrixA_i = new double[N*(Nx-1)];
	const int Ldh_MatrixA_i = 3;	
	MatrixA_i[2] = -1;
	for (int i = 1; i < N; ++i) {
		if (i%(Ldh-1) != 0){
			A[i*Ldh_MatrixA_i    ] =  1;
			A[i*Ldh_MatrixA_i + 2] = -1;				
		}				
    }

	// build the martix A to solve W_(i,j+1) - W_(i,j-1) using y = A* W
	double *MatrixA_j = new double[N*(Nx-1)];
	const int Ldh_MatrixA_j = 2*(Nx-1)-1;
	MatrixA_j[159] = -1;
	for (int i = 1; i < N; ++i) {
		if (i < Nx-2){			
			MatrixA_j[i*Ldh_MatrixA_j + Ldh_MatrixA_j-1] = -1;	
		}
		else {			
				MatrixA_j[i*Ldh_MatrixA_j + Ldh_MatrixA_j-1] = -1;
				MatrixA_j[i*Ldh_MatrixA_j        ] = 1;								
		}
	}
	double *RhsAns = new double [N]; // the answer of the right hand side of the equation 11
	double *LhsStreamAns_i = new double[N]; 
	double *LhsStreamAns_j = new double[N];
	double *LhsVortiAns_i = new double[N];
	double *LhsVortiAns_j = new double[N];
	
	
	double det_y = 1.0/double (Ny-1);
	double det_x = 1.0/double (Nx-1);
	const int KLU_i = Ldh_MatrixA_i -2;  // sub-diagonal & super-diagonal
	const int KLU_j =  Nx-2;
	const double alpha_i = 0.5*det_y*det_x*det_y;
	const double alpha_j = 0.5*det_y*det_x*det_x;

	//double eps;
	//int k;               //iteration counter
	double tol = 1e-08;
	//do {
		
		//cout << "Iteration " << k << endl;
		cblas_dsbmv(CblasColMajor,CblasUpper,N, k, 1/Re,A,Ldh,VortiInter,Inc,beta,RhsAns,Inc);
		cblas_dgbmv(CblasColMajor,CblasNoTrans,N,N,KLU_i,KLU_i,alpha_i,MatrixA_i,Ldh_MatrixA_i,
						VortiInter,Inc,beta,LhsVortiAns_i,Inc);
						
		cblas_dgbmv(CblasColMajor,CblasNoTrans,N,N,KLU_i,KLU_i,alpha_i,MatrixA_i,Ldh_MatrixA_i,
						streamInter,Inc,beta,LhsStreamAns_i,Inc);
						
		cblas_dgbmv(CblasColMajor,CblasNoTrans,N,N,KLU_j,KLU_j,alpha_j,MatrixA_j,Ldh_MatrixA_j,
						VortiInter,Inc,beta,LhsStreamAns_j,Inc);
						
		cblas_dgbmv(CblasColMajor,CblasNoTrans,N,N,KLU_j,KLU_j,alpha_j,MatrixA_j,Ldh_MatrixA_j,
						streamInter,Inc,beta,LhsStreamAns_j,Inc);
		
		double *LHS_1 = new double[180];
		double *LHS_2 = new double[180];
		for(int i = 0; i < 180; ++i){
			LHS_2[i] = 3;
			LHS_1[i] = 2;
		}
		cblas_ddot(180,LHS_2,Inc,LHS_1,Inc); //LHS 1
		
		cout<<LHS_1[1]<<endl;
		cout<<LHS_2[1]<<endl;
		
		cblas_ddot(N,LhsStreamAns_j,Inc,LhsVortiAns_i,Inc); //LHS 1
		cblas_ddot(N,LhsStreamAns_i,Inc,LhsVortiAns_j,Inc); //LHS 2
		cblas_daxpy(N, 1, RhsAns, 1, LhsVortiAns_i, 1);
		cblas_daxpy(N, -1, LhsVortiAns_i, 1, LhsVortiAns_j, 1);
		
//		eps = cblas_dnrm2(N, r, 1);
//        cout << "eps: " << sqrt(eps) << " tol=" << tol << endl;
//        if (sqrt(eps) < tol) {
//            break;
//        }
//		
	//	k++;
	//}while (k<500);
	
}

void LidDrivenCavity::Integrate()
{
	
}
