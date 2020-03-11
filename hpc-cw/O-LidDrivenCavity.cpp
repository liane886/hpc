#include "LidDrivenCavity.h"
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <fstream>
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
	 N = 0;
	 v = {};
	 s = {};
	 CalVI_A = {};
	 BCforV_i = {};
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

void LidDrivenCavity::SetGridSize(int nx, int ny,int NumberofPoints)
{
	Nx = nx;
	Ny = ny;
	N = NumberofPoints;
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

void LidDrivenCavity::Initialise(double* omag, double* fi,double* V_i,double*V_j)
{
	v = omag;
	s = fi;
	BCforV_i = V_i;
	BCforV_j = V_j;
	//double* topBC = new double[Nx-2];
	//double* bottomBC = new double[Nx-2];
	//double* leftBC = new double[Nx-2];
	//double* rightBC = new double[Nx-2];
	double U = 1.0;
	double det_y = Ly/double (Ny-1);
	double det_x = Lx/double (Nx-1);
	
	
//	for (int i=0; i<Nx-2; i++){
//		topBC[i] = fi[i]-fi[i-1]*(2.0/(det_y*det_y)) - (2.0*U/det_y);
//		bottomBC[i] = (fi[i]-fi[i+1])*(2.0/(det_y*det_y));
//		leftBC [i]=  (fi[i]-fi[i+2*(Nx-2-i)+1])*(2.0/(det_x*det_x));
//		rightBC [i]= (fi[i]-fi[N-2*(Nx-2)+1+(N-i)])*(2/(det_x*det_x));
//	}
	

	//Initialise boundary condition
	for (int i=0; i< N; i++){
		
		if (i < Nx-2){
			if (i == 0) {
				omag[i] = 1; //Bottom	
				
			}
			else if (i == Nx-3) {
				omag[i] =1;//Bottom	
				//cout<<"11111111111111111"<<endl;
			}
			else {
				omag[i] = 1;
				//cout<<i<<endl;
			}
			
		}	
		else if (i == Nx-2) {
			
			omag[i] = 1;
		}
		else if (i == (Nx-2)*(Ny-3)) {
			omag[i] = 1;
		}
		else if (i >  (Nx-2)*(Ny-3)){
			
			 if (i ==  (Nx-2)*(Ny-3) +1){
				 omag [i] =1;
			 }
			else if(i == N-1) {
				omag[i] = 1;
			}
			else {
				omag[i] = 1;
			}
		}
		else omag[i] = 0;
			
	}
		
	for (int i=0; i<N; i++){
		
		if (i% Nx-2 ==0){
			V_i[i]  = 1;			
		}	
		else if (i% Nx-2 ==1) {
			V_i[i] = 1;
			//cout<<i<<endl;
		}
		
		else V_i [i] = 0;	
	}
	
	
	for (int i = 0; i<N; i++){
		
		if (i < Nx-2){
			V_j[i] = 1;			
		}	
		else if (i > N-(Nx-2)-1) {
			V_j[i] =1;
		}
		else V_j[i] = 0;
			
	}
	

//	for (int i=0; i< N; i++){
//		if (i < Nx-2){
//			if (i == 0) omag[i] = (fi[i]-fi[i+2*(Nx-2-i)+1])*(2.0/(det_x*det_x))-(fi[i]-fi[i+1])*(2.0/(det_y*det_y)); //Bottom	
//			if (i == Nx-3) omag[i] =(fi[i]-fi[i+2*(Nx-2-i)+1])*(2.0/(det_x*det_x)) -  fi[i]-fi[i-1]*(2.0/(det_y*det_y)) - (2.0*U/det_y);//Bottom	
//			else omag[i] = (fi[i]-fi[i+2*(Nx-2-i)+1])*(2.0/(det_x*det_x));
//			}	
//		if (i == Nx-2) omag[i] = (fi[i]-fi[i+1])*(2.0/(det_y*det_y));
//		if (i == (Nx-2)*(Ny-3)) omag[i] = fi[i]-fi[i-1]*(2.0/(det_y*det_y)) - (2.0*U/det_y);
//		if (i >  (Nx-2)*(Ny-3)){
//			if (i ==  (Nx-2)*(Ny-3) +1) omag [i] = (fi[i]-fi[N-2*(Nx-2)+1+(N-i)])*(2/(det_x*det_x)) - (fi[i]-fi[i+1])*(2.0/(det_y*det_y));
//			if(i == N-1) omag[i] = (fi[i]-fi[N-2*(Nx-2)+1+(N-i)])*(2/(det_x*det_x)) - fi[i]-fi[i-1]*(2.0/(det_y*det_y)) - (2.0*U/det_y);
//			else omag[i] = (fi[i]-fi[N-2*(Nx-2)+1+(N-i)])*(2/(det_x*det_x));
//		}
//			
//	}
//	for (int i=0; i<N; i++){
//		if (i% Nx-2 == 0){
//			V_i[i]  = (fi[i]-fi[i+1])*(2.0/(det_y*det_y));			
//		}	
//		else if (i% Nx-2 ==1) {
//			V_i[i] = fi[i]-fi[i-1]*(2.0/(det_y*det_y)) - (2.0*U/det_y);
//		}
//			
//	}
//	
//	for (int i = 0; i<N; i++){
//		if (i < Nx-2){
//			V_j[i] = (fi[i]-fi[i+2*(Nx-2-i)+1])*(2.0/(det_x*det_x));			
//		}	
//		else if (i > N-(Nx-2)-1) {
//			V_j[i] = (fi[i]-fi[N-2*(Nx-2)+1+(N-i)])*(2/(det_x*det_x));
//		}
//			
//	}
//
//
//
//
//	
//	for(int i=0; i<Nx*Ny; i++){
//		if (i < Nx){			
//			omag[i] =(fi[i]-fi[Ny+i])*(2.0/(det_y*det_y)); //Bottom			
//			}			
//		else if(i>Nx*(Ny-1)){
//			omag[i] = fi[i-Ny]-fi[i-2*Nx]*(2.0/(det_y*det_y)) - (2.0*U/det_y); //Top		
//		}
//		else if (i%Nx ==1){		
//			omag[i] = (fi[i]-fi[i+1]); //left
//		}
//		else if (i%Nx == 0){			
//			omag[i] = (fi[i]-fi[i-1])*(2/(det_x*det_x)); // right
//		}	
// 	
//	}
	
}

/**  Solve inter vortisity 
 * @brief Solve inter vortisity by y = A*x --> y    && using cblas-dsbmv (band& symmetric matrix)
 * @param 		A 			Matrix  A in the equation y = A*x 
 * @param 	VortiInter 		inter vorticity (y in the equation y = A*x)
 *
 */
 
void LidDrivenCavity::CalVorticityT(double* A, double* omag,double* fi)
{
	
	CalVI_A = A;	
//generate matrix A 
	double det_y = Ly/double (Ny-1);
	double det_x = Lx/double (Nx-1);
	//const int N = (Nx-2)*(Ny-2);    // 	Matrix dimension
	const int Ldh = Nx-1;       // leading diamention
	A[159] = 2*(det_y*det_y+det_x*det_x);
	for (int i = 1; i < N; ++i) {
		if (i < Nx-2){			
			A[i*Ldh + Ldh-2] = det_x*det_x;
			A[i*Ldh + Ldh-1] = 2*(det_y*det_y+det_x*det_x);
		}
		else {
			if (i%(Ldh-1) == 0){
				A[i*Ldh + Ldh-1] = 2*(det_y*det_y+det_x*det_x);
				A[i*Ldh        ] = det_y*det_y;					
			}
			else {
				A[i*Ldh + Ldh-1] = 2*(det_y*det_y+det_x*det_x);	
				A[i*Ldh + Ldh-2] = det_x*det_x;
				A[i*Ldh        ] = det_y*det_y;	
				
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
	double* Ans = new double[N];
	cblas_dcopy(N, omag, 1, Ans, 1); 
	cblas_dsbmv(CblasColMajor,CblasUpper,N, k, alpha,A,Ldh,fi,Inc,beta,Ans,Inc);
	
}

/**  Solve the equation 11
 * @brief Solve inter vortisity by y = A*x --> y    && using cblas-dsbmv (band& symmetric matrix)
 * @param 		A 			Pointer to vector of length N*(Nx-1). Matrix A in the equation y = A*x is stored in this way.
 * @param 	VortiInter 		Pointer to vector of length N & inter vorticity (y in the equation y = A*x)
 *
 */
 
void LidDrivenCavity::CalVorticityTplus(double* A,double* omag,double* fi,double* V_i,double*V_j)
{	
	//double *A2 = new double[(3*Nx+1)*(Ny-2)*(Nx-2)];
	const double beta = 0.0;
	const int Inc = 1; 	
	const int Ldh = Nx-1;
	const int k = Ldh -2;			
	// build the martix A to solve W_(i+1,j) - W_(i+1,j) using y = A* W
	double *MatrixA_i = new double[N*3];
	const int Ldh_MatrixA_i = 3;	
	MatrixA_i[2] = -1;
	for (int i = 1; i < N; ++i) {
		if (i%(Ldh-1) != 0){
			MatrixA_i[i*Ldh_MatrixA_i    ] =  1;
			MatrixA_i[i*Ldh_MatrixA_i + 2] = -1;				
		}				
    }
	// build the martix A to solve W_(i,j+1) - W_(i,j-1) using y = A* W
	const int Ldh_MatrixA_j = 2*(Nx-2)+1;
	double *MatrixA_j = new double[N*Ldh_MatrixA_j];
	
	MatrixA_j[Ldh_MatrixA_j-1] = -1;
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

	// Initialise bc
	double *S_i = new double[N];
	double *S_j = new double[N];	
	double *LHS_1 = new double[N];
	double *LHS_2 = new double[N];	
	double det_y = Ly/double (Ny-1);
	double det_x = 	Lx/double (Nx-1);
	const int KLU_i = Ldh_MatrixA_i -2;  // sub-diagonal & super-diagonal
	const int KLU_j =  Nx-2;
	const double alpha_i = 0.5*det_y*det_x*det_y;
	const double alpha_j = 0.5*det_y*det_x*det_x;

	//double eps;
	for (int j =0;j<(Nx-2)*(Ny-2); ++j){
		cout<<v[j];
		cout<<" ";
	}
	cout<<endl;	
	//double tol = 1e-08;
		cblas_dcopy(N, omag, 1, RhsAns, 1); 

		//cout << "Iteration " << k << endl;
		cblas_dsbmv(CblasColMajor,CblasUpper,N, k, 1/Re,A,Ldh,omag,Inc,beta,RhsAns,Inc);  //calculate  right hand side 
		
		cblas_dgbmv(CblasColMajor,CblasNoTrans,N,N,KLU_i,KLU_i,alpha_i,MatrixA_i,Ldh_MatrixA_i, //
						omag,Inc,beta,V_i,Inc);
		
		cblas_dgbmv(CblasColMajor,CblasNoTrans,N,N,KLU_i,KLU_i,alpha_i,MatrixA_i,Ldh_MatrixA_i,
						fi,Inc,beta,S_i,Inc);
						
		cblas_dgbmv(CblasColMajor,CblasNoTrans,N,N,KLU_j,KLU_j,alpha_j,MatrixA_j,Ldh_MatrixA_j,
						omag,Inc,beta,V_j,Inc);
						
		cblas_dgbmv(CblasColMajor,CblasNoTrans,N,N,KLU_j,KLU_j,alpha_j,MatrixA_j,Ldh_MatrixA_j,
						fi,Inc,beta,S_j,Inc);
		
				
		double da = (det_x*det_x*det_y*det_y);
		for (int i=0;i<N;i++){
			LHS_2 [i] = dt*(RhsAns[i]- S_j[i]*V_i[i]+ S_i[i]*V_j[i]) + da*omag[i];
		}

		//cblas_daxpy(N, 1, RhsAns, 1, LHS_2, 1);
	cblas_dswap(N,omag,Inc,LHS_2,Inc);
	

	
		//cblas_daxpy(N, -1, LHS_1, 1, LHS_2, 1);
		//cblas_daxpy(N, da, omag, 1, LHS_2, 1);

	//ofstream stream;
	//stream.open("/home/li/Desktop/hpc/hpc-cw/data.txt");
	//int Ldh = 3*Nx+1;
	//for(int i =0;i<(Nx-1);++i){

	//}
//	stream.close();

//		eps = cblas_dnrm2(N, r, 1);
//        cout << "eps: " << sqrt(eps) << " tol=" << tol << endl;
//        if (sqrt(eps) < tol) {
//            break;
//        }
	//Initialise(omag,fi,V_i,V_j);


	delete[] LHS_2;
	delete[] LHS_1;
	delete[] S_j;
	delete[] S_i;
	//delete[] V_i;
	//delete[] V_j;
	delete[] RhsAns;	
}

void LidDrivenCavity::Integrate()
{
	
}
