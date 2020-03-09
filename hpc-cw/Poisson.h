
class Poisson
{
	
private:
    int    Nx;
    int    Ny;
	double Lx;
    double Ly;
	double* PoissonMA = nullptr;
	double* Pv = nullptr;
    double* Ps = nullptr;
	
public: 
	Poisson();
	~Poisson();
	void ComputeStreamFunction(double* fi, double* omag,double* A2);
	void SetGridSize (int nx,int ny);
	void SetDomainSize (double xlen,double ylen);
};