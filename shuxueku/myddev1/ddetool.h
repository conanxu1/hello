#include <iostream>
#include <chrono>
#include <string>
// #include <stdlib.h>
#include <cstring>
#include <Eigen/Dense>



#include <cmath>
#include <iomanip>

using namespace std;
using namespace Eigen;



typedef double* Rn;

typedef Rn (*C_tau)(double t);

typedef Rn (*C_R)(double t);

typedef VectorXd (*Vec_C_R)(double t);
typedef VectorXd   (*Vec_Cphi)(double t0,double t   );


typedef Rn (*Cphi)(double t0,double t);


class RK4_ddesolver ;

typedef Rn ( RK4_ddesolver  :: *C_tau2)(double t);


typedef Rn (*functional)(double t,Rn x ,      Rn* xdelay      ,C_R fu_t        );







/* dx=f(t,x,u,x(~),u(~))
*	dx x in R^dimx u in R^dimu
*	x(~) in [-tau,0] -> R^dimx Cab  C_tau
*  u(~) in [-tau,0] -> R^dimu
*
*
*/



// class abs_ddesolver{
	// int dimx;
	// int dimu;
	
	

	// public:
		
		
		// // hcout();
		// // virtual void  solve(  dde_f f )=0;
	
// };


 
class RK4_ddesolver 
{
	private:
		int dimx;
		int dimu;
		int tauN;
		double taumax;
		double t0,tf;	
		Cphi phi;
		C_R    ut;
		functional   Ftx;
		Rn tauk;
		
		Rn solt;
		Rn* solx;
		
		
		
		
		
		
		Vec_Cphi vec_phi;
		Vec_C_R vec_u;

		int stepN=100000+1;
		double h0;
		
		
		Matrix<double, Dynamic, Dynamic> *Ai;
		Matrix<double, Dynamic, Dynamic> *Bi;
	
		VectorXd *Vec_solx;
	
	
	
		
	public:
		RK4_ddesolver(int dix,int diu):dimx(dix),dimu(diu){}
		
		~RK4_ddesolver(){
			// cout<<"here";
			// for(int i=0;i<this->stepN;i++)
			// {
				// delete (this->solx)[i];
			// }
			// delete solx;			
		}
		
		void display();
		void setTauN(int tauN);
		void setTauk(Rn tauk);
		void set_t0_tf(double t0,double tf);
	
		
		
		void set_phi(Cphi phi);
		void set_ut(C_R ut);
		
		void set_eq(functional   Ftx );
	
		
		void set_Ai( Matrix<double, Dynamic, Dynamic>* tAi  );
		void set_Bi( Matrix<double, Dynamic, Dynamic>* tBi  );
		
		void set_vecphi(Vec_Cphi mphi);
		void set_vecut(Vec_C_R mut);
		
		
		void display_AiBi();

		
	
		void solve(  );
		void solveAiBi( );
		void display_sol();
		void display_solVec();
		VectorXd fAiut(double t , VectorXd x ,   VectorXd *xdelay    , Vec_C_R fu);
		
		// VectorXd f2(double t  ,VectorXd x, Vec_C_R fu);
		VectorXd f2(double t   ,VectorXd x,VectorXd* xdelay, Vec_C_R fu);
		
		
	private:
		
		void calx( C_R f,double t,Rn as        );
		void calu( C_R f,double t,Rn as        );
		Rn* chi_t(double t);
		VectorXd*   vec_chi_t(double t);
		
		
		
		void onestep(Rn xn1,Rn xn ,double hh,Rn ff   );

	
		
		
			
 };







// class RK4_ddesolver
// {
	// private:
		// int g;
	
	// public:
		 // RK4_ddesolver(int i) : g(i)
		// {
			
			
		// }
// };


