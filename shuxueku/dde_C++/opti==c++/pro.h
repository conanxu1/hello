#include <iostream>
#include "time.h"
#include <chrono>
#include <string>
// #include <stdlib.h>
#include <stdio.h>

#include <cstring>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "nlopt.h"



using namespace std;


//#define noeigen 1
#ifdef noeigen
	typedef double* Rn;
	typedef double R;
	double* Rninit(int nnnn);
#else
	#include <eigen3/Eigen/Dense>
	using namespace Eigen;
	typedef VectorXd Rn;
	typedef double R;
	VectorXd Rninit(int nnnn);
#endif


typedef Rn (*dae_f)(Rn x ,      Rn u     ,double        t );

typedef R (*opt_int)(Rn x ,      Rn u     ,double        t );
typedef R (*opt_phi)(Rn xt0 ,     double t0 ,Rn xtf ,     double tf     );
typedef Rn ( *pfpx)(Rn x ,      Rn u     ,double        t ,int dimx,int dimu);


R L_zero(Rn x ,      Rn u     ,double        t );


struct bolza_con
{
	opt_phi psi;
	opt_int L;
};

struct point_con
{
	dae_f c;
};



class Pro;
class Sol;



struct individual
{
	Rn state;
	R   val;
	R   fitness;
	R   prob;
	R   sum_prob;
	int cho_times;
};





// typedef R (*fff)(Rn x ,      Rn u     ,double        t );

class Sol
{
	private:
		int dimx;
		int dimu;
		double t0;
		double tf;
		
		
		//原方程
		dae_f fxut;
		opt_int opt_L;
		opt_phi opt_Phi;
		pfpx pLpu;
		pfpx pLpx;
		

		
		//各种约束
		int numB;
		//终端约束
		int numC;
		//点约束
		
		bolza_con	 	*b_c;
		point_con 	*p_c;
		
		int N;
		double  h;
		
		

		
		
		Pro *tPro;
		
		Rn tk;
		

		
		
		
		
		
		
	public:
		Sol(Pro* ttPro):tPro(ttPro){}
		
		void set();

		
		
		void display();
		void sol(int N);
		
		R testf(Rn);
		
		R objfun( Rn theta, Rn grad);
		
		
		
		Rn x_c(dae_f temf,Rn xk,Rn uk,double tk,Rn xk1,Rn uk1,double tk1);
		
		//约束函数
		Rn hk_th( Rn theta,int k);	
		R gi_th(  Rn theta);
		R c1jk(    Rn theta);
	
	
};



class Pro
{
	private:
		int dimx;
		int dimu;
		double t0;
		double tf;
		dae_f fxut;
		opt_int j_opt;
		opt_phi phi_opt;
		
		
		
		int numB;
		//终端约束
		int numC;
		//点约束
		
		bolza_con	 	*b_c;
		point_con 	*p_c;
		
		
		// dae_f *nn;
		// nn=new dae_f[5];
		
		pfpx pLpu;
		pfpx pLpx;
	
	
	public:
		Pro(int tdimx,int tdimu,double tt0,double ttf ,int nnumB,int nnumC  ):dimx(tdimx),dimu(tdimu),t0(tt0),tf(ttf), numB(nnumB),numC(nnumC) {
			b_c=new 	bolza_con[numB];
			p_c=new 	point_con[numC];
			
			
		}
			
		~Pro(){}
			
		
		void set_daef(dae_f fxut );
		void set_opt_int(opt_int jj_opt );
		void set_opt_phi(opt_phi tphi_opt );
		
		void set_opt_pLp(	pfpx pLpx,	pfpx pLpu);

		friend void Sol::set(  );
};




class ga_opti
{
private:
	int mpop_size;
	int max_itr;
	int dim;
	
 
	Rn lb;
	Rn ub;
	
	struct individual *group;
	double sum_fit;
	
	
public:	

	ga_opti(int tdim,int tmpop_size,int tmax_itr): dim(tdim),mpop_size(tmpop_size),	 max_itr(tmax_itr){}

	void initialize();
	void set_lbub(Rn llb,Rn uub);
	void gadisplay();
	double myfit(double tval);
	
	void comp_val_fitness_prob();
	void chooose();
	 
};









// class Alg
// {
	// private:
		// int dimx;
		// int dimu;
		// double t0;
		// double tf;

// }






