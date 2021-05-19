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




#define USE_OPENCV

#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#endif




//=======================
// #define noeigen 1	
 
//定义     不使用Eigen 库
//不定义 使用
// cmake add_definitions
//#define USE_EIGEN 1

#ifndef USE_EIGEN
	typedef double* Rn;
	typedef double R;
	typedef double* Rmn;
	
#endif



#ifdef USE_EIGEN
	#include <Eigen/Dense>
	using namespace Eigen;
	typedef VectorXd Rn;
	typedef MatrixXd Rmn;
	typedef double R;

#endif

#ifdef USE_EIGEN
void double2mat(double *d_a,Rmn &m_a,int m,int n);
void mat2double(Rmn &m_a,double *d_a,int m,int n);
void double2vex(double *d_a,Rn &m_a,int m);
void vec2double(Rn &m_a,double *d_a,int m);

#endif









Rn Rninit(int nnnn);
Rmn Rmninit(int m,int n);




void Rmn_copy(double* tm, Rmn (&ptM),int m,int n);
void Rn_copy(double* tm,Rn  (&tM),int n);
 


typedef Rn (*dae_f)(Rn x ,      Rn u     ,double        t );

typedef R (*opt_int)(Rn x ,      Rn u     ,double        t );
typedef R (*opt_phi)(Rn xt0 ,     double t0 ,Rn xtf ,     double tf     );
typedef Rn ( *pfpx)(Rn x ,      Rn u     ,double        t ,int dimx,int dimu);

typedef Rn (*Rn_f)(double t);


//工具函数声明
template <typename T_Rn>
void print_vec(T_Rn x,int n);




//具体问题中的一些数据结构


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







//抽象类：微分方程求解器
class DAE_Solver
{
	
public:
	// // virtual void sol()=0;//不需要在基类中给出函数的函数体
	virtual void sol()=0;
	
};




//欧拉法
class Euler_Ode_Sol  : public DAE_Solver
{
private:
		int dimx;
		int dimu;
		double t0;
		double tf;
		int lenN;
		double h;
		dae_f fxut;	
		Rmn tA;
		Rmn tB;
		Rn_f tu;
		Rn x0;
		
		


		#ifdef SMALL_SCALE
		Rn* solx;
		R*   solt;
		#endif
	 	
		#ifdef LARGE_SCALE
		R*   solx;
		R*   solt;
		#endif
	 	
    //solx=[x1,x2,...,xn].';

		
		
		
public:
		
		
		Euler_Ode_Sol(int tdimx,int tdimu,double tt0,double ttf,int tlenN):dimx(tdimx),dimu(tdimu),t0(tt0),tf(ttf),lenN(tlenN){}
		
		void set(Rmn tA,Rmn tB,Rn_f tu,Rn x0);
		void sol( );
		
		
		// // // void display();
		// // // void sol(int N);
		


};





//微分方程求解


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





//约束优化问题

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
		
		void set_opt_pLp(pfpx pLpx,	pfpx pLpu);

		friend void Sol::set(  );
};

























//遗传算法




struct individual
{
	Rn state;
	R   val;
	R   fitness;
	R   prob;
	R   sum_prob;
	int cho_times;
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



void test();





// class Alg
// {
	// private:
		// int dimx;
		// int dimu;
		// double t0;
		// double tf;

// }







