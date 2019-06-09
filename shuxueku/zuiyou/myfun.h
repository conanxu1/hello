#pragma once

#define AAA 1



#ifdef  AAA
#include <cblas.h>
#include <lapacke.h>
#include <lapacke_config.h>
#include <lapacke_utils.h>



#else
#include <cblas.h>
#include <lapacke_config.h>
#include <lapacke_utils.h>
#include "fftw3.h"


#include <lapacke.h>

#endif





















double jifen(double a,double b,int n,double (*f1)(double));

double _Complex lamw(double _Complex a,int k);


int lamwm(double _Complex* A,int k,int dim,double _Complex*	 w);

void lamwhb(double _Complex* A,double _Complex* w,int dim);

int expm(double _Complex* A,int dim);



double mubiao(double _Complex* A,double _Complex* Ad,double _Complex* Q,double h,int dim);

int st(int p,int k);



double *lianxu(double t,double *x);



double  *ode1(double a,double b,int n,double *(*fff)(double ,double *),double *x0,int dim,double *jieguo);

int  jiecheng(int n);



void shuchuf(float *p,int a,int b);
void shuchud(double *p,int a,int b);
void shuchuz(double _Complex *p,int a,int b);
void shuchui(int *p,int a,int b);

double f1(double x);


double _Complex lnk(double _Complex a,int k);
 
double* danwei(int m,int n);

void sduijiao(double* p,int m,int n,double* Q2,double* V1,double* B);

void svd(double* p,int m,int n);



void qrweiyi(double* A,int m,int n);

void QR(double* A,int m,int n,double* Q1,double* R1);




double dot(double* a,double* b,int dim);
void house(double* a,int dim,int i);

void xqu(double* A,double* u,int i,int m,int n);
void xqu2(double* A,double* u,int i,int m,int n);

double* fuzhi1(double h1,int i,int m);


double* duqu(char *p,int n);


void inrk(double ***Yni,double ***Kni,double **yn,int r,int s,int m,double tau);

double* dfini(double t);
double* fini(double t);

void drk();

// void pq();
void svd2();




void qiuQ(double _Complex* A,double _Complex* Ad,double _Complex* Q,double h,int dim);

void nugrad(double _Complex* A,double _Complex* Ad,double _Complex* Q,double h,int dim,double _Complex* gra);



void jiangjie();

int m0618();
int phi(double* f,int N,double t,double T);

int initxishu(double* D1,double* D0,double* D,double* S,int N,double T,double tau);

double *exmf(double t,double *x,double *sigma,int dim);
double shixing(double t,double a,double b);
int changema(double* x,int size_x);
void mypinv();





int fft(double _Complex  * ai,double _Complex  * ao,int N);
void ni(double *A,int dim);



int ooo();


/*
int erciw(
		double *H,		//hessian
		double *h,		//原问题grad
		double *be,		//b   等式约束
		double *Ae,		//系数
		double *bi,		////b   不等式约�?
		double *Ai,
		int dim,		//问题的维�?
		int e,			//等式个数
		int ie			//不等式个�?
);
*/


int erci(
		double *H,		//hessian
		double *h,		//原问题grad
		double *be,		//b   等式约束
		double *Ae,		//系数
		double *bi,		////b   不等式约�?
		double *Ai,
		int dim,		//问题的维�?
		int e,			//等式个数
		int ie,
		double *xk);
		
int myqp(
		double *G,		//hessian
		double *A,		//grad
		double *gk,
		double *b,
		int dim,		//G维数
		int   e	);		//A的列�?
		

int mychol(
		double *L,		//输入 返回 下三�?
		int dim,
		double *Lz		//转好
	);		//维数

int qxt(
		double *H,		//hessian
		double *h,		//原问题grad
		double *be,		//b   等式约束
		double *Ae,		//系数
		double *bi,		////b   不等式约�?
		double *Ai,

		int dim,		//问题的维�?
		int e,		
		int ie			//不等式个�?原问题的  2e+ie)
		);

//Rn->R
typedef double (*RnR)(double *x,int n);
typedef double (*RR)(double x);
typedef double *(*RnRn)(double *x,int n);

//维数作为全局�?不出�?
typedef double (*yRnR)(double *x);
typedef double (*yRR)(double x);
typedef double* (*yRnRn)(double *x);



//二级指针向量列表

typedef double (*lyRnR)(double **x);

typedef double (*lyRR)(double x);

typedef double** (*lyRnRn)(double **x);

typedef double* (*lyRnRnk)(double **x);



//计算需要终端约�?分对xu 偏导的雅可比 和对时间的偏�?
typedef double* (*piandao)(double **x,int type);


double fxutRnR(lyRnR f,double *X,int *listx,int gex);



//指标转换

int v2tt(int v,int dim,int dimu,int N);
int v2qp(int v,int dim,int dimu,int N);
int v2qx(int v,int dim,int dimu,int N);
int v2s(int v,int dim,int dimu,int N);




//求解非线性约束优化问?
int yueshu(RnR f,RnR *fk);

//梯度
double gk(double *x,int n);







int xxwg(double *AM,double *A,int m,int n);








//勒让德高斯积分的  勒让德节?
int lgd(double *x,int n);


//勒让德高斯系数Ak  Dki   A0...An   导数系数   d Li(xk)

int lg_AD(double *xk,double *Ak,double *Dki,int n);





//追赶法解方程

int zhui(double *A,		//系数
		double *d,		//右侧
		int n,
		double *jie);



int zuoyongyu(int **a);








 

int zuoyongyu(int **a);



/*--------------------------*/
int xishu(lyRnR PHI,       //目标函数中的终端
		lyRnR g,     //目标函数中的被积函数
		lyRnRnk f,
		

		lyRnR *cek,
		int numcek,		//逐点等式约束
		
		lyRnR *phii,
		int numphii,	//mayer终端等式约束的个?
		
		lyRnR *cik,
		int numcik,		//逐点不等式约? cik<=0
		
		
		lyRnR *psii,
		int numpsii,	//不等式终?


		piandao gPHI,			//PHI(x,t) 返回PHI_x,PHI_t
		
		
		piandao gg,	 //被积函数的梯?
		//gg(**(x,u,t))  返回(**( gg_kesi,gg_t ))		
		
		piandao gf,	 //状态方程函数的梯度
		
		piandao *gcek,	//逐点的梯度函
		
		piandao *gphi,
					
		piandao *gcik,
		
		piandao *gpsi,				//mayer 的梯
		
	
		
		int dimx,
		int dimu,
		int Ntau,  //几阶勒让德方
		double *tauk, //勒让德点
		double *wk, //高斯勒让德积分系
		double *Dki, //导数系数
		double *zuiyouX, //最优结
	                          	//  double *H, H由这些和子问题等等进行修? Bk+1
		double *h,
		double *be,	
		double *Ae,	
		double *bi,
		double *Ai,
		double *XUk,
		double tk,
		double *x0,
		double t0
);

int X2xutk(double *X,double **x,int dimx,int dimu,double *tk,int k);
int cshi(double *x,int n);



double kzPHI(double **xftf);
double kzg(double **xut);
double* kzf(double **xut);

double kzphii1(double **xftf);
double kzphii2(double **xftf);
double kzcik1(double **xut);
double kzcik2(double **xut);
double* kzgPHI(double **xftf,int flag);

double* kzgg(double **xut,int flag);	
double* kzgf(double **xut,int flag);	
double* kzgphi1(double **xftf,int flag);	
double* kzgphi2(double **xftf,int flag);	
double* kzgcik1(double **xut,int flag);
double* kzgcik2(double **xut,int flag);	



