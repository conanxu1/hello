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

int erciw(
		double *H,		//hessian
		double *h,		//原问题grad
		double *be,		//b   等式约束
		double *Ae,		//系数
		double *bi,		////b   不等式约束
		double *Ai,
		int dim,		//问题的维数
		int e,			//等式个数
		int ie			//不等式个数
);



int erci(
		double *H,		//hessian
		double *h,		//原问题grad
		double *be,		//b   等式约束
		double *Ae,		//系数
		double *bi,		////b   不等式约束
		double *Ai,
		int dim,		//问题的维数
		int e,			//等式个数
		int ie,
		double *xk);
		
int myqp(
		double *G,		//hessian
		double *A,		//grad
		double *gk,
		double *b,
		int dim,		//G维数
		int   e	);		//A的列数
		

int mychol(
		double *L,		//输入 返回 下三角
		int dim,
		double *Lz		//转好
	);		//维数