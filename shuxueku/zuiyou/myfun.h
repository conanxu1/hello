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



double  *ode1(double a,double b,int n,double *(*fff)(double ,double* ),double *x0,int dim);

int  jiecheng(int n);



void shuchuf(float *p,int a,int b);
void shuchud(double *p,int a,int b);
void shuchuz(double _Complex *p,int a,int b);


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


int phi(double* f,int N,double t,double T);

int initxishu(double* D1,double* D0,double* D,int N);
