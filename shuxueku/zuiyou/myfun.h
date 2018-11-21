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

#define PI 3.14159265358979323846
double jifen(double a,double b,int n,double (*f1)(double));

double _Complex lamw(double _Complex a,int k);

int st(int p,int k);



double *lianxu(double t,double *x);



double  *ode1(double a,double b,int n,double *(*fff)(double ,double* ),double *x0,int dim);

int  jiecheng(int n);



void shuchuf(float *p,int a,int b);


void shuchud(double *p,int a,int b);


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


void svd2();
