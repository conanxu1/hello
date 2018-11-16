#pragma once

#define AAA 1



#ifdef  AAA
#include <openblas/cblas.h>


#else
#include <cblas.h>


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



