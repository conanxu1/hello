#include <stdio.h>


#include "myfun.h"



double jifen(double a,double b,int n,double (*fff)(double))
{	double h = (b - a)/n, T = 0; 
	for(int i = 1; i < n; i++) 
		T += (*fff)(a + i * h); 
	
	return  h/2*(fff(a)+T*2+fff(b));
	
}


/*
double  ode1(double a,double b,int n,double* (*fff)(double t,double *x),double *x0,int dim)
{/*a初值,b终值,n分点,fff连续函数   


double h = (b - a)/n;

double t0=a;
double q1[dim],q2[dim],q3[dim],q4[dim];






for(int i = 1; i < n; i++)
{
q1=f(tk,yk)
q2=f(tk+b/2,yk+b/2*q1)
q3=f(tk+b/2,yk+b/2*q1)
q4=f(tk+b,yk+b*q1)
yk+1=yk+b/6*(q1+2*q1+2*q3+q4)
}





}


*/











void shuchuf(float *p,int a,int b)
{
	for(int i=1;i<=a;i++)
	{	for(int j=1;j<=b;j++)
		{printf("%f,",p[(i-1)*b+j-1]);}
	printf("\n");
	}
	
}

void shuchud(double *p,int a,int b)
{
	for(int i=1;i<=a;i++)
	{	for(int j=1;j<=b;j++)
		{printf("%f,",p[(i-1)*b+j-1]);}
	printf("\n");
	}
	
}