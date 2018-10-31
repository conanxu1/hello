#include <stdio.h>


#include "myfun.h"



double jifen(double a,double b,int n,double (*fff)(double))
{	double h = (b - a)/n, T = 0; 
	for(int i = 1; i < n; i++) 
		T += (*fff)(a + i * h); 
	
	return  h/2*(fff(a)+T*2+fff(b));
	
}

double ode(double a,double b,int n,double (*fff)(double),)
{/*a初值,b终值,n分点,fff连续函数   */



}














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