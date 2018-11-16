#include <stdio.h>

#include <stdlib.h>
#include "myfun.h"
#include <complex.h>

#include <math.h>



double jifen(double a,double b,int n,double (*fff)(double))
{	double h = (b - a)/n, T = 0; 
	int i;
	for(i = 1; i < n; i++) 
		T += (*fff)(a + i * h); 
	
	return  h/2*(fff(a)+T*2+fff(b));
	
}



double  *ode1(double a,double b,int n,double *(*fff)(double ,double *),double *x0,int dim)
{/*a初值,b终值,n分点,fff连续函数*/   

FILE *logf;
double h = (b - a)/n;

double t0=a,tk,jieguo[dim*n],II[dim];



double* yk =(double *)malloc(dim*sizeof(double));
double* q1 =(double *)malloc(dim*sizeof(double));
double* q2 =(double *)malloc(dim*sizeof(double));
double* q3 =(double *)malloc(dim*sizeof(double));
double* q4 =(double *)malloc(dim*sizeof(double));


double* tem=(double *)malloc(dim*sizeof(double));







//*a要分配内存 或者定义时候分配
//double *dy =(double *)malloc(size*sizeof(double));






int i;

for(i = 0; i < dim; i++)
{
	II[i]=1;	
}

logf=NULL;


if((logf = fopen("log.txt" , "w+")) == NULL)
		{
			printf("Cannot create/open file");
			exit(1);		
			
		}
	


int weizhi=0;


tk=t0;
yk=x0;


for(i = 0; i < n; i++)
{
q1=lianxu(tk,yk);




cblas_daxpby(dim, 1, yk, 1, 0, tem, 1);
cblas_daxpby(dim,h/2, q1, 1, 1, tem, 1);
q2=lianxu(tk+h/2,tem);
//lianxu(tk+h/2,yk+h/2*q1);




cblas_daxpby(dim, 1, yk, 1, 0, tem, 1);
cblas_daxpby(dim,h/2, q2, 1, 1, tem, 1);
//q3=lianxu(tk+h/2,yk+h/2*q1);
q3=lianxu(tk+h/2,tem);



cblas_daxpby(dim, 1, yk, 1, 0, tem, 1);
cblas_daxpby(dim,h, q3, 1, 1, tem, 1);

//q4=lianxu(tk+h,yk+h*q1);
q4=lianxu(tk+h,tem);




cblas_daxpby(dim, 2, q2, 1, 1, q1, 1);
cblas_daxpby(dim, 2, q3, 1, 1, q1, 1);
cblas_daxpby(dim, 1, q4, 1, 1, q1, 1);



//q1+2*q2+2*q3+q4





cblas_daxpby(dim, h/6, q1, 1, 1, yk, 1);
tk+=h;
//yk





	
int j;
for(j=0;j<dim;j++)


//yk1=yk+h/6*(q1+2*q1+2*q3+q4);



for(j=0;j<dim;j++)
{jieguo[weizhi+j]=yk[j];

	
	
		
		fprintf(logf , "%lf,%lf" , yk[j],tk);
		fprintf(logf,"\n");
		


}



weizhi+=dim;

}
fclose(logf);


return q1;
}






double _Complex lnk(double _Complex a,int k)
{	double _Complex y;
	y=clog(a)+2*k*I*PI;
	return y;
	
	
}


//第一类斯特林数
int st(int p,int k)
{
	if (p==k)
	{return 1;}
	else if (k==0&&p>0)
	{return 0;}
	else 
	{return st(p-1,k-1)+(p-1)*st(p-1,k);}
	
	
	
}



int  jiecheng(int n)
{
    int i;
    int res=1;
    for(i=2; i<=n; ++i)
        res = res * i;
    return res;
}


double _Complex lamw(double _Complex a,int k)
{	double _Complex y;
	y=lnk(a,k)-lnk(lnk(a,k),0);
	for(int l=0;l<10;l++)
	for(int m=1;m<10;m++)
	{y+=cpow(lnk(lnk(a,k),0),m)/cpow(lnk(a,k),l+m)*st(l+m,l+1)/jiecheng(m)*cpow(-1,l);
	

	}
	
	
	
	return y;
	
	
}



void shuchuf(float *p,int a,int b)
{     int i,j;
	for(i=1;i<=a;i++)
	{	for(j=1;j<=b;j++)
		{printf("%f,",p[(i-1)*b+j-1]);}
	printf("\n");
	}
	
}

void shuchud(double *p,int a,int b)
{     int i,j;
	for(i=1;i<=a;i++)
	{	for( j=1;j<=b;j++)
		{printf("%f,",p[(i-1)*b+j-1]);}
	printf("\n");
	}
	
}
