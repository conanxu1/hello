#include <stdio.h>
#include <string.h>

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
		{printf("%.15f,",p[(i-1)*b+j-1]);}
	printf("\n");
	}
	
}





void svd(double* p,int m,int n)
{
	double *Q,*V,*B;
	
	B=(double *)malloc(m*m*sizeof(double));
	
	

	sduijiao(p,m,m,Q,V,B);
	
	qrweiyi(B,n,n);
	
}













void sduijiao(double* p,int m,int n,double* Q2,double* V1,double* B)
{
	//双对角约化
	double *Q;
	double	*PPP,*tem,*tt1,*tt2,*ttt,*Q1;
	
	Q=(double *)malloc(m*n*sizeof(double));
	tem=(double *)malloc(m*n*sizeof(double));
	
	Q1=danwei(m,m);
	V1=danwei(n,n);
	tt1=danwei(m,m);
	tt2=danwei(n,n);
	
	
	memcpy(Q,p,m*n*sizeof(double));
	double h1;
	double *ui,*uj;
	ui=(double *)malloc(m*sizeof(double));
	uj=(double *)malloc(n*sizeof(double));
	PPP=danwei(1,1);

	//双对角化
	for(int i=0;i<m-1;i++)
	{	free(PPP);
		
		
		
		
		xqu(Q,ui,i,m,n);
		house(ui,m,i);
		PPP=danwei(m,m);
		
		
		if(cblas_ddot(m,ui,1,ui,1)>1e-17)
		{cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasTrans, m, m,1, -2/cblas_ddot(m,ui,1,ui,1),ui, 1,ui,1, 1,PPP, m);
		}
		
		
		tem=danwei(m,n);
		cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, m, n,m, 1,PPP, m,Q,n, 0,tem, n);
		free(Q);
		Q=tem;
		
		
		
		
		
		cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, m, m,m, 1,PPP, m,Q1,m, 0,tt1, n);
						ttt=Q1;
						Q1=tt1;
						tt1=ttt;
		
		
		
		xqu2(Q,uj,i,m,n);
		house(uj,n,i+1);
		free(PPP);
		PPP=danwei(m,m);
		if(cblas_ddot(m,ui,1,ui,1)>1e-17)
		{
		cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasTrans, n, n,1, -2/cblas_ddot(n,uj,1,uj,1),uj, 1,uj,1, 1,PPP, n);
		}
		
		tem=danwei(m,n);
		cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, m, n,n, 1,Q, n,PPP,n, 0,tem, n);                                                  
		free(Q);
		Q=tem;
		cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, m, m,m, 1,PPP, m,V1,m, 0,tt2, n);
						ttt=V1;
						V1=tt2;
						tt2=ttt;
		
		
	}
	

	tem=danwei(n,n);
	cblas_dgemm(CblasRowMajor, CblasTrans,CblasNoTrans, n, n,m, 1,Q, n,Q,n, 0,tem, n);  
	
	
	
	
	
	free(Q);
	Q=tem;
	
	
	
	memcpy(B,Q,m*n*sizeof(double));
	
		
	
	
	//Q1,V1  左右
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}



double* danwei(int m,int n)
{
	
	
	
	double* y =(double *)malloc(m*n*sizeof(double));
	
	memset(y,0,m*n*sizeof(double));
	
	
	if(m>n)
	{for(int i=0;i<n;i++)
	y[i*n+i]=1.0;
//(i,i)(i,j)
	}
	else
	{for(int i=0;i<m;i++)
	y[i*n+i]=1.0;
	//(i,i)(i,j)
	}
	return y;



	
	
}



void house(double* a,int dim,int i)
{
	
	
	
	if(a[i]<0)
	a[i]=a[i]-sqrt(cblas_ddot(dim,a,1,a,1));
	else
	a[i]=a[i]+sqrt(cblas_ddot(dim,a,1,a,1));
	
	// cblas_daxpy(dim,1/sqrt(dot(a,a,dim))-1,a,1, a, 1);
	
	
	
	
	
}







double dot(double* a,double* b,int dim)
{
double y[0],d;
cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasTrans, 1, 1,dim, 1,a, dim,b,dim, 0,y, 1); 
d=y[0];
return d;
}



//配合双对角约化
void xqu(double* A,double* u,int i,int m,int n)
{

//从0

memset(u,0,m*sizeof(double));
	

for(int j=i;j<m;j++)
{
	u[j]=A[j*n+i];
		
}
}

void xqu2(double* A,double* u,int i,int m,int n)
{

//从0
memset(u,0,n*sizeof(double));

for(int j=i+1;j<n;j++)
{
	u[j]=A[i*n+j];
		
}
}






void qrweiyi(double* A,int m,int n)
{
	
	double* Q1;
	double* R1;
	double* PPP;
	double* H;
	double* EE;
	
	
	
	EE=(double *)malloc(m*m*sizeof(double));//jiaohuan
	
	
	Q1=(double *)malloc(m*m*sizeof(double));
	R1=(double *)malloc(m*n*sizeof(double));
	
	
	H=(double *)malloc(m*n*sizeof(double));
	
	memcpy(H,A,m*n*sizeof(double));
	
	
	
	PPP=danwei(m,m);
	double h,de;
	
	for(int k=0;k<20;k++)
	{	de=(H[(m-2)*(n+1)]-H[(m-1)*(n+1)])/2;
		if(de>0)
			{h=H[(m-1)*(n+1)]+de-sqrt(de*de+H[(m-1)*(n+1)-1]*H[(m-1)*(n+1)-1]);
			}
		else
			{h=H[(m-1)*(n+1)]+de+sqrt(de*de+H[(m-1)*(n+1)-1]*H[(m-1)*(n+1)-1]);
			}
	
		for(int i=0;i<m;i++)
		{H[i*n+i]-=h;}
			
			
			
		//H-sigI
		
		//局部变量内存
				 
				 
				 
				
			QR(H,m,n,Q1,R1);
				
				
			
			cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, m, m,m, 1,R1, m,Q1,m, 0,H, m);
		
		
		for(int i=0;i<m;i++)
		{H[i*n+i]+=h;}
	
	
	
	}
	
	 
	for(int i=0;i<m;i++)	
	{
		printf("%0.13f,\n",sqrt(fabs(H[i*n+i])));

}
	printf("\n-----\n\n\n ");
				
	
	//abs 整数
	}



void QR(double* A,int m,int n,double* Q1,double* R1)
{
	
	
	
	
	double *EE,*OO,*ui,*tem,*ttt,*QQ,*V1,*Q,*PPP;
	
	Q=(double *)malloc(m*n*sizeof(double));
	ui=(double *)malloc(m*sizeof(double));
	
	QQ=danwei(m,m);
	
	memcpy(Q,A,m*n*sizeof(double));
	EE=danwei(m,m);
//交换用

tem=danwei(m,n);



	
for(int i=0;i<m-1;i++)
{

xqu(Q,ui,i,m,n);
house(ui,m,i);

//ui

PPP=danwei(m,m);





if(cblas_ddot(m,ui,1,ui,1)>1e-17)
{	
cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasTrans, m, m,1, -2/cblas_ddot(m,ui,1,ui,1),ui, 1,ui,1, 1,PPP, m);
}
//P  此处多次申请内存 可优化

 
cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasTrans, m, m,m, 1,QQ, m,PPP,m, 0,EE, m);
 
					ttt=EE;
					EE=QQ;
					QQ=ttt;
 

 
 
 
 
 
 
cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, m, n,m, 1,PPP, m,Q,n, 0,tem, n);
//PA
					ttt=Q;
					Q=tem;
					tem=ttt;



					
free(PPP);
}		
//




free(EE);
free(ui);
free(tem);



memcpy(Q1,QQ,m*m*sizeof(double));
	
memcpy(R1,Q,m*n*sizeof(double));













	
	
		
	
}


double* fuzhi1(double h1,int i,int m)
{
	
	double* P1=(double *)malloc(m*m*sizeof(double));
	memset(P1,0,m*m*sizeof(double));
	
	for(int j=0;j<m;j++)
	{P1[j*m+j]=1;}



	for(int j=i;j<m;j++)
	P1[j*m+j]+=1+h1;

	return P1;
	
	
}



double* duqu(char *p,int n)
{
	
	FILE *fp;
	double* xx;
	xx=(double *)malloc(n*sizeof(double));
	
	fp=fopen(p,"r");
	
	if(fp!=NULL)
	{
		for(int j=0;j<n;j++)
		{fscanf(fp,"%lf\n",&xx[j]);
		
		}
		
	}
	fclose(fp);
	return xx;
}












	// OO=(double *)malloc(n*n*sizeof(double));
	// memset(OO,0,n*n*sizeof(double));



