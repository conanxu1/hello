#include <stdio.h>
#include <string.h>
#include "fftw3.h"


#include<stdarg.h>




#include <stdlib.h>
#include "myfun.h"
#include <complex.h>

#include <math.h>

#define PI acos(-1.0)
// void pq()
// {
// double B[4*2]={1,0,2,0,3,0,4,0};
// double B1[4*2]={1,0,2,0,3,0,4,0};
// double *B2;
// int s=4;
// double A[4*4*2]={1,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,-4,0};
// //-A -A'







// double *Al,*Ar,*Alc,*Arc,*tem,*ttt;
// tem=(double *)malloc(16*2*sizeof(double));
// Al=(double *)malloc(16*2*sizeof(double));
// Ar=(double *)malloc(16*2*sizeof(double));

// Alc=(double *)malloc(16*sizeof(double));
// Arc=(double *)malloc(16*sizeof(double));





// int ipiv[s];
// int info;

// double T=50000;
// int n=1000000;
// double h=2*T/n;


// B2=(double *)malloc(16*2*sizeof(double)); 

// double yi[2]={1.0,0};
// double ling[2]={0,0};
// double er[2]={2,0};

// double shh[2]={h/2/2/3.141592653589793,0};

// cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasTrans, s, s,1,yi,B, 1,B1,1, ling,B2, s);









// double *Sum;
// Sum=(double *)malloc(s*s*2*sizeof(double));

// for(int j = 1; j < s*s*2; j++)
// Sum[j]=0;


// ////////////////////////////////////////



// double tt;
// for(int j = 1; j < n; j++)
// {
// tt=-T+j*h;					
// cblas_zaxpby(s*s, yi, A, 1, ling, Al, 1);
// cblas_zaxpby(s*s, yi, A, 1, ling, Ar, 1);

// for(int oo=0;oo<s;oo++)
// {
// Al[(oo*s+oo)*2+1]+=tt;
// Ar[(oo*s+oo)*2+1]+=-tt;
// }








// info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR,s,s,Al,s,ipiv);
// info = LAPACKE_zgetri(LAPACK_ROW_MAJOR,s,Al,s,ipiv);



	
// info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR,s,s,Ar,s,ipiv);
// info = LAPACKE_zgetri(LAPACK_ROW_MAJOR,s,Ar,s,ipiv);



	
// cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, s, s,s,yi,Al, s,B2,s, ling,tem, s);

	
// //////////
	
// cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, s, s,s,yi,tem, s,Ar,s, ling,Al, s);




	

// cblas_zaxpby(s*s, er, Al, 1, yi, Sum, 1);




// }
// //////////////////////////////////////////
// tt=-T;					
// cblas_zaxpby(s*s, yi, A, 1, ling, Al, 1);
// cblas_zaxpby(s*s, yi, A, 1, ling, Ar, 1);

// for(int oo=0;oo<s;oo++)
// {
// Al[(oo*s+oo)*2+1]+=tt;
// Ar[(oo*s+oo)*2+1]+=-tt;
// }






// info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR,s,s,Al,s,ipiv);
// info = LAPACKE_zgetri(LAPACK_ROW_MAJOR,s,Al,s,ipiv);

	
// info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR,s,s,Ar,s,ipiv);
// info = LAPACKE_zgetri(LAPACK_ROW_MAJOR,s,Ar,s,ipiv);

	
// cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, s, s,s,yi,Al, s,B2,s, ling,tem, s);



	
// cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, s, s,s,yi,tem, s,Ar,s, ling,Al, s);

// cblas_zaxpby(s*s, yi, Al, 1, yi, Sum, 1);


// ////////////////////////////////////////
// tt=T;					
// cblas_zaxpby(s*s, yi, A, 1, ling, Al, 1);
// cblas_zaxpby(s*s, yi, A, 1, ling, Ar, 1);

// for(int oo=0;oo<s;oo++)
// {
// Al[(oo*s+oo)*2+1]+=tt;
// Ar[(oo*s+oo)*2+1]+=-tt;
// }
// info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR,s,s,Al,s,ipiv);
// info = LAPACKE_zgetri(LAPACK_ROW_MAJOR,s,Al,s,ipiv);

	
// info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR,s,s,Ar,s,ipiv);
// info = LAPACKE_zgetri(LAPACK_ROW_MAJOR,s,Ar,s,ipiv);

	
// cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, s, s,s,yi,Al, s,B2,s, ling,tem, s);

	

	
// cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, s, s,s,yi,tem, s,Ar,s, ling,Al, s);



// cblas_zaxpby(s*s, yi, Al, 1, yi, Sum, 1);







// cblas_zaxpby(s*s, shh, Sum, 1, ling, tem, 1);







// printf("the matrix ALLLLLLLLLL is:\n");
	// for (int i = 0; i<s; i++)
	// {
		// for(int j=0;j<s;j++)
		// {
			// printf("%f+%f*I\t",tem[(i*s+j)*2],tem[(i*s+j)*2+1]);
		// }
		// printf("\n");
	// }



// // tt=-T
// // cblas_zaxpby(s*s, 1, Al, 1, s*s, Sum, 1);
// // tt=T

// // h/2*
	

	
// }













double jifen(double a,double b,int n,double (*fff)(double))
{	double h = (b - a)/n, T = 0; 
	int i;
	for(i = 1; i < n; i++) 
		T += (*fff)(a + i * h); 
	return  h/2*(fff(a)+T*2+fff(b));
}



//改写为带控制�?


double  *ode1(double a,double b,int n,double *(*fff)(double ,double *),double *x0,int dim,double *jieguo)
{/*a初�?b终�?n分点,fff连续函数*/   

FILE *logf;
//记录�?


double h = (b - a)/n;
//划分

double t0=a,tk;



double* yk =(double *)malloc(dim*sizeof(double));
double* q1 =(double *)malloc(dim*sizeof(double));
double* q2 =(double *)malloc(dim*sizeof(double));
double* q3 =(double *)malloc(dim*sizeof(double));
double* q4 =(double *)malloc(dim*sizeof(double));

double* tem=(double *)malloc(dim*sizeof(double));


//*a要分配内�?或者定义时候分�?
//double *dy =(double *)malloc(size*sizeof(double));

int i;
int j;



logf=NULL;
if((logf = fopen("log.txt" , "w+")) == NULL)
{printf("Cannot create/open file");
exit(1);}
//打开文件



tk=t0;
yk=x0;


for(j=0;j<dim;j++)
{jieguo[j]=yk[j];
}
jieguo[dim]=t0;
int weizhi=dim+1;




for(i = 0; i < n; i++)
{
q1=fff(tk,yk);


cblas_daxpby(dim, 1, yk, 1, 0, tem, 1);
cblas_daxpby(dim,h/2, q1, 1, 1, tem, 1);
q2=fff(tk+h/2,tem);
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








//yk1=yk+h/6*(q1+2*q1+2*q3+q4);



for(j=0;j<dim;j++)
{jieguo[weizhi+j]=yk[j];
fprintf(logf , "%lf," , jieguo[weizhi+j]);
}




jieguo[weizhi+dim]=tk;
fprintf(logf , "??%lf" , tk);
fprintf(logf,"\n");


weizhi=weizhi+dim+1;

}
fclose(logf);

shuchud(jieguo,n+1,dim+1);



return q1;
}

























double _Complex lnk(double _Complex a,int k)
{	double _Complex y;
	y=clog(a)+2*k*I*PI;
	return y;
	
	
}


//第一类斯特林�?
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
{	

//用渐进公�? 截断多项式给初�? 通过牛顿迭代 梯度下降   求精确解





double _Complex y0,y;




if(cabs(a)<1/exp(1))
{y0=a-cpow(a,2)+3/2*cpow(a,3)-8/3*cpow(a,4)+125/24*cpow(a,5);
}
else if(cabs(a)>1/exp(1)&&cabs(a)<10)
{	
y0=lnk(a,k)*(0.05+0*I);
// printf("\n99\n");
}
else
{
	// printf("\n\n\n3\n\n\n");

	y0=lnk(a,k)-lnk(lnk(a,k),0);
	for(int l=0;l<10;l++)
	for(int m=1;m<10;m++)
	{y0+=cpow(lnk(lnk(a,k),0),m)/cpow(lnk(a,k),l+m)*st(l+m,l+1)/jiecheng(m)*cpow(-1,l);}
}

int pp;
double w1,w2,z1,z2,w1k,w2k,v1,v2,norm,test;

z1=creal(a);
z2=cimag(a);

w1=creal(y0);
w2=cimag(y0);

 //printf("预估�?.6f,%.5f\n\n\n",w1,w2);



y=y0;



while(cabs(a-y*cexp(y))>1e-12)
{	
w1=creal(y);
w2=cimag(y);
// printf("当前�?lf,%lf\n",w1,w2);




v1=(2*w1*exp(2*w1)+2*(w1*w1+w2*w2)*exp(2*w1)-2*exp(w1)*(w1*z1*cos(w2)+w2*z2*cos(w2)-w2*z1*sin(w2)+w1*z2*sin(w2)+z1*cos(w2)+z2*sin(w2)));
	
	
v2=(2*w2*exp(2*w1)-2*exp(w1)*(-w1*z1*sin(w2)+z2*cos(w2)-w2*z2*sin(w2)-z1*sin(w2)-w2*z1*cos(w2)+w1*z2*cos(w2)));


norm=sqrt(v1*v1+v2*v2);
v1=v1/norm;
v2=v2/norm;


//线性搜�?
w1k=w1-1e-15*v1;	
w2k=w2-1e-15*v2;	
	
y=w1k+w2k*I;
	
test=cabs(a-y*cexp(y));
pp=15;


// printf("&&oooo梯度----|%lf,%lf\n",v1,v2);

// printf("误差%lf\n",cabs(a-y*cexp(y)));



while(pp>-1)
{
	
	--pp;	

w1k=w1-cpow(0.1,pp)*v1;	
w2k=w2-cpow(0.1,pp)*v2;	
	
y=w1k+w2k*I;
// printf("\npp%d\n",pp);





// printf("---------%.15f,%.15f\n",w1k,w2k);

if(cabs(a-y*cexp(y))>=test)
{	pp++;
	w1k=w1-cpow(0.1,pp)*v1;	
	w2k=w2-cpow(0.1,pp)*v2;	
	// printf("000");
	y=w1k+w2k*I;
	
	
	
	break;}
	

	test=cabs(a-y*cexp(y));
	// printf("误差%lf\n",test);
	
}	
// printf("%.13f,%.13f\n",100*(w1-w1k),100*(w2-w2k));
	
	
	
	



// printf("%.6f,%.5f\n",w1,w2);

	
}






	// printf("w%lf,%lf",creal(y),cimag(y));
	
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
{   

int i,j;
	for(i=1;i<=a;i++)
	{	for( j=1;j<=b;j++)
		{printf("%lf,",p[(i-1)*b+j-1]);}
	printf("\n");
	}
	
}

void shuchui(int *p,int a,int b)
{   

int i,j;
	for(i=1;i<=a;i++)
	{	for( j=1;j<=b;j++)
		{printf("%d,",p[(i-1)*b+j-1]);}
	printf("\n");
	}
	
}


void shuchuz(double _Complex *p,int a,int b)
{printf("\n");
for(int i=0;i<a;i++)
{for(int j=0;j<b;j++)
{printf("%f+%fI\t,",creal(p[i*b+j]),cimag(p[i*b+j]));
}
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
	//双对角约�?
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



//配合双对角约�?
void xqu(double* A,double* u,int i,int m,int n)
{

//�?

memset(u,0,m*sizeof(double));
	

for(int j=i;j<m;j++)
{
	u[j]=A[j*n+i];
		
}
}

void xqu2(double* A,double* u,int i,int m,int n)
{

//�?
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
		
		//局部变量内�?
				 
				 
				 
				
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
//交换�?

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
//P  此处多次申请内存 可优�?

 
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

int xxwg(double *AM,double *A,int m,int n)
{
//只求行满�?


double *TEM=(double *)malloc(m*n*sizeof(double));
double *TEM2=(double *)malloc(m*n*sizeof(double));


	memset(TEM,0,m*n*sizeof(double));
	memcpy(TEM,A,m*n*sizeof(double));
	int matrix_order = LAPACK_ROW_MAJOR;
	char jobu = 'A';
	char jobvt = 'A';
	
	int ff=m,lda=n;


	

		
	double *u=(double *)malloc(m*m*sizeof(double));
	
	int ldu = m;
	double *vt=(double *)malloc(n*n*sizeof(double));
	int ldvt = n;
	
	
	int ot=(m+n)/2+abs(m-n)/2;
	double *superb=(double *)malloc((m+n)*sizeof(double));
	double *s=(double *)malloc((m+n)*sizeof(double));

	
	//printf("ok!!!!!\n");
	//printf("ok!!!!!\n");
 //shuchud(TEM,m,n);


	LAPACKE_dgesvd(matrix_order,jobu, jobvt, m, n, TEM,lda, s, u, ldu, vt, ldvt, superb);
 	
	//printf("uuuuuuuu\n\n");
	//shuchud(u,m,m);
	 	
	//printf("vvvvvvvvv\n\n");
	//shuchud(vt,n,n);
	
	
	
	int wz=ff-1;
	if(s[wz]>1e-15)
			{return ff;   } //hangmanzhi


memset(TEM,0,m*n*sizeof(double));

memcpy(TEM,A,1*n*sizeof(double));

int zhi=1;

for(int pp=1;pp<m;pp++)
   {	
	memcpy(TEM2,TEM,zhi*n*sizeof(double));

		
	for(int jj=0;jj<n;jj++)	
		TEM2[zhi*n+jj]=A[pp*n+jj];
	
	//printf("pppp\n");
	//shuchud(TEM2,zhi+1,n);
	
	//printf("\n\n\n00000pppp\n");
	
	
	
	LAPACKE_dgesvd(matrix_order,'A', 'A', zhi+1, n, TEM2,n, s, u, ldu, vt, ldvt, superb);
	shuchud(s,zhi+1,1);
	
	
	
	//printf("\n\n\n0\n\n\n");
 	if(s[zhi]>1e-15)
		{
						for(int jj=0;jj<n;jj++)	
						TEM[zhi*n+jj]=A[pp*n+jj];
					zhi+=1;
		}

	}

//printf("zhi\n");



AM=(double *)malloc(zhi*n*sizeof(double));

memcpy(AM,TEM,zhi*n*sizeof(double));

free(TEM);

free(TEM2);

return zhi;
   

   

}












void svd2()
{
	
	
	



    int matrix_order = LAPACK_ROW_MAJOR;//修改过排列方式了
    char jobu = 'A';
    char jobvt = 'A';
	
	double *A;

	int n=1000;
	int m=n;
	A=duqu("wind.txt",n*n);
	int lwork;
	double wk[6*n];
	double s[n];
    int lda = n;
    int info;
	
	
    double u[n*n];
    int ldu = n;
    double vt[n*n];
    int ldvt = n;
	double superb[n];
//    int info = LAPACKE_dgeev(matrix_order,jobvl,jobvr,n,A,lda,wr,wi,vl,ldvl,vr,ldvr);
	LAPACKE_dgesvd(matrix_order,jobu, jobvt, m, n, A,n, s, u, ldu, vt, ldvt, superb);



	// for(int i=0;i<n;i++)
	// {printf("%.14g\n",s[i]);}

/*   if(info==0){
        int i = 0;
        int j = 0;
        for(i=0;i<n;i++){
            printf("eigenvalue %d:\n",i);
            printf("%.6g + i %.6g \t",wr[i],wi[i]);
            printf("right eigenvector: ");
            for(j=0;j<ldvr;j++)
                printf("%.6g \t",vr[i*4+j]);
            printf("\n");
        }
        printf("SUCCESS\n");

	}
	*/
}


//破坏性逆，也可以先复制

void ni(double *A,int dim)
{
int s=dim;
int ipiv[s];
int info;
info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,s,s,A,s,ipiv);



info = LAPACKE_dgetri(LAPACK_ROW_MAJOR,s,A,s,ipiv);
	
}	
	
	














void mypinv()
{
int matrix_order = LAPACK_COL_MAJOR;
char jobu = 'A';
char jobvt = 'A';
double *A;
int n=1;
int m=n;
A[0]=1;

double s[n];
int lda = n;
int info;
double u[n*n];
int ldu = n;
double vt[n*n];
int ldvt = n;
double superb[n];
LAPACKE_dgesvd(matrix_order,jobu, jobvt, m, n, A,n, s, u, ldu, vt, ldvt, superb);

shuchud(A,n,n);
}










void drk()
{FILE *logf;
logf=NULL;
if((logf = fopen("log.txt" , "w+")) == NULL)
{
	printf("Cannot create/open file");
	exit(1);		
}



	int r=2,s=2,m=40,n=60000;
	//r个y   s个K Y m-tau h分割   n步迭�?
	
	int dim=2;
	//维数
	double tau=1.1;
	
	
	double L[4]={-2.0,0,0,-0.9};
	double M[4]={-1,0,-1,-1};
	double N[4]={0.9,0.45,0,0.05};
	
	// double *A,*B,*the,*gam;
	double A[4] = {5.0/14,9.0/14,-1.0/2,3.0/2};
	
	double B[4]={15.0/14,-5.0/14,-1.0/2,3.0/2};
	double gam[2]={21.0/20,3.0/20};
	double the[2]={1.0/10,9.0/10};
	
	double ***Yni,***Kni,**yn;
	
	double h=tau/m;
	
	
	double *C;
	
	C=danwei(s,s);	
	
	
	for(int i=0;i<s;i++)
	{
		C[i*s+i]-=B[i*s+i]*h;
	}
	
	for(int i=0;i<s;i++)
	for(int j=0;j<s;j++)
	{
		if(i!=j)
		C[i*s+j]=-B[j*s+i]*h;
	}
	
	shuchud(C,s,s);
	
	
	int ipiv[s];
	int info;
	
	
	
 info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,s,s,C,s,ipiv);
info = LAPACKE_dgetri(LAPACK_ROW_MAJOR,s,C,s,ipiv);
	shuchud(C,s,s);
	
	
	
	
	double Cs[s];
	
	
	for(int i=0;i<s;i++)
	{	Cs[i]=0;
		
		for(int q=0;q<s;q++)
		Cs[i]+=C[q*s+i];
	}
	
	
	
	
	
	
	// A=(double*)malloc(sizeof(double)*s*r);  
    // B=(double*)malloc(sizeof(double)*s*r);  

	
	Yni=(double***)malloc(sizeof(double**)*s);  
    for(int i=0;i<s;i++)  
    {Yni[i]=(double**)malloc(sizeof(double*)*(m+1));  
		for(int j=0;j<m+1;j++)  
		Yni[i][j]=(double*)malloc(sizeof(double)*dim); 
	}
/****************************/

	Kni=(double***)malloc(sizeof(double**)*s);  
    for(int i=0;i<s;i++)  
    {Kni[i]=(double**)malloc(sizeof(double*)*(m+1));  
		for(int j=0;j<m+1;j++)  
		Kni[i][j]=(double*)malloc(sizeof(double)*dim); 
	}  
// /****************************/


	


	yn=(double**)malloc(sizeof(double*)*(r+1));  
    // for(i=0;i<s;i++)  
    // yn[i]=(double*)malloc(sizeof(double)*dim);  



/****  求A   B  theta gamma  ***/
	
	
	
	
/*******    初始�?  **Yn,**Kn,**yn,**Kn_m,**Yn_m   **********/
//初始化是相同
	
	double *tem,*tem2;
	tem=(double*)malloc(sizeof(double)*dim);  
    int wei;
	
	
	
	
inrk(Yni,Kni,yn,r,s,m,tau);

// for(int oo=0;oo<s;oo++)
// {for(int uu=0;uu<m+1;uu++)
// {shuchud(Kni[oo][uu],dim,1);
// printf("====\n");
// }
// }	
	
	

for(int tt=0;tt<n;tt++)
{
	
	for(int i=0;i<s;i++)
	{
		
		
		
		
	
//减少整体平移故作位置变换在原数组上覆盖数�?
			
	wei=(tt)%(m+1);
	cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, dim, 1,dim, 1,M, dim,Yni[i][wei],1, 0,tem,1 );
	
	wei=(tt)%(m+1);
	cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, dim, 1,dim, 1,N, dim,Kni[i][wei],1, 1,tem,1 );
	
	//n-m最左边
	
	
			for(int j=0;j<r;j++)
		{	
	
		wei=(tt+j)%(r+1);
			
		
		
		cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, dim, 1,dim, A[i*(r)+j],L, dim,yn[wei],1, 1,tem,1 );
	
			
			
		}
	
	////////////////////beta C�?
	//【Kn1,Kn2,....,Kns】C=[beta,beta,...,beta]
	//【Kn1,Kn2,....,Kns�?[beta,beta,...,beta]*(C^-1)
	//   Kni_oo=    beta(oo) sum_j Cji
	
	
	wei=(tt)%(m+1);
	

	cblas_daxpby( dim,Cs[i],tem, 1,0,Kni[i][wei],1);
	
	
	
	}
	//交换地址Kni
	
	
	
	
	for(int i=0;i<s;i++)
	{
	
	
					for(int jj=0;jj<dim;jj++)
					tem[jj]=0;
	
	
	///////////////////////////////////////////////////
	//区分前后
		for(int j=0;j<r;j++)
		{	
	
		wei=(tt+j)%(r+1);
		
		
		
		cblas_daxpby( dim,A[i*(r)+j],yn[wei], 1,1,tem,1);
		
		
		// printf("---%d\n\n",wei);
		// shuchud(tem,dim,1);
		
		
		}
	//yn特殊最左边为当前时�?
	
		for(int j=0;j<s;j++)
		{	
	
	
	//Kni已经更新�?
		wei=(tt)%(m+1);
		cblas_daxpby( dim,h*B[i*(s)+j],Kni[j][wei], 1,1,tem,1);
	
		}	
	//////////////////////////////////////////////////////
	wei=(tt)%(m+1);
	tem2=tem;
	tem=Yni[i][wei];
	Yni[i][wei]=tem2;
	//n+r+1,1,2,......
	
	
	}
	
	
	
	
	for(int jj=0;jj<dim;jj++)
	tem[jj]=0;
	//////////////////////////////////////////////////
	
	
			for(int j=0;j<r;j++)
		{	
	
		wei=(tt+j)%(r+1);
		cblas_daxpby( dim,the[j],yn[wei],1, 1,tem,1);
	
		}
	
	
			for(int j=0;j<s;j++)
		{	
	//已更�?
		wei=(tt)%(m+1);
		cblas_daxpby( dim,h*gam[j],Kni[j][wei],1, 1,tem,1);
	
		}
	
	
	
	
	
wei=(tt+r)%(r+1);
tem2=tem;
tem=yn[wei];
yn[wei]=tem2;
	
	
fprintf(logf , "%lf,%lf" , yn[wei][1],tt*h);
fprintf(logf,"\n");
	
	
	
	// printf("\n");
	// shuchud(tem2,dim,1);
	
	
	// printf("\n");
	
	
	
}
	

		
		// Kn[i]=L*Yn[i]+M*Yn_m[i]+N*Kn_m[i];
		  
		// Yn[i]=L*Yn[i]+M*Yn_m[i]+N*Kn_m[i];
		
		
		
		
	// }
	
	
	
	
	
	
	
	fclose(logf);
	
	
	
}













void inrk(double ***Yni,double ***Kni,double **yn,int r,int s,int m,double tau)
{double h=tau/m;
	
	for(int i=0;i<s;i++)  
    { 
		for(int j=0;j<m+1;j++)  
		{
			
		Yni[i][j]=fini((j-m)*h); 
	
		Kni[i][j]=dfini((j-m)*h); 
		
		}
	} 
	
	
	
	
	
	for(int j=0;j<r+1;j++)
	yn[j]=fini(0);
	
	//默认不变
	
	
	
	
	
	
	
	
}





//考虑特征值互不相同的情况

void jiangjie()
{
	

	printf("\n");
	
	
	
	
}

//为计算方�?特征值互不相�?
int lamwm(double _Complex* A,int k,int dim,double _Complex*	 w)
{
double _Complex u[dim*dim];
double _Complex v[dim*dim],tem1[dim*dim];

	char jobvl='V';
	char jobvr='V';
	int n=dim;
	int lda=dim;
	int ldvl=dim;
	int ldvr=dim;
	double _Complex	 tem2[n];
	//返回特征值向�?
	
	
	 
	int matrix_order = LAPACK_ROW_MAJOR;
	
	
	LAPACKE_zgeev(matrix_order ,jobvl, jobvr, n, A, lda,w, u , ldvl , v, ldvr);
	
	// printf("--\n");
	 // //shuchuz(w,1,dim);//Av=vD
	// printf("--\n");
	
	
int ipiv[dim];
int info;

	
//左右特征值矩阵的值应该是相同�?至少可逆性一�?考察右特征值是否可�?

LAPACKE_zgeev(matrix_order ,jobvl, jobvr, n, u, lda,tem2, tem1 , ldvl , A, ldvr);
	
	

		
		
memcpy(A,v,dim*dim*sizeof(double _Complex));		
		

	for(int i=0;i<dim;i++)
	{
		if(cabs(tem2[i])<1e-14)
			return 1;
	
		
			
	//不可逆失�?
	
	
	
	
	
	}	
	
	 // printf("+++++++++++++++++++\n");
	
	
	// shuchuz(w,1,dim);
	
	for(int i=0;i<dim;i++)
	{
		// printf("\nD\n%lf,%lf\n",creal(w[i]),cimag(w[i]));
		w[i]=lamw(w[i],k);
		// printf("......\n");
	
	}
	return 0;
	
}







void lamwhb(double _Complex* A,double _Complex* w,int dim)
{
double _Complex A2[dim*dim],jjj[dim*dim];
int ipiv[dim];
int info;

double  yi[2]={1,0},ling[2]={0,0};




memcpy(A2,A,dim*dim*sizeof(double _Complex));		

//右列
for(int i=0;i<dim;i++)	
for(int j=0;j<dim;j++)
{
	A[i*dim+j]*=w[j];
}

info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR,dim,dim,A2,dim,ipiv);
info = LAPACKE_zgetri(LAPACK_ROW_MAJOR,dim,A2,dim,ipiv);


cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, dim, dim,dim,yi,A, dim,A2,dim, ling,jjj, dim);

	
	
	
memcpy(A,jjj,dim*dim*sizeof(double _Complex));		


	
	
	
}

















double mubiao(double _Complex* A,double _Complex* Ad,double _Complex* Q,double h,int dim)
{
	double _Complex tem1[dim*dim],tem2[dim*dim],tem3[dim*dim],w[dim];
	
	double xishuh[2]={h,0},ling[2]={0,0},yi[2]={1,0},fuyi[2]={-1,0},xifuh[2]={-h,0};
	
	cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, dim, dim,dim,xishuh,Ad, dim,Q,dim, ling,tem1, dim);

	// printf("W1\n");
	
	
	
	
	
	int flag;
	
	flag=lamwm(tem1,0,dim,w);
	
	// printf("*****\n");
	
	
	// printf("%d\n\n",flag);
	
	
		
	
	// shuchuz(w,1,dim);
	
	lamwhb( tem1,w,dim);
	// shuchuz(tem1,dim,dim);	

	
	
	
	
	memcpy(tem2,tem1,dim*dim*sizeof(double _Complex));		

	cblas_zaxpby(dim*dim,xishuh,A,1,yi,tem2,1 );
	//G+Ah
	
	expm(tem2,dim);
	
	cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, dim, dim,dim,yi,tem1, dim,tem2,dim, ling,tem3, dim);

	
	
	memcpy(tem2,Ad,dim*dim*sizeof(double _Complex));
	cblas_zaxpby(dim*dim,yi,tem3,1,xifuh,tem2,1 );
	//-adh
	
	// printf("////\n");
	// shuchuz(tem2,dim,dim);
	
	return cblas_dznrm2(dim*dim,tem2,1);
	
}






void qiuQ(double _Complex* A,double _Complex* Ad,double _Complex* Q,double h,int dim)
{double _Complex gra[dim*dim],tem[dim*dim];
double xishu[2]={1,0},yi[2]={1,0};
double wucha,mo,bu,test;
int pp;
wucha=mubiao(A,Ad,Q,h,dim);
while(wucha>1e-14)	
{	FILE *fp;
	fp=fopen("Q.txt","w");
	if(fp!=NULL)
	{
		for(int j=0;j<dim*dim;j++)
		{fprintf(fp,"%lf+%lf*I\n",creal(Q[j]),cimag(Q[j]));}
	}
fclose(fp);




wucha=mubiao(A,Ad,Q,h,dim);
	 shuchuz(Q,dim,dim);
	
	
	nugrad(A,Ad,Q,h,dim, gra);
	memcpy(tem,Q,dim*dim*sizeof(double _Complex));		
	
	mo=cblas_dznrm2(dim*dim,gra,1);
	pp=-3;
	
	xishu[0]=-pow(10,pp);
	cblas_zaxpby(dim*dim,xishu,gra,1,yi,tem,1 );
	
	printf("))))))))");
	test=mubiao(A,Ad,tem,h,dim);
	for(;pp<3;)
	{pp++;

	memcpy(tem,Q,dim*dim*sizeof(double _Complex));		
	

	xishu[0]=-pow(10,pp);
	cblas_zaxpby(dim*dim,xishu,gra,1,yi,tem,1 );
	printf("pp%d\n",pp);
	printf("gra\n");
	shuchuz(gra,dim,dim);
	
	
	mubiao(A,Ad,tem,h,dim);
	printf("ooo\n");

			if(test>mubiao(A,Ad,tem,h,dim))
			{	pp--;
				printf("oo222o\n");
	
	
				break;
			}
			
			test=mubiao(A,Ad,tem,h,dim);
			printf("oo33366So\n");
			
	}
	printf("oo333So\n");
	
	
	xishu[0]=-pow(10,pp)/mo;
	cblas_zaxpby(dim*dim,xishu,gra,1,yi,Q,1 );
	
printf("误差%lf\n",wucha);
	
	
	
}	
	
	
	
	
}








void nugrad(double _Complex* A,double _Complex* Ad,double _Complex* Q,double h,int dim,double _Complex* gra)
{

double _Complex Qh[dim*dim];



double op=1,tem;

for(int i=0;i<dim*dim;i++)
{
	memcpy(Qh,Q,dim*dim*sizeof(double _Complex));		
	
	Qh[i]+=op;
	tem=(mubiao(A,Ad,Qh,h,dim)-mubiao(A,Ad,Q,h,dim))/h;
	
	Qh[i]+=-op+op*I;
	
	gra[i]=(mubiao(A,Ad,Qh,h,dim)-mubiao(A,Ad,Q,h,dim))/h*I+tem;
	
	
	
	
	
	// printf("xxxxxxxxxxxxxxxxxxx%lf\n",cimag(gra[i]));
	
}






}


int expm(double _Complex* A,int dim)
{
	//shuchuz(A,dim,dim);//Av=vD
	
	
	
double _Complex u[dim*dim];
double _Complex v[dim*dim],tem1[dim*dim],w[dim];

	char jobvl='V';
	char jobvr='V';
	int n=dim;
	int lda=dim;
	int ldvl=dim;
	int ldvr=dim;
	double _Complex	 tem2[n];
	//返回特征值向�?
	
	
	 
	int matrix_order = LAPACK_ROW_MAJOR;
	
	
	LAPACKE_zgeev(matrix_order ,jobvl, jobvr, n, A, lda,w, u , ldvl , v, ldvr);
	
	
int ipiv[dim];
int info;

	
//左右特征值矩阵的值应该是相同�?至少可逆性一�?考察右特征值是否可�?

LAPACKE_zgeev(matrix_order ,jobvl, jobvr, n, u, lda,tem2, tem1 , ldvl , A, ldvr);
	
	

		
		
memcpy(A,v,dim*dim*sizeof(double _Complex));		
		

	for(int i=0;i<dim;i++)
	{
		if(cabs(tem2[i])<1e-14)
				return 1;
	//不可逆失�?
	}	
	
	
	
	for(int i=0;i<dim;i++)
	{
		w[i]=cexp(w[i]);
	}
	
	
	
	
double _Complex A2[dim*dim],jjj[dim*dim];


double  yi[2]={1,0},ling[2]={0,0};




memcpy(A2,A,dim*dim*sizeof(double _Complex));		

//右列
for(int i=0;i<dim;i++)	
for(int j=0;j<dim;j++)
{
	A[i*dim+j]*=w[j];
}

info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR,dim,dim,A2,dim,ipiv);
info = LAPACKE_zgetri(LAPACK_ROW_MAJOR,dim,A2,dim,ipiv);


cblas_zgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, dim, dim,dim,yi,A, dim,A2,dim, ling,jjj, dim);

	
	
	
memcpy(A,jjj,dim*dim*sizeof(double _Complex));		


	
	
	//shuchuz(A,dim,dim);
	
	
	
	
	
	return 0;
	
}

/*******傅里叶延时辨�?*/
//傅里叶基函数列向�?N+1
int phi(double* f,int N,double t,double T)
{
	
	printf("%30.30f",PI);
	
		f[0]=1;
	for(int i=1;i<=N;i++)
		{f[i]=cos(2*PI*t*i/T);}
	for(int i=1;i<=N;i++)
		{f[i+N]=sin(2*PI*t*i/T);}


return 1;


//tem=(double *)malloc(16*2*sizeof(double))

}

//D1

int initxishu(double* D1,double* D0,double* D,double* S,int N,double T,double tau)
{
D1=(double *)malloc((2*N+1)*(2*N+1)*sizeof(double));
D0=(double *)malloc((2*N+1)*(2*N+1)*sizeof(double));
S=(double *)malloc((2*N+1)*(2*N+1)*sizeof(double));

memset(D1,0,(2*N+1)*(2*N+1)*sizeof(double));
memset(D0,0,(2*N+1)*(2*N+1)*sizeof(double));
memset(S,0,(2*N+1)*(2*N+1)*sizeof(double));

int jie=(2*N+1);

//D1
for(int i=1;i<=N;i++)	
{D1[(i-1+N)*jie+(i+N)]=-2*PI/T*i;

D1[(i-1+2*N)*jie+(i)]=2*PI/T*i;
}

//D0
D0[0]=0.5;
for(int i=1;i<=N;i++)	
{D1[(i)*jie+(i)]=1;
}


S[0]=1;
for(int i=1;i<=N;i++)	
{
S[(i)*jie+(i)]=cos(2*PI*i*tau/T);
S[(i)*jie+(i+N)]=sin(2*PI*i*tau/T);
S[(i+N)*jie+(i)]=-sin(2*PI*i*tau/T);
S[(i+N)*jie+(i+N)]=cos(2*PI*i*tau/T);
}

return 0;	
}





/*
int m0618()
{double ep=1e-8;
	double tau=(sqrt(5)-1)/2;
	double a=-10;
	double b=10;
	double alpha=0.5*(b-a);
	double al=a;
	double ar=b;

	while(b-a>ep)
	{
		al=a+(1-tau)*(b-a);
		ar=a+tau*(b-a);
		alpha=0.5*(a+b);
		if(f(al)<f(ar))
		{
			b=ar;
			ar=al;
			al=a+(1-tau)*(b-a);
		}
		else
		{
			a=al;
			al=ar;
			ar=a+tau*(b-a);
		}
		printf("计算中a=%lf,b=%lf\n",a,b);
	}
	alpha=0.5*(a+b);
	printf("结果%lf\n",alpha);
	return 0;
}



*/












////////////////////////////////////

double shixing(double t,double a,double b)
{
	if(t>=a&&t<b)
		return 1;
	else
		return 0;
	
	
}


// double armijo()


// int mysqp()
// {
	
	
	
	
	
	
// erci()

//追赶�?
int zhui(double *A,double *d,int n,double *jie)
{
//不覆�?


double *ci,*di;
ci=(double *)malloc((n-1)*sizeof(double));
di=(double *)malloc((n-1)*sizeof(double));




// cp=A[(p-1)*n+p];
// bp=A[(p-1)*n+(p-1)];
// ap=A[(p-1)*n+(p-2)];
//cp*=ci[p-1];


int i,p;
ci[0]=A[1]/A[0];
di[0]=d[0]/A[0];


for(p=2;p<=n-1;p++)
{
ci[p-1]=A[(p-1)*n+p]/(A[(p-1)*n+(p-1)]-A[(p-1)*n+(p-2)]*ci[p-2]);	
	
di[p-1]=(d[p-1]-A[(p-1)*n+(p-2)]*di[p-2])/(A[(p-1)*n+(p-1)]-A[(p-1)*n+(p-2)]*ci[p-2]);	
}
p=n;
di[p-1]=(d[p-1]-A[(p-1)*n+(p-2)]*di[p-2])/(A[(p-1)*n+(p-1)]-A[(p-1)*n+(p-2)]*ci[p-2]);


jie[n-1]=di[n-1];
for(p=n-2;p>=0;p--)
{  jie[p]=di[p]-ci[p]*jie[p+1];
}


free(ci);
free(di);

ci=NULL;
di=NULL;


shuchud(jie,n,1);


return 1;
}
	
	


int zuoyongyu(int **a)
{
	
/*
主函数中
int *a;
zuoyongyu(&a);
printf("\n%d\n",*a);

在子函数里用malloc给参数变量分配空间，变量赋值后，主函数的值不会变�?
原因：malloc出来的地址跟main中声明的变量的地址是不一样的，子函数中的赋值语句只是给malloc出来的那个空间付了�?
解决方法：在主函数定义变量时，定义成指针变量。调用时�?，在子函数的参数里用**�?
*/	
	
	
	
 (*a)=(int *)malloc(sizeof(a));
 **a=100;
 return 1;
}	
	
	
	
	
	
	
	
	
// }

int erci(
		double *H,		//hessian
		double *h,		//原问题grad
		double *be,		//b   等式约束
		double *Ae,		//系数
		double *bi,		////b   不等式约�?
		double *Ai,
		int dim,		//问题的维�?
		int e,			//等式个数
		int ie,			//不等式个�?
		double *xk)
{
//等式约束

//[G,-A;-A' 0] A等式



double *G=(double *)malloc(dim*dim*sizeof(double));


cblas_daxpby(dim*dim, 2, H, 1, 0, G, 1);
int qinum=0;
double ep=1e-14;
//指标�?自动要求等式约束 
int *A0=(int *)malloc((ie)*sizeof(int));
int *tp=(int *)malloc((ie)*sizeof(int));
//A0>=0

double *tg=(double *)malloc(dim*sizeof(double));
double *zuoyong=(double *)malloc((e+ie)*dim*sizeof(double));
double *tb=(double *)malloc((e+ie)*sizeof(double));
double alpha=0;
double *dk=(double *)malloc(dim*sizeof(double));
double *tbi=(double *)malloc(ie*sizeof(double));




double test;
double *ait=(double *)malloc(dim*sizeof(double));
double zuixiao;
int index;

memset(A0,-1,ie*sizeof(int));

memcpy(tb,be,e*sizeof(double));
memcpy(zuoyong,Ae,dim*e*sizeof(double));
memcpy(tg,h,dim*sizeof(double));






memcpy(tbi,bi,ie*sizeof(double));






cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, ie, 1,dim, 1,Ai, dim,xk, 1,-1,tbi,1 );



index=0;  //借用


//printf("Ai\n");
//shuchud(Ai,ie,dim);
//printf("xk\n");
//shuchud(xk,dim,1);

//printf("bi\n");
//shuchud(bi,ie,1);

printf("tbi\n");
shuchud(tbi,ie,1);
//printf("H\n");
//shuchud(H,dim,dim);




for(int i=0;i<ie;i++)
{
	if(tbi[i]<ep) 
			{	
		
				qinum+=1;
				A0[qinum-1]=i;
				for(int tt=0;tt<dim;tt++)
				{zuoyong[dim*e+(qinum-1)*dim+tt]=Ai[dim*i+tt];
			
			
			
				}
				
			}
}



//定义其作用集 函数 由初始点起作用集





 

//xk=tg

/////***//d


int ho=0;

//////diedai
while(ho<30)
{
ho++;
printf("\n%d>>>>>>>>>>>>>>>>>>>>>>>>>>>>d>>>>>>>>>>>>\n",ho);




printf("xk\n");
shuchud(xk,dim,1);

printf("G %lf g %lf dangqiandianxk\n",G[dim*dim-1],h[dim-1]);
shuchud(xk,dim,1);	
	printf("\n\n\n\n\n\n\nqiz\n");
	shuchui(A0,ie,1);
	
	
	
//	printf("..\n");
	printf("\n\n\nqinum\n");
	printf("%d",qinum);


	cblas_daxpby(dim, 1,h, 1, 0, tg, 1);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, dim, 1,dim, 1,G, dim,xk, 1,1,tg,1 );	
	


	memset(tb,0,(e+ie)*sizeof(double));

	printf("\n++...G \n");
	shuchud(G,dim,dim);
	
	printf("++...zuoyong \n");
	shuchud(zuoyong,e+qinum,dim);
	printf("++...tg\n");
	shuchud(tg,dim,1);
	printf("++... %d\n",e+qinum);
	printf("..%lf..\n",tb[0]);




	myqp(H,zuoyong,tg,tb,dim,(e+qinum));

	memcpy(dk,tg,dim*sizeof(double));
		
	printf("++...dk\n");
	shuchud(tg,dim,1);	
	printf("++...lam_ki\n");
	shuchud(tb,e+qinum,1);	
	printf("++...zuoyong \n");
	shuchud(zuoyong,e+qinum,dim);

		



	//dk=0
	if(cblas_dasum(dim, tg,1)<ep)
	{	
		zuixiao=0;
		index=-1;
		for(int i=0;i<qinum;i++)
		{	
			if(tb[e+i]<zuixiao)
				{index=i;
				zuixiao=tb[e+i];
				}	
		}
		if(index>=0)
		{printf("index%d\n",index);
			if(index<qinum-1)
			{
			A0[index]=A0[qinum-1];
			A0[qinum-1]=-1;
			for(int tt=0;tt<dim;tt++)
				{zuoyong[e*dim+index*dim+tt]=zuoyong[e*dim+ (qinum-1)*dim+tt];
				//用最后的覆盖
				}
			qinum=qinum-1;
			}
		 else
			{A0[qinum-1]=-1;
			qinum=qinum-1;}	
		}
		else{
///gai x0weihao			
			
			printf("ffinishd\n");
			return 1;
		}
	}
	else
	{	//alpha xk

		index=-1;
		zuixiao=1;
		alpha=1;
		memset(tp,0,ie*sizeof(int));
		for(int j=0;j<qinum;j++)
		{
			if(A0[j]>=0)
			{tp[A0[j]]=1;}
		}
		//printf("bukexingji....\n");
		//shuchui(tp,ie,1);
		//alpha测试一�?
		for(int j=0;j<ie;j++)
		{	
			//不属于的
			if(tp[j]<1)
			{for(int tt=0;tt<dim;tt++)
					{ait[tt]=Ai[dim*j+tt];}
			
			

			if(cblas_ddot(dim, ait, 1, dk,1)<0)
			{

				test=(bi[j]-cblas_ddot(dim, ait, 1, xk,1));
				test=test/cblas_ddot(dim, ait, 1, dk,1);
				//printf("%lf\n",test);		
					if((test<zuixiao))
					{zuixiao=test;
					index=j;}
			}}
		printf("\n.dangqian.%lf..alpha",zuixiao);
		}


		printf("\n....dk");
		shuchud(dk,dim,1);
		printf("\n");
		printf("\n..%lf..alpha",zuixiao);
		printf("\n");
		
		
		cblas_daxpby(dim,zuixiao,dk,1,1,xk,1);
	/////////



		if(index>-1)
			{qinum=qinum+1;
			A0[qinum-1]=index;
			for(int tt=0;tt<dim ;tt++)
			{zuoyong[e*dim+(qinum-1)*dim+tt]=Ai[index*dim+tt];}
			
			}
	}

	
}
		





free(ait);
free(dk);
free(tb);
	
	
free(tg);
free(zuoyong);
free(A0);
free(G);


}

















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
		)
{
double *bbe=(double *)malloc(0*sizeof(double));
double *aae=(double *)malloc(0*sizeof(double));

double *BM=(double *)malloc((2*e+ie+1)*sizeof(double));
double *AM=(double *)malloc((dim+1)*(2*e+ie+1)*sizeof(double));
double *ait=(double *)malloc(dim*sizeof(double));



/*扩张*/
for(int i=0;i<e;i++)
{	for(int j=0;j<dim;j++)
	{
		AM[i*(dim+1)+j]=Ae[i*dim+j];
	}
}
for(int i=0;i<e;i++)
{	AM[i*(dim+1)+dim]=1;
}

for(int i=0;i<e;i++)
{	BM[i]=be[i];
}

///
for(int i=0;i<e;i++)
{	for(int j=0;j<dim;j++)
	{
		AM[e*(dim+1)+i*(dim+1)+j]=-Ae[i*dim+j];
	}
}
for(int i=0;i<e;i++)
{	AM[e*(dim+1)+i*(dim+1)+dim]=1;
}

for(int i=0;i<e;i++)
{	BM[e+i]=-be[i];
}

///
for(int i=0;i<ie;i++)
{	for(int j=0;j<dim;j++)
	{
		AM[2*e*(dim+1)+i*(dim+1)+j]=Ai[i*dim+j];
	}
}
for(int i=0;i<ie;i++)
{	AM[2*e*(dim+1)+(i+1)*(dim+1)-1]=1;
}

for(int i=0;i<ie;i++)
{	BM[2*e+i]=bi[i];
}


/////////
for(int i=0;i<dim;i++)
{	AM[(2*e+ie)*(dim+1)+i]=0;
}

AM[(2*e+ie+1)*(dim+1)-1]=1;
BM[(2*e+ie)]=0;


double *x0=(double *)malloc((dim+1)*sizeof(double));



/////////////////////////////////////
double *hw=(double *)malloc((dim+1)*sizeof(double));
double *Hw=(double *)malloc((dim+1)*(dim+1)*sizeof(double));
double test;
double ep=1e-15;


double M=1e10;
double t=0;


for(int i=0;i<dim;i++)
{	for(int j=0;j<dim;j++)
	{	Hw[i*(dim+1)+j]=H[i*dim+j];}
}

//Hw �?   M t2+M t


//调用erci 改变了什�?



//shuchud(Ae,e,dim);
//printf("Ae..\n");
//shuchud(Ai,ie,dim);
//printf("Ai..\n");
shuchud(BM,(2*e+ie+1),1);
printf("BM..\n");
shuchud(AM,2*e+ie+1,dim+1);
printf("AM..\n");

shuchud(Hw,dim+1,dim+1);
printf("HW\n");
shuchud(H,dim,dim);
printf("H\n");




for(int i=0;i<e;i++)
{if(fabs(BM[i])>t)
	{t=abs(BM[i]);}}

for(int i=0;i<ie;i++)
{if(bi[i]>t)
	{t=bi[i];}}

//dengshi

if(t>0)
{
x0[dim]=t;	//每次迭代t是满足的
}
else
{x0[dim]=0;
}
///


memset(hw,0,dim*sizeof(double));
hw[dim]=1;


while(x0[dim]>1e-13&&M<1e16)
{

//迭代完成 hw 即x0自动满足约束






M=M*10;


Hw[(dim+1)*(dim+1)-1]=2*M;			//Hw做一下改�?



   //每次都要赋�?
hw[dim]=M;



printf("--------------x0\n");
shuchud(x0,dim+1,1);

erci(Hw,hw,bbe,aae,BM,AM,dim+1,0,(2*e+ie+1),x0);	


}

printf("hw\n\n");
shuchud(hw,dim+1,1);
printf("x0%lf\n\n",log(M)/log(10));
shuchud(x0,dim+1,1);



free(ait);
free(BM);
free(AM);
free(x0);









}













/*


int erciwM(
		double *H,		//hessian
		double *h,		//原问题grad


		double *bi,		////b   不等式约�?
		double *Ai,
		int dim,		//问题的维�?
	
	
		int ie		//不等式个�?原问题的  2e+ie)
)
{
//后续   采用硬盘读取载入内存  此处已经修正	
// double *H,		//hessian
// double *h,		//原问题grad
// double *be,		//b   等式约束
// double *Ae,		//系数
// double *bi,		////b   不等式约�?
// double *Ai,
//增加t
//f(x)+Mt2+Mt
//1e4�? 800M


double t=0;


for(int i=0;i<ie;i++)
{if(abs(bi[i])>t)
	{t=abs(bi[i]);}}




double *hw=(double *)malloc(dim*sizeof(double));

double ep=1e-15;;
double *x0=(double *)malloc(dim*sizeof(double));
memset(x0,0,dim*sizeof(double));




while(t>ep)
{

double *AM=(double *)malloc((2*e+ie)*dim*sizeof(double));
double *BM=(double *)malloc((2*e+ie)*sizeof(double));

double *bbi=(double *)malloc(0*sizeof(double));
double *aai=(double *)malloc(0*sizeof(double));
//无等�?



memcpy(hw,h,dim*sizeof(double));

for(int i=0;i<e*dim;i++)
{AM[i]=Ae[i];
}
for(int i=0;i<e*dim;i++)
{AM[dim*e+i]=-Ae[i];
}
for(int i=0;i<ie*dim;i++)
{AM[dim*e*2+i]=Ai[i];
}


for(int i=0;i<e;i++)
{BM[i]=be[i]-t;
}
for(int i=0;i<e;i++)
{BM[e+i]=-be[i]-t;
}
for(int i=0;i<ie;i++)
{BM[e+i]=bi[i]-t;
}



erci(H,hw,BM,AM,bbi,aai,dim,0,(2*e+ie),x0); //有误

shuchud(hw,dim,1);

memcpy(x0,hw,dim*sizeof(double));
t=t/2;
}






erci(H,h,be,Ae,bi,Ai,dim,e,ie,hw);

}
*/


















int mychol(
		double *L,		//输入 返回 下三�?
		int dim,
		double *Lz		//转好
		)		//维数
{

int s=dim;
int ipiv[s];
int info;
info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,s,s,L,s,ipiv);
//上三�?
//抹去上面  对角线开根号


printf("\n\n\n");
for(int i=0;i<dim;i++)
{
	for(int j=i+1;j<dim;j++)
	{
	L[(i)*dim+j]=0;
	Lz[(j)*dim+i]=0;
	}

	for(int j=0;j<i;j++)
	{
	L[(i)*dim+j]=L[(i)*dim+j]*sqrt(L[(j)*dim+j]);
	Lz[(j)*dim+i]=L[(i)*dim+j];
	}

}


//最后再变化对角�? 不能影响  下三角时以列为主循环
for(int i=0;i<dim;i++)
{L[(i)*dim+i]=sqrt(L[(i)*dim+i]);
Lz[(i)*dim+i]=L[(i)*dim+i];

}


	
//得下三角

	
}


int myqp(
		double *H,		//hessian
		double *Aw,		//yueshu
		double *gk,
		double *b,
		int dim,		//G维数
		int   ge			//A的列�? 应该改过方向�?
		)
{
double *AM;
double *A=(double *)malloc(ge*dim*sizeof(double));
memcpy(A,Aw,ge*dim*sizeof(double));

int tte=xxwg(AM,A, ge, dim);
printf("000000000000000000000000000000000000000000000000000000000000step1\n");

if(tte<ge)
{	
	AM=(double *)malloc(ge*(dim+1)*sizeof(double));
	for(int qq=0;qq<ge;qq++)
	{	AM[qq*(dim+1)+dim]=b[qq];

		for(int kk=0;kk<dim;kk++)
		{
			AM[qq*(dim+1)+kk]=Aw[qq*(dim)+kk];
			
		}
		
		
		
	}
	
	
	
	if(xxwg(A,AM, ge, dim+1)>tte)
	{printf("wujie");
		return 1;
	}
	
	
	
	memcpy(AM,Aw,ge*dim*sizeof(double));

	tte=xxwg(A,AM, ge, dim);




	
	
	
	
	//xishu zhen  增广阵的�?
}


shuchud(A,tte,dim);
int e=tte;




//qiu Ajidaxianxingwuguanzu 




	
//复制G
double *G=(double *)malloc(dim*dim*sizeof(double));
double *tgk=(double *)malloc(dim*sizeof(double));

cblas_daxpby(dim*dim, 2, H, 1, 0, G, 1);
if(e==0)
{	//无约束还要修�?  通过ll分解判断正定
	ni(G,dim);
	memcpy(tgk,gk,dim*sizeof(double));
	cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, dim, 1,dim, -1,G, dim,tgk, 1,0,gk,1 );
	free(G);
	free(tgk);
	return 1;}



double *L=(double *)malloc(dim*dim*sizeof(double));
double *LW=(double *)malloc(e*e*sizeof(double));
double *Lz=(double *)malloc(dim*dim*sizeof(double));
double *LWz=(double *)malloc(e*e*sizeof(double));


double *TEM=(double *)malloc(dim*dim*sizeof(double));
double *GI=(double *)malloc(dim*dim*sizeof(double));
double *TEM2=(double *)malloc(dim*e*sizeof(double));

double *TEM3=(double *)malloc(e*e*sizeof(double));

double *V=(double *)malloc(e*e*sizeof(double));


memcpy(L, G, dim*dim*sizeof(double));	
mychol(L,dim,Lz);
//求下三角


	
memcpy(GI, G, dim*dim*sizeof(double));	
ni(GI,dim);




 

//求�?


cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasTrans, dim, e,dim, 1,GI, dim,A, dim,0,TEM2,e );


cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, e, e,dim, 1,A, dim,TEM2, e,0,TEM3,e );




memcpy(V, TEM3, e*e*sizeof(double));
memcpy(LW, TEM3, e*e*sizeof(double));


mychol(LW,e,LWz);
	

	
	
//求V L波浪
	

double *u=(double *)malloc(dim*sizeof(double));
	
cblas_daxpby(dim, -1, gk, 1, 0, u, 1);
//u=-gk


int info;
int ipiv[dim],ipive[e];

//dim   yigelie
//要转一起转�?
//dgesv会破坏原矩阵



memcpy(TEM,L,dim*dim*sizeof(double));
LAPACKE_dgesv(LAPACK_ROW_MAJOR,dim,1,TEM,dim,ipiv,u,1);
memcpy(TEM,Lz,dim*dim*sizeof(double));
LAPACKE_dgesv(LAPACK_ROW_MAJOR,dim,1,TEM,dim,ipiv,u,1);
//-gk u w 




//bw v lam
double *bw=(double *)malloc(e*sizeof(double));
cblas_daxpby(e, 1, b, 1, 0, bw, 1);
cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, e, 1,dim, -1,A, dim,u, 1,1,bw,1 );	


memcpy(TEM3,LW,e*e*sizeof(double));
LAPACKE_dgesv(LAPACK_ROW_MAJOR,e,1,LW,e,ipive,bw,1);
memcpy(TEM3,LWz,e*e*sizeof(double));
LAPACKE_dgesv(LAPACK_ROW_MAJOR,e,1,LWz,e,ipive,bw,1);



//-gk hw y x
cblas_daxpby(dim, -1, gk, 1, 0, u, 1);





cblas_dgemm(CblasRowMajor, CblasTrans,CblasNoTrans, dim, 1,e, 1,A, dim,bw, 1,1,u,1 );


memcpy(TEM,L,dim*dim*sizeof(double));
LAPACKE_dgesv(LAPACK_ROW_MAJOR,dim,1,TEM,dim,ipiv,u,1);
memcpy(TEM,Lz,dim*dim*sizeof(double));
LAPACKE_dgesv(LAPACK_ROW_MAJOR,dim,1,TEM,dim,ipiv,u,1);

//printf("lambda\n");
//shuchud(bw,e,1);

//printf("x*\n");
//shuchud(u,dim,1);
	
memcpy(gk,u,dim*sizeof(double));		
memcpy(b,bw,dim*sizeof(double));		
	
	
	
free(bw);
free(u);	

free(G);
free(GI);

free(L);
free(LW);
free(Lz);
free(LWz);
free(TEM);
free(TEM2);
free(TEM3);
free(V);

	
	
return 0;	
}




//有约束优�?局部sqp    假设梯度已知  拟牛�?powell

// void(*pfunarr[3])();
// 函数指针数组
// void(*(*pfunarr2)[3])();
// 指向函数指针数组的指�?

//








//f(X1,X2,X3...X_gex;P1,P2,P3...,P_gep)
//[x;u;t]    f(x,u;t)
//
//f(向量指针列表)
//listx[i]   Xi deweishu
//f({1,2,    3})





//当函数给定后 分量的排列是确定�?




// x=(double **)malloc(gex*sizeof(double *));


// for(i=0;i<gex;i++)
// {
// x[i]=(double *)malloc(listx[i]*sizeof(double));



int X2x(double *X,double **x,int *listx,int gex)
{
	
int i,j,dangqian=0;

for(i=0;i<gex;i++)
{
for(j=0;j<listx[i];j++)
x[i][j]=X[dangqian+j];

dangqian+=listx[i];
}
	
	
return 1;	
	
	
}


//有带待充系数   q第几�?
//第一个l=1   取高斯点�?t_0 为左端点

//////int XP2xpl(double *X,double *P,double **xp,int *listx,int gex,int *listp,int gep,int l)
//{
//	
//int i,j,dangqian=0;
//int zonx=0,zonp=0;
//
//for(i=0;i<gex;i++)
//zonx+=listx[i];
//
//for(i=0;i<gep;i++)
//zonp+=listp[i];
//
//
//
//dangqian=zongx*(l-1);
//
//
//
//
//
//
//
//for(i=0;i<gex;i++)
//{
//for(j=0;j<listx[i];j++)
//xp[i][j]=X[dangqian+j];
//
//dangqian+=listx[i];
//}	
//	
//	
//dangqian=zongp*(l-1);
//
//	
//	
//for(i=0;i<gep;i++)
//{
//for(j=0;j<listp[i];j++)
//xp[i+gex][j]=P1[dangqian+j];
//
//dangqian+=listp[i];
//}	
//	
//	
//return 1;	
//}





double fxutRnR(lyRnR f,double *X,int *listx,int gex)
{

int i,j,k;

int dangqian=0;

double **x,rs;

x=(double **)malloc(gex*sizeof(double *));


for(i=0;i<gex;i++)
{
x[i]=(double *)malloc(listx[i]*sizeof(double));
for(j=0;j<listx[i];j++)
x[i][j]=X[dangqian+j];

dangqian+=listx[i];
}

rs=f(x);

free(x);
x=NULL;

return rs;



}




////指标转换
//
//int v2tt(int v,int dim,int dimu,int N)
//{
//int tt=ceil(v/(dimx+dimu));
//return tt;	
//}
//
//int v2qp(int v,int dim,int dimu,int N)
//{
//int qp=v%(dimx+dimu+1);
//return qp;	
//}
// 
//
//int v2s(int v,int dim,int dimu,int N)
//{
//int qp=v%(dimx+dimu+1);
//
//int s=1;
//
//if qp>dimu;
//s=2;
//
//return s;	
//}
//
//int v2qx(int v,int dim,int dimu,int N)
//{
//int qp=v%(dimx+dimu+1);
//int qx=qp;
//
//
//int s=1;
//if qp>dimx;
//{
//	qp=qp-dimx;
//	
//	
//}
//
//
//
//return qx;	
//}
//
//
//
//
//
//
//
//
////partial yita F   partial f
//
////partial keis_k F   partial f
//
////逐点函数到全决策变量函数（xi,ui,tf�?
//
//
//// RnR
////3�?	gg()
//
//int gg2gk(lyRnRnk gxu,double *XUtfk, int dimx,int dimu,   double *tk,int N,int k,double *gk)
//{
//	
//double**  x=(double **)malloc(3*sizeof(double *));
//int i,j;
//
//x[0]=(double *)malloc(dimx*sizeof(double ));
//
//for(i=0;i<dimx;i++)
//{
//x[0][i]=XUtfk[(dimx+dimu)*(k-1)+i];
//
//}
//
//
//for(i=0;i<dimu;i++)
//{
//x[1][i]=XUtfk[(dimx+dimu)*(k-1)+dimx+i];
//
//}
//x[2][0]=tk[k-1];
//
//gk=gxu(x);
//free(x);
//x=NULL;
//
//
//return 1;
//}
//
//
//
////tf 偏导
//int gt2gitf(lyRnRnk gg,double *XUtfk, int dimx,int dimu,   double *tk, double *tauk, int N,int k,double *gk)
//{  
// 
//  
//double**  x=(double **)malloc(3*sizeof(double *));
//int i,j;
//
//x[0]=(double *)malloc(dimx*sizeof(double ));
//x[1]=(double *)malloc(dimu*sizeof(double ));
//x[2]=(double *)malloc(sizeof(double ));
//
//
//
//for(i=0;i<dimx;i++)
//{
//x[0][i]=XUtfk[(dimx+dimu)*(k-1)+i];
//
//}
//
//
//for(i=0;i<dimu;i++)
//{
//x[1][i]=XUtfk[(dimx+dimu)*(k-1)+dimx+i];
//
//}
//x[2][0]=tk[k-1];
//
//gk=gt(x)*(tauk[k]/2+1.0/2);
//free(x);
//x=NULL;
//
//
//return 1;
//   
//    
//}  
//  
//   
//
//
////性能指标的偏�?
////常微离散化偏�?
//
////逐点等式
////mayer 等式
////逐点 不等�?
////mayer 不等式个�?
////tf 偏导
//
//

int X2xutk(double *X,double **x,int dimx,int dimu,double *tk,int k)
{
int i;
x[0]=(double *)malloc(dimx*sizeof(double ));

for(i=0;i<dimx;i++)
{
x[0][i]=X[(dimx+dimu)*(k-1)+i];

}


for(i=0;i<dimu;i++)
{
x[1][i]=X[(dimx+dimu)*(k-1)+dimx+i];

}
x[2][0]=tk[k-1];

return 1;
}



int xishu(lyRnR PHI,       //目标函数中的终端
		lyRnR g,     //目标函数中的被积函数
		lyRnRnk f,
		

		lyRnR *cek,
		int numcek,		//逐点等式约束
		
		lyRnR *phii,
		int numphii,	//mayer终端等式约束的个�?
		
		lyRnR *cik,
		int numcik,		//逐点不等式约�? cik<=0
		
		
		lyRnR *psii,
		int numpsii,	//不等式终�?


		piandao gPHI,			//PHI(x,t) 返回PHI_x,PHI_t
		
		
		piandao gg,	 //被积函数的梯�?
		//gg(**(x,u,t))  返回(**( gg_kesi,gg_t ))		
		
		piandao gf,	 //状态方程函数的梯度
		
		piandao *gcek,	//逐点的梯度函�?
		
		piandao *gphi,
					
		piandao *gcik,
		
		piandao *gpsi,				//mayer 的梯�?
		
	
		
		int dimx,
		int dimu,
		int Ntau,  //几阶勒让德方�?
		double *tauk, //勒让德点
		double *wk, //高斯勒让德积分系�?
		double *Dki, //导数系数
		double *zuiyouX, //最优结�?
		 
	                          	//  double *H, H由这些和子问题等等进行修�? Bk+1
		double *h,
		double *be,	
		double *Ae,	
		double *bi,
		double *Ai,
		double *XUk,
		double tk,
		double *x0,
		double t0
)
{
	
double **xftf;








//xf=x0+....;
double  tf=XUk[Ntau*(dimx+dimu)];	
double tt2=(tf-t0)/2;	
	
	
	
	
//不转化成f(x)  减少不需要的赋�? dairu(f,xuk ,....)	
	
int i,j,k,lo,mu; 
int zong=(dimx+dimu)*Ntau+1;
double **xk,ttem;
double *temg,*temg2,*temf,*temf2,*temnew,*temPHI1,*temPHI2;
double *tem1,*tem2;

cshi(temg,1);
cshi(temg2,1);
cshi(temf,1);
cshi(temf2,1);
 
cshi(temnew,1);
 
cshi(temPHI1,1);
cshi(temPHI2,1);
 
cshi(tem1,1);
cshi(tem2,1);
 









double *jiluxf,*jiluliang;
//初始�?


double *temkesi,*temxu,*temxu1;




















int L;






xk=(double **)malloc(3*sizeof(double *));
xk[0]=(double *)malloc(dimx*sizeof(double ));
xk[1]=(double *)malloc(dimu*sizeof(double ));
xk[2]=(double *)malloc(sizeof(double ));















//状态方程的等式和梯度的系数

//梯度�?
//
//对角  O	
//				-P_kesil f
//对角  O	
//

//形式







//be维数

//Ntau*dimx
//+Ntau*numce
//+numphii	
	
	
//bi维数
	
//Ntau*numci
//+numpsii	

	


//cek
//当前�?
int dangqianhe=Ntau*dimx;

//不等式的�?
int dangqianhi=0;

//尽可能复用一些运�?
int	maxnine=numcek;



h[Ntau*(dimx+dimu)]=0;






for(i=0;i<Ntau;i++)
{
X2xutk(XUk,xk,dimx,dimu, tauk,i);




/*因为可导  指标 仅需要知道梯�?/
//积分部分的梯�?
//两部�?对状态控�?  对tf
//mayer的梯度和函数值涉�?xf  需要计算所有的   fk




//kesi 
temg=gg(xk,1);
cblas_daxpby((dimx+dimu), tt2*wk[i],temg , 1, 0, temxu, 1);
for(L=0;L<dimx+dimu;L++)
{
h[i*(dimx+dimu)+L]+=temxu[L];
}
//wk gk 部分



free(temg);
temg=NULL;
temg=gg(xk,2);
//p_tf  jifen   对tf的导�?要对 g_t修正

h[Ntau*(dimx+dimu)]+=(0.5*g(xk)+tt2*temg[0]*(tauk[i]/2+0.5))*wk[i];



///////////////////////////////////////////////////////////////////////


//状态方�?
// 终端函数的偏导数里有�?
//列向�?
//temf在循环外先定义一个东�?

free(temf);
temf=NULL;



temf=f(xk);

cblas_daxpby(dimx, tt2*wk[i],temf , 1, 1, jiluxf, 1);
	


//第一组等式约�?
for(L=0;L<dimx;L++)
{
be[i*(dimx)+L]+=tt2*temf[L];
	
	
	
	for(mu=0;mu<Ntau;mu++)
	{be[mu*(dimx)+L]+=-Dki[mu*(Ntau)+i]*xk[0][L]�?
	}

}



//求时间的偏导�?

free(temg)
temg=NULL;
temg=gf(xk,2);

for(L=0;L<dimx;L++)
{
be[i*(dimx)+L]+=(tauk[i]+1)/2*temg[L];
}

cblas_daxpby(dimx, tt2*wk[i],temg , 1, 1, jiluliang, 1);
	
//tt2*wk[i] 这个也可以优化减少乘次数  增加时间复杂�?


//xu偏导数在下一个循环里




//////////////////////////////////////////////////////////////////
//第二 四组
//系数和梯度都可以计算
dangqianhe=dimx*Ntau;


for(j=0;j<numcek;j++)
{
free(temg)
temg=NULL;



temg=gcek[j](xk,1);
be[dangqianhe+i*numcek+j]=-cek(xk);
//返回 P_kesi P_t

for(k=0;k<(dimx+dimu);k++)	
Ae[(dangqianhe+i*numcek+j)*zong+  i*(dimx+dimu)+k]=temg[k];


free(temg)
temg=NULL;
temg=gcek[j](xk,2);

//一维数�?
Ae[(dangqianhe+i*numcek+j)*zong+  N*(dimx+dimu)]=temg[0]*(tauk[i]/2+0.5);
}







dangqianhi=dimx*Ntau+numci*Ntau+numphii;

for(j=0;j<numcik;j++)
{
free(temg)
temg=NULL;



temg=gcik[j](xk,1);
bi[dangqianhi+i*numcik+j]=-cik(xk);
//返回 P_kesi P_t

for(k=0;k<(dimx+dimu);k++)	
Ai[(dangqianhi+i*numcik+j)*zong+  i*(dimx+dimu)+k]=temg[k];


free(temg)
temg=NULL;
temg=gcik[j](xk,2);

//一维数�?
Ai[(dangqianhi+i*numcik+j)*zong+  N*(dimx+dimu)]=temg[0]*(tauk[i]/2+0.5);
}

}
/*--------------------------------------------------------------------------------------------------------*/


memcpy(x0,temf,dimx*sizeof(double));
cblas_daxpby(dimx, 1,jiluxf , 1, -1, temf, 1);
cblas_daxpby(dimx, 1/(tf-t0),temf , 1, 1, jiluliang, 1);


//jiluliang+=jiluxf-x0/tf-t0;


//终端关于时间  的偏导数有几项可以先计算
//



//终端约束


 //目标函数的梯度修�?  
	free(temPHI1);
	temPHI1=NULL;
	temPHI1=gPHI(xftf,1);
	
	free(temPHI2);
	temPHI2=NULL;
	temPHI2=gPHI(xftf,2);
 
	
	





dangqianhe=Ntau*dimx+Ntau*numcek;
for(k=0;k<numphii;k++)
{
	ttem=phii[k](xftf);
	be[dangqianhe +k  ]=ttem;	
	free(temg);
	temg=NULL;
	temg=gphi[k](xftf,1);
		for(mu=0;mu<dimx;mu++)
			Ae[(dangqianhe+k)*zong+mu]=temg[mu];
	free(temg);
	temg=NULL;
	temg=gphi[k](xftf,2);
	be[(dangqianhe+k)*zong+dimx]=temg[0];
	
}	
	
 

dangqianhi= Ntau*numcik;
for(k=0;k<numpsii;k++)
{
	ttem=psii[k](xftf);
	bi[dangqianhi +k  ]=ttem;	
	free(temg);
	temg=NULL;
	temg=gpsi[k](xftf,1);
		for(mu=0;mu<dimx;mu++)
			Ai[(dangqianhi+k)*zong+mu]=temg[mu];
	free(temg);
	temg=NULL;
	temg=gphi[k](xftf,2);
	bi[(dangqianhe+k)*zong+dimx]=temg[0];
	
}	
	











/*--------------------------------------------------------------------------------------------------------*/
	
//计算遗留的系�?
//终端约束 终端约束的梯�?
//偏f
	
for(k=0;k<Ntau;k++)
{
	
X2xutk(XUk,xk,dimx,dimu, tauk,i);
	
temg=gf(xk,1);


for(mu=0;mu<(dimx);mu++)	
{	Ae[(k*dimx+mu)*zong+  mu]+=Dki[k*(Ntau)+k];

	for(lo=0;lo<(dimx+dimu);lo++)	
		Ae[(k*dimx+mu)*zong+  lo]+=-tt2*temg[mu*(dimx+dimu)+lo];
}



//先算状态方�?然后其他终端的的约束计算过程应该类似

 
	 
	
//指标 只要梯度 偏PHI	
//状�?

free(tem1);
tem1=NULL;

tem1=(double *) malloc(dimx*sizeof(double));
free(tem2);
tem2=NULL;

tem2=(double *) malloc( sizeof(double));

 
 
 
 
 
 
 cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,1, (dimx+dimu), dimx,tt2*wk[k],temPHI1, dimx,temg,dimx+dimu, 0,temxu, dimx+dimu);
 
 cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,1, 1,dimx, 1,temPHI1, dimx,jiluliang,1, 1,temPHI2, 1);
	
	for(mu=0;mu<dimx;mu++)	
		h[k*(dimx+dimu)+mu]+=temxu[mu];
	
	h[Ntau*(dimx+dimu)]+=temPHI2[0];
	
 
 
 
 
 
 
 
 
 
for(lo=0;lo<numphii;lo++)
{	for(mu=0;mu<dimx;mu++)
		tem1[mu]=Ae[(dangqianhe+lo)*zong+mu];
	tem2[0]=Ae[(dangqianhe+lo)*zong+dimx];
	
		/*..............................................................*/

		cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,1, (dimx+dimu), dimx,tt2*wk[k],tem1, dimx,temg,dimx+dimu, 0,temxu, dimx+dimu);
	for(mu=0;mu<(dimx+dimu);i++)
	Ae[(dangqianhe+lo)*zong+mu]=temxu[mu];
	
		cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,1, 1,dimx, 1,tem1, dimx,jiluliang,1, 1,tem2, 1);
		
	Ae[(dangqianhe+lo)*zong+dimx+dimu]=tem2[0];
}

//P kesiL  f   这个雅可比矩阵在mayer 指标 和mayer 等式不等式中都用�?
//按照雅可比阵�? 行列关系
//P_1 PHI 为行矩阵
//mayer 的输入为（xf,tf�?
//返回

//f 的第一个是雅可比阵
//gf 返回一个二级指�?指向两个一级指针第一个是雅可比阵d�?(dx*du）列


 
for(lo=0;lo<numpsii;lo++)
{	for(mu=0;mu<dimx;i++)
		tem1[mu]=Ai[(dangqianhi+lo)*zong+mu];
	tem2[0]=Ai[(dangqianhi+lo)*zong+dimx];
	
		/*..............................................................*/

		cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,1, (dimx+dimu), dimx,tt2*wk[k],tem1, dimx,temg,dimx+dimu, 0,temxu, dimx+dimu);
	for(mu=0;mu<(dimx+dimu);i++)
	Ai[(dangqianhi+lo)*zong+mu]=temxu[mu];
	
		cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,1, 1,dimx, 1,tem1, dimx,jiluliang,1, 1,tem2, 1);
		
	Ai[(dangqianhi+lo)*zong+dimx+dimu]=tem2[0];
}

	
	
	
	
	
}









}




/*。。。。。。。。。。。。。。。。。。。。。。。。。。。。。�?/
//计算H B k+1












// // // //sum Dki xi-
// // // for(mu=0;mu<dimx;mu++)
// // // {
	// // // for(k=0;k<Ntau;k++)
	// // // {
		// // // be[i*dimx+mu]-=Dki[i*Ntau+k]*XUk[k*(dimx+dimu)+mu];	
	// // // }
// // // }






















   



//控制问题优化  自由tf 和x_tf


//优化变量(x1,u1,...,xn,un,tf)



//f目标函数
//gf 梯度
//cek  等式约束
//cik  不等约束
//gcek  等式约束梯度
//gcik  不等约束梯度

//先不考虑积分形过程约�? PHI 及phi 都是 xtf tf的函�?
//


int kzzy(lyRnR PHI,       //目标函数中的终端
		lyRnR g,     //目标函数中的被积函数
		lyRnR f,
		
		lyRnR *phii,
		int numphii,	//mayer终端等式约束的个�?
		lyRnR *cek,
		int numcek,		//逐点等式约束
		
		
		lyRnR *psii,
		int numpsii,	//不等式终�?

		lyRnR *cik,
		int numcik,		//逐点不等式约�? cik<=0
		
		lyRnRn gPHI,			//PHI(x,t) 返回PHI_x,PHI_t
		
		
		lyRnRn gg,	 //被积函数的梯�?
		//gg(**(x,u,t))  返回(**( gg_kesi,gg_t ))		
		
		lyRnRn gf,	 //状态方程函数的梯度
		
		lyRnRn *gphi,
		
		
		lyRnRn *gcek,	//逐点的梯度函�?
		
		lyRnRn *gpsi,				//mayerde 梯度
		
		
		lyRnRn *gcik,
		
		int dimx,
		int dimu,
		int Ntau,  //几阶勒让德方�?
		double *tauk, //勒让德点
		double *wk, //高斯勒让德积分系�?
		double *Dki, //导数系数
		double *zuiyouX, //最优结�?
		 
		double *H,
		double *h,
		double *be,	
		double *Ae,	
		double *bi,
		double *Ai,
		double *XUk
		
		)
{
int i;
int dimxu=dimx+dimu;
int dimxu_tf=(dimx+dimu)*Ntau+1;
	
	
	
double *Bk;
//Bk=(double *)malloc(dim*dim*sizeof(double)); 
//memset(Bk,0,dim*dim*sizeof(double));
//for(i=0;i<dim;i++)
//Bk[i*dimxu_tf+i]=1;
	
	
	
//s1  初值为0
//s2 子问�?
//















	
	return 0;
	
	
}






								














/*







//min f(x_i^n_j)     f(x,u,v,x^1,u^1,v^1)  i代表 x,u,v    n代表不同的时�? j 代表每个列向�? x,u,v  的分�?


// ni多少种列向量  即nf的维�?  nt 分割的时�?  nf  每种列向量的维数构成数组 
//x0 输入初�?  维数（各种向量维数的和）*时刻
//输出最优点�?

//f目标函数
//gf 梯度
//cek  等式约束
//cik  不等约束
//gcek  等式约束梯度
//gcik  不等约束梯度
//









int yueshu(lyRnR f,lyRnRn gf,lyRnR *cek,lyRnR *cik,lyRnRn *gcek,lyRnRn *gcik,int nt,int *nf,int ni,double *x0)
{
	

int dim=0,i;

for(i=0;i<ni;i++)
{
	dim+=nf[i];	
}
dim*=nt;
//x0的维�?

double *Bk;


Bk=(double *)malloc(dim*dim*sizeof(double)); 
memset(Bk,0,dim*dim*sizeof(double));
for(i=0;i<dim;i++)
Bk=1;
//单位�?








//转化成系系数�?







erci(double *H,	double *h,	double *be,	double *Ae,	double *bi,	double *Ai,		int dim,int e,	int ie,		double *xk);
		

//jiaozheng








	



	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
double
	
	
	
	
	
	
	
double x[2]={1.0,0};	
 
	
	
	
return 0;
	
}



//系数转化





//int *lf,int np,int *lp    (X1,...X_ni;P_1,...P_np)  每个的维数lf   



//tk 不含端点  TN总共的分�?
//有大M法求初�? krylov 法算线性方程组的解

//先固定端�?考虑逐点约束
//分别有分离的逐点约束 过程�?泛函约束  有些终端条件

//int gshuce,int gshuci,int *lxce,int lxci	
//等式不等式约束的总个�?    逐点约束的个�?过程终端 约束的个�?      定义梯度函数时可以分�?


//先做固定首末端问�?


//假定只有终端状态未�?



int yueshuxishu(lyRnR f,lyRnRn gf,lyRnR *cek,lyRnR *cik,lyRnRn *gcek,lyRnRn *gcik,int ni, int *lf,int np,int *lp  ,    double *tauk,int tauN,double *wk,double t0,double tf,
double *H,	double *h,	double *be,	double *Ae,	double *bi,	double *Ai,	double *XUk)
{
	
	
	
//x u t所有的维数
//x,u;t   np=1
 	
int i,j,k;

int zonx=0,zonp=0,dangqian=0;
double **xut,*temg;

xut=(double **)malloc((ni+np)*sizeof(double *));

for(i=0;i<ni;i++)
{zonx+=lf[i];
xut[i]=(double *)malloc(lf[i]*sizeof(double ));
}


for(i=0;i<np;i++)
{zonp+=lp[i];
xut[ni+i]=(double *)malloc(lf[i]*sizeof(double ));
}




temg=(double *)malloc((zonx+zonp)*sizeof(double ));






int dim=(zonx+zonp)*TN;

	
int i;
for(i=0;i<tauN;i++)
{
//若只有逐点约束只要求某一时刻的量



XP2xpl(XUk,Tk,xut,lf,ni,lp,np,i);
 
 
 
 
 
 
 
temg=gf(xut);



for  
   XUk...=



be=-c(xut)
bi=-ci(xut)

be=-phi(XUk)
bi=-phi(XUk)




Ae
Ai



}	
	
	
	
	
	
	
}

*/
























int fft(double _Complex  * ai,double _Complex  * ao,int N)
{
	 fftw_complex *in, *out;
    fftw_plan p;
   
    int i;
    int j;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	
	
        for( i=0; i < N; i++)
    {
            in[i][0] = creal(ai[i]);
            in[i][1] = cimag(ai[i]);
            printf("%6.2f ",in[i][0]);
    }
	
    printf("\n");
	
    p=fftw_plan_dft_1d(N,in,out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p); /* repeat as needed*/
        
	for(j = 0;j < N;j++)
    {
			ao[j]=out[j][0]+out[j][1]*I;
		shuchuz(ao,N,1);
    }
    printf("\n");
    
	
	fftw_destroy_plan(p);
    fftw_free(in); 
    fftw_free(out);
    return 0;
	
	
 }


















/*


//拉格朗日插值导数Dki 一块算 减少重复计算

int dlagx(double x,double  *c,dim NN,double * Dki)
{

double fenmui=1;

int j;
for(j=0;j<NN;j++)
{if(i!=j)
{fenmui*=(c[i]-c[j]);}}





for(k=0;k<NN;k++)
{for(i=0;i<NN;i++)
{
if(i!=k)	
{
	for(j=0;j<NN;j++)
	{
		
	}
}
else
{}

}
}	
	
	
	
	
	
	
	
}






 
//拉格朗日插值基函数






int lagx(double x,int i,double  *c,int N)
{
	
	
	
	
	
	
	
	
	
	
}




double ttau(double tau,double t0,double tf)
{
	return  (tf-t0)/2*tau+(tf+t0)/2
}






//高斯勒让德积�?
int glint(double *(*fff)(double ,double *),double* Xk,double* Uk,double tauk,int N)
{
	
	
	wk*fff+
	
	
	
	
	
	return 0;
}





//切比雪夫高斯�?

int cgd(double *x,int K)
{
int i;
for(i=1;i<=K;i++)
{
x[i-1]=cos((K+1-i)/(K+1)*PI);
	
	
}

return 0;
}	








*/

//勒让德高斯点
//x0,x1...,xn    初值为0
//2n+1次代数精�?
//x n+1�?

//输入节点指针和n
//输出节点不含端点


int lgd(double *d,int n)
{



int i;
//-1+i*h

//s1 求yk
//
 
memset(d,0,(n+1)*sizeof(double));

double *e=(double *)malloc((n)*sizeof(double));


double *z=(double *)malloc((n+1)*(n+1)*sizeof(double));
int ldz=n+1;



for(i=1;i<=n;i++)
{



e[i-1]=sqrt((2*i-1.0)/(2*i+1.0))*i/(2*i-1.0);


}

printf("\n");

 
//方向  job ’N‘只求特征�?’v�? 还要特征向量

//对称三对角特征�?
//z 特征向量


char job = 'N';

LAPACKE_dstevd( CblasRowMajor, job, n+1,  d,e, z, ldz );

 


 


return 0;
}	





//勒让德高斯系数Ak  Dki
//A0...An   导数系数   d Li(xk)



int lg_AD(double *xk,double *Ak,double *Dki,int n)
{
double dp_1=0,dp_2=1,x,dp,p_1,p_2,p;




int i,j,k;
if(n==0)
{
Ak[0]=2;
Dki[0]=0;



return 0;
}
else if(n==1)
{
Ak[0]=1;
Ak[1]=1;


Dki[0]=xk[0]/(1-xk[0]*xk[0]);
Dki[1*(n+1)+(1)]=xk[1]/(1-xk[1]*xk[1]);


Dki[1*(n+1)+(0)]=1.0/(xk[1]-xk[0]);

Dki[0*(n+1)+(1)]=1.0/(xk[0]-xk[1]);



return 0;}













for(i=0;i<=n;i++)
{	x=xk[i];

	dp_2=0;
	dp_1=1;
	p_2=1;
	p_1=x;
	
	//复用dp
	
	


	
	for(j=2;j<=n+1;j++)	
		{
			//dp=p'j
			dp=((2*j-1.0)*x*dp_1-(j-1.0)*dp_2+(2*j-1.0)*p_1)/j;
			dp_2=dp_1;
			dp_1=dp;
			
			p=((2*j-1.0)*x*p_1-(j-1)*p_2)/j;
			p_2=p_1;
			p_1=p;
		}

	Ak[i]=2/(1-x*x)/(dp*dp);
	
	
	
	//保存 Pn'(xi)
	Dki[i]=dp;
	
	
	
}


//由Pn'(xk)  Pn'(xk)  算下三角

for(i=0;i<=n;i++)
{



//k 行i�?
for(k=i+1;k<=n;k++)	
{
Dki[k*(n+1)+i]=Dki[k]/Dki[i]/(xk[k]-xk[i]);
}





//对角�?
Dki[i*(n+1)+i]=xk[i]/(1-xk[i]*xk[i]);
}




//下三角对称覆盖上三角
for(i=0;i<=n;i++)
for(k=i+1;k<=n;k++)	
Dki[i*(n+1)+k]=Dki[k*(n+1)+i];






return 0;
}	





int cshi(double *x,int n)
{
	x=(double *)malloc(n*sizeof(double));
	memset(x,0,n*sizeof(double));
	return 1;
}









//
// 
//
//
///*变址计算，将x(n)码位倒置*/  
//int changema(double* x,int size_x)        
//{  
//  double  temp;  
//  unsigned short i=0,j=0,k=0;  
//  double t;  
//  for(i=0;i<size_x;i++)  
//  {  
//shuchud(x,size_x,1);
//
//
//    k=i;j=0;  
//    t=(log(size_x)/log(2));  
//  while((t--)>0 )    //利用按位与以及循环实现码位颠�? 
//  {  
//    j=j<<1;  
//    j|=(k & 1);  
//    k=k>>1;  
//  }  
//  if(j>i)    //将x(n)的码位互�? 
//  {  
//  printf("-----------\n");
//  shuchud(x,size_x,1);
//    temp=x[i]; 
//	 printf("%lf,\n",temp);
// 
//    x[i]=x[j];  
//printf("%lf,\n",temp);
//    x[j]=temp;
//printf("####%lf,%lf####%d",temp,x[3],j);
//
//
//shuchud(x,size_x,1);	
//  printf("++++++%d++\n",i);
//   } 
//   
//   
//   
//    
//  
//  
//  
//   shuchud(x,size_x,1);
//   return 1;
//
//  
//}  
//
//
//
//
//
//
//















// for(int i=0;i<s;i++)
// {shuchud(Kn_m[i],dim,1);	
// printf("\n");
// }
	








// for(i=0;i<m;i++)  
    // {
        // for(j=0;j<n;j++)  
        // {
            // printf("%p\n",&a[i][j]);     //输出每个元素地址，每行的列与列之间的地址时连续的，行与行之间的地址不连�?
        // }
    // }
    // for(i=0;i<m;i++)  
    // free(a[i]);
 
    // free(a);  
// --------------------- 
// 作者：阿阿阿阿阿阿�?
// 来源：CSDN 
// 原文：https://blog.csdn.net/fengxinlinux/article/details/51541003 
// 版权声明：本文为博主原创文章，转载请附上博文链接�?








	// OO=(double *)malloc(n*n*sizeof(double));
	// memset(OO,0,n*n*sizeof(double));  


