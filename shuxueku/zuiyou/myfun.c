#include <stdio.h>
#include <string.h>
#include "fftw3.h"

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



//改写为带控制的


double  *ode1(double a,double b,int n,double *(*fff)(double ,double *),double *x0,int dim,double *jieguo)
{/*a初值,b终值,n分点,fff连续函数*/   

FILE *logf;
//记录用


double h = (b - a)/n;
//划分

double t0=a,tk;



double* yk =(double *)malloc(dim*sizeof(double));
double* q1 =(double *)malloc(dim*sizeof(double));
double* q2 =(double *)malloc(dim*sizeof(double));
double* q3 =(double *)malloc(dim*sizeof(double));
double* q4 =(double *)malloc(dim*sizeof(double));

double* tem=(double *)malloc(dim*sizeof(double));


//*a要分配内存 或者定义时候分配
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
{	

//用渐进公式  截断多项式给初值  通过牛顿迭代 梯度下降   求精确解





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

 //printf("预估值%.6f,%.5f\n\n\n",w1,w2);



y=y0;



while(cabs(a-y*cexp(y))>1e-12)
{	
w1=creal(y);
w2=cimag(y);
// printf("当前点%lf,%lf\n",w1,w2);




v1=(2*w1*exp(2*w1)+2*(w1*w1+w2*w2)*exp(2*w1)-2*exp(w1)*(w1*z1*cos(w2)+w2*z2*cos(w2)-w2*z1*sin(w2)+w1*z2*sin(w2)+z1*cos(w2)+z2*sin(w2)));
	
	
v2=(2*w2*exp(2*w1)-2*exp(w1)*(-w1*z1*sin(w2)+z2*cos(w2)-w2*z2*sin(w2)-z1*sin(w2)-w2*z1*cos(w2)+w1*z2*cos(w2)));


norm=sqrt(v1*v1+v2*v2);
v1=v1/norm;
v2=v2/norm;


//线性搜索
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



void svd2()
{
	
	
	



    int matrix_order = LAPACK_COL_MAJOR;
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

shuchud(A,s,s);


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
	//r个y   s个K Y m-tau h分割   n步迭代
	
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
	
	
	
	
/*******    初始化   **Yn,**Kn,**yn,**Kn_m,**Yn_m   **********/
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
		
		
		
		
	
//减少整体平移故作位置变换在原数组上覆盖数据	
			
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
	
	////////////////////beta C逆/
	//【Kn1,Kn2,....,Kns】C=[beta,beta,...,beta]
	//【Kn1,Kn2,....,Kns】=[beta,beta,...,beta]*(C^-1)
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
	//yn特殊最左边为当前时刻
	
		for(int j=0;j<s;j++)
		{	
	
	
	//Kni已经更新过
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
	//已更新
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

//为计算方便 特征值互不相同
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
	//返回特征值向量
	
	
	 
	int matrix_order = LAPACK_ROW_MAJOR;
	
	
	LAPACKE_zgeev(matrix_order ,jobvl, jobvr, n, A, lda,w, u , ldvl , v, ldvr);
	
	// printf("--\n");
	 // //shuchuz(w,1,dim);//Av=vD
	// printf("--\n");
	
	
int ipiv[dim];
int info;

	
//左右特征值矩阵的值应该是相同的 至少可逆性一致 考察右特征值是否可逆

LAPACKE_zgeev(matrix_order ,jobvl, jobvr, n, u, lda,tem2, tem1 , ldvl , A, ldvr);
	
	

		
		
memcpy(A,v,dim*dim*sizeof(double _Complex));		
		

	for(int i=0;i<dim;i++)
	{
		if(cabs(tem2[i])<1e-14)
			return 1;
	
		
			
	//不可逆失败 
	
	
	
	
	
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
	//返回特征值向量
	
	
	 
	int matrix_order = LAPACK_ROW_MAJOR;
	
	
	LAPACKE_zgeev(matrix_order ,jobvl, jobvr, n, A, lda,w, u , ldvl , v, ldvr);
	
	
int ipiv[dim];
int info;

	
//左右特征值矩阵的值应该是相同的 至少可逆性一致 考察右特征值是否可逆

LAPACKE_zgeev(matrix_order ,jobvl, jobvr, n, u, lda,tem2, tem1 , ldvl , A, ldvr);
	
	

		
		
memcpy(A,v,dim*dim*sizeof(double _Complex));		
		

	for(int i=0;i<dim;i++)
	{
		if(cabs(tem2[i])<1e-14)
				return 1;
	//不可逆失败
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

/*******傅里叶延时辨识**/
//傅里叶基函数列向量2N+1
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

//int mysqp()


/*
int erci(
		double *H,		//hessian
		double *gk,		//grad
		double *be,		//b   等式约束
		double *ae,		//系数
		double *bi,		////b   不等式约束
		double *ai,
		int dim,		//问题的维数
		int e,			//等式维数
		int ie)
{
//等式约束

//[G,-A;-A' 0] A等式


double *lag=(double *)malloc((dim+e)*(dim+e)*sizeof(double));
double *you=(double *)malloc((dim+e)*sizeof(double));




memset(lag,0,(dim+e)*(dim+e)*sizeof(double));

for(int i=0;i<dim;i++)
{for(int j=0;j<dim;j++)
{lag[(i)*(dim+e)+j]=G[(i)*(dim)+j];}}
//lag赋值G

for(int i=0;i<dim;i++)
{for(int j=0;j<e;j++)
{lag[(i)*(dim+e)+j+dim]=-A[(i)*(e)+j];}}
//lag赋值G
//高立  p217
	
for(int i=0;i<e;i++)
{for(int j=0;j<dim;j++)
{lag[(i+dim)*(dim+e)+j]=-A[(j)*(e)+i];}}


for(int j=0;j<dim;j++)
{you[j]=-gk[j];}

for(int j=0;j<e;j++)
{you[j]=-be[j];}



}

*/

int mychol(
		double *L,		//输入 返回 下三角
		int dim,		//维数
		)
{

int s=dim;
int ipiv[s];
int info;
info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,s,s,L,s,ipiv);
//上三角
//抹去上面  对角线开根号

printf("\n\n\n");
for(int i=0;i<dim;i++)
{L[(i)*dim+i]=sqrt(L[(i)*dim+i]);
for(int j=i+1;j<dim;j++)
{
L[(i)*dim+j]=0;	
}}



	
//得下三角

	
}


int myqp(
		double *G,		//hessian
		double *A,		//grad
		double *gk,
		double *b,
		int dim,		//G维数
		int   e			//A的列数
		)
{
//复制G
double *GC=(double *)malloc(dim*dim*sizeof(double));
double *L=(double *)malloc(dim*dim*sizeof(double));
double *LW=(double *)malloc(e*e*sizeof(double));


double *TEM=(double *)malloc(dim*dim*sizeof(double));
double *GI=(double *)malloc(dim*dim*sizeof(double));
double *TEM2=(double *)malloc(dim*e*sizeof(double));

double *TEM3=(double *)malloc(e*e*sizeof(double));

double *V=(double *)malloc(e*e*sizeof(double));


memcpy(L, G, dim*dim*sizeof(double));	
mychol(L,dim);
//求下三角	
	
memcpy(GI, G, dim*dim*sizeof(double));	
ni(GI,dim)
//求逆


cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, dim, e,dim, 1,GI, dim,1,A, e,0,TEM2,e );	
cblas_dgemm(CblasRowMajor, CblasTrans,CblasNoTrans, e, dim,dim, 1,A, e,1,TEM2, e,0,TEM3,e );	
memcpy(V, TEM3, e*e*sizeof(double));
memcpy(LW, TEM3, e*e*sizeof(double));


mychol(LW,e);
	
	
	

double *u=(double *)malloc(dim*sizeof(double));
	
cblas_daxpby(dim, -1, gk, 1, 0, u, 1);

int info;
	LAPACKE_dgesv(LAPACK_ROW_MAJOR,dim,1,L,dim,ipiv,u,1);

//dim   yigelie
	
LAPACKE_dgesv(LAPACK_COL_MAJOR,dim,1,L,dim,ipiv,u,1);





double *bw=(double *)malloc(e*sizeof(double));

cblas_dgemm(CblasRowMajor, CblasTrans,CblasNoTrans, e, 1,dim, -1,A, e,u, e,1,bw,1 );	


LAPACKE_dgesv(LAPACK_ROW_MAJOR,e,1,LW,e,ipiv,bw,1);

LAPACKE_dgesv(LAPACK_COL_MAJOR,e,1,LW,e,ipiv,bw,1);

cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, dim, 1,e, 1,A, e,bw, e,1,gk,1 );	



LAPACKE_dgesv(LAPACK_ROW_MAJOR,e,1,LW,e,ipiv,bw,1);

LAPACKE_dgesv(LAPACK_COL_MAJOR,e,1,LW,e,ipiv,bw,1);

cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, dim, 1,e, 1,A, e,bw, e,1,gk,1 );	


	
	
		
	
	
	
	
	
	
return 0;	
}










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




 


/*变址计算，将x(n)码位倒置*/  
int changema(double* x,int size_x)        
{  
  double  temp;  
  unsigned short i=0,j=0,k=0;  
  double t;  
  for(i=0;i<size_x;i++)  
  {  
shuchud(x,size_x,1);


    k=i;j=0;  
    t=(log(size_x)/log(2));  
  while((t--)>0 )    //利用按位与以及循环实现码位颠倒  
  {  
    j=j<<1;  
    j|=(k & 1);  
    k=k>>1;  
  }  
  if(j>i)    //将x(n)的码位互换  
  {  
  printf("-----------\n");
  shuchud(x,size_x,1);
    temp=x[i]; 
	 printf("%lf,\n",temp);
 
    x[i]=x[j];  
printf("%lf,\n",temp);
    x[j]=temp;
printf("####%lf,%lf####%d",temp,x[3],j);


shuchud(x,size_x,1);	
  printf("++++++%d++\n",i);
   } 
   
   
   
  }  
  
  
  
   shuchud(x,size_x,1);
   return 1;

  
}  






















// for(int i=0;i<s;i++)
// {shuchud(Kn_m[i],dim,1);	
// printf("\n");
// }
	








// for(i=0;i<m;i++)  
    // {
        // for(j=0;j<n;j++)  
        // {
            // printf("%p\n",&a[i][j]);     //输出每个元素地址，每行的列与列之间的地址时连续的，行与行之间的地址不连续
        // }
    // }
    // for(i=0;i<m;i++)  
    // free(a[i]);
 
    // free(a);  
// --------------------- 
// 作者：阿阿阿阿阿阿鑫 
// 来源：CSDN 
// 原文：https://blog.csdn.net/fengxinlinux/article/details/51541003 
// 版权声明：本文为博主原创文章，转载请附上博文链接！








	// OO=(double *)malloc(n*n*sizeof(double));
	// memset(OO,0,n*n*sizeof(double));



