#include <math.h>
#include <stdlib.h>


#ifdef  AAA

#include <openblas/cblas.h>

#else
#include <cblas.h>


#endif





double f1(double x)
{
	
	return sin(x);
	
	
}


double *lianxu(double t,double *x)
{	

//可在此加入控制

int size=3;
double *dy =(double *)malloc(size*sizeof(double));
	dy[0]=x[0];
	dy[1]=x[1];
	dy[2]=x[2];
	return dy;
}




double gk(double *x,int n)
{
	int i;
	double s=0;
	for(i=0;i<n;i++)
	{
		s+=x[i]*x[i];
		
	}
	
	return s;
}













/*
double *kongzhi(double t)
{
dim=2;
	
double *u =(double *)malloc(dim*sizeof(double));

u[0]=0;
u[1]=sigmak[kk_]*shixing(t,TTk[kk_-1],TTk[kk_]);
	
return u;	
}



/////////余长君ppt p27  f(x,sigma) [tk-1,tk)
//带控制

double *exmf(double t,double *x)
{	
dim=2
	
double *dy =(double *)malloc(dim*sizeof(double));
double *u =(double *)malloc(dim*sizeof(double));
u=kongzhi(t);

dy[0]=x[1];
dy[1]=-x[0]+x[1]*(1.4-0.14*x[1]*x[1])+4*u[1];



return dy;
}



int *pfpx(double t,double *x,double *u,double *A)
{
dim=2;



A[0]=0;
A[1]=1;
A[2]=-1;
A[3]=1.4-3*0.14*x[1]*x[1];

return 0;
}

int *pfpu(double t,double *x,double *u,double *A)
{
dim=2;



A[0]=0;
A[1]=0;
A[2]=0;
A[3]=4;

return 0;
}



//余长君ppt p25
//D ga=P f/P x *ga +P f/P u 
//gamak相当于x而    原x由方程已知 gamak为矩阵分dim次求


double *gamak(double t,double *x)
{int dim;
double *A;
double *Tk;  //确定分点

double *tem=(double *)malloc(dim*dim*sizeof(double));

double *FPX=(double *)malloc(dim*dim*sizeof(double));
double *FPU=(double *)malloc(dim*dim*sizeof(double))
double *dy =(double *)malloc(dim*sizeof(double));
//dy即gamak分一列求




pfpx(t,x,u,FPX);
pfpu(t,x,u,FPU);

//pfpx(t,A)* gamak+pfpu


cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasTrans, dim, dim,dim,1,FPX,dim, x,1, 1,FPU, dim);
	
for(int i=0;i<dim;i++)
dy[i]=FPU[(i-1)*dim];


return dy;
}















*/










double* fini(double t)
{	int dim=2;
	double *y =(double *)malloc(dim*sizeof(double));
	y[0]=sin(t)-2;
	y[1]=t+2;


	
	
	
    return y;
	
}

double* dfini(double t)
{	int dim=2;
	double *dy =(double *)malloc(dim*sizeof(double));
	dy[0]=cos(t);
	dy[1]=1;

	
    return dy;
	
}













/*


#include<stdio.h>
#define N 10
int main(){
    int *addOne(int a[]);
    int a[N]={1,2,3,4,5,6,7,8,9,10};
    int i;
    printf("\n调用函数之后:\n");
    int *b = addOne(a);
    for(i=0;i<N;i++)
        printf("%d\t",b[i]);
    
}
int *addOne(int a[]){
    int i;
    for(i=0;i<N;i++)//全部加一
        a[i] += 1;
    return a;
}








int *addOne(int a[]){
    int *b = (int *)malloc(N*sizeof(int));//定一个int型的指针b，并申请N*sizeof（int）个字节的存储空间，即N*4个字节
    int i;
    for(i=0;i<N;i++)//全部加一
        b[i] = a[i] + 1;
    return b;
}



#include<stdio.h>
#define M 3
#define N 2
int main(){
    int **addOne(int a[M][N]);
    int a[M][N]={{1,1},{2,2},{3,3}};
    int i,j;
    printf("\n调用函数之后:\n");
    int **b = addOne((int **)a);
    for(i=0;i<M;i++)
        for(j=0;j<N;j++)
            printf("%d\t",*((int *)b+N*i+j));
    
}
int **addOne(int a[M][N]){
    int **b =(int **)malloc(M*sizeof(int *));//先申请M个指针型字节的空间
    for (int i=0;i<M;i++)
        b[i]=(int *)malloc(N*sizeof(int));//然后依次按一维申请
    int i,j;
    for(i=0;i<M;i++)
        for(j=0;j<N;j++)
            *((int *)b+N*i+j) = *((int *)a+N*i+j) + 1;//找地址，如同数据结构中矩阵找地址相同，首地址+（次数行数-1）*总列数+次数列数-1
                                                    //i，j都是从0开始，可以不用减1
    return b;
}







*/














