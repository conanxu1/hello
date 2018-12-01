#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stddef.h>
#include <malloc.h>

#include "myfun.h"

#include <complex.h>









// extern lapack_int LAPACKE_dgeev( int matrix_order, char jobvl, char jobvr,
                          // lapack_int n, double* a, lapack_int lda, double* wr,
                          // double* wi, double* vl, lapack_int ldvl, double* vr,
                          // lapack_int ldvr );


void matrix_vector();



void main()
{
	
	
	/*
	double *y,x0[3]={1,1,1},*p;
	
	
	double y0[3]={1,1,1};
time_t timer;//time_t就是long int 类型
struct tm *tblock;
timer = time(NULL);
tblock = localtime(&timer);
printf("Local time is: %s\n", asctime(tblock));	




	
y=ode1(0,10,10000,lianxu,x0,3);
	
printf("ok");	


 timer = time(NULL);
 tblock = localtime(&timer);
printf("Local time is: %s\n", asctime(tblock));

*/


time_t timer;//time_t就是long int 类型
struct tm *tblock;
timer = time(NULL);
tblock = localtime(&timer);
printf("Local time is: %s\n", asctime(tblock));	



drk();

// svd2();




/*
double* p;
double* C;
double *U,*V;

C=danwei(2,3);

double *q;

int n=1000;

q=duqu("wind.txt",n*n);



printf(";;;;;;;;;\n");
svd(q,n,n);

// house(q,4);
// shuchud(q,4,1);

// shuchud(q,2,2);

// p=danwei(2,3);


// printf("\n");
// p[3]=7;
// shuchud(p,2,3);


// printf("\n");


// cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans, 2, 3,2, 1,q, 2,p,3, 0,C, 3); 


// // cblas_sgemv(CblasRowMajor, CblasNoTrans, m, n, alpha, array, n, x, 1, beta, y, 1); 
// // cblas_sgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,M,N,K,alpha,A,A的列数,B,B的列数,beta,C,C的列数)

// // CblasRowMajor表示数组时行为主，相应矩阵大小为(M*K)乘以(K*N)，可以得到M，K，N的值

// // cblas_sgemm(order,transA,transB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC);



// shuchud(C,2,3);
*/



 timer = time(NULL);
 tblock = localtime(&timer);
printf("Local time is: %s\n", asctime(tblock));


}




/*

double _Complex d;
for(int ii=0;ii<1000;ii++)
d=lamw(100,100);


 printf(" %.15f,i%.15f\n", creal(d),cimag(d));





void matrix_vector(){
 float array[6] = { 1,2,3,4,5,6 };
 float x[3] = { 1,2,3 };
 int m = 2; 
 int n = 3;
 
 float alpha=1,beta=1;
 printf("\n\n\n");
 
 
 shuchuf(array,2,3);
 
 cblas_sgemv(CblasRowMajor, CblasNoTrans, m, n, alpha, array, n, x, 1, beta, y, 1); 
 
	 printf("\n\n\n");
	shuchuf(x,3,1);
	 printf("\n\n\n");
	
	shuchuf(y,2,1);

 //cblas_sgemv(CblasColMajor, CblasTrans, n, m, alpha, array, n, x, 1, beta, y, inc_y); //两个结果一致 
 }
 
 
 
 
 
 
 
 
 
 

 
 double x[2]={0,1},y[2]={0,1};
	
	printf("\n\nhello%f\n",jifen(0,3,100,f1));
	//matrix_vector();
	
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 
 
 
 
 
 
 
 
 
 
 
 	p=y0;
	
	cblas_daxpby(3, 1, x0, 1, 1, y0, 1);
	shuchud(x0,1,3);
	shuchud(y0,1,3);

	shuchud(p,1,3);

	
	
	
	
	
	
	
	memset（a，val，sizeof(a));快速赋值 
	//
 
 */
 
 
 
 
/*extern lapack_int LAPACKE_dgeev( int matrix_order, char jobvl, char jobvr,
                          lapack_int n, double* a, lapack_int lda, double* wr,
                          double* wi, double* vl, lapack_int ldvl, double* vr,
                          lapack_int ldvr );
*/

 
 
 
 
 
 
 
 
 
