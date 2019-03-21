#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stddef.h>
#include <malloc.h>
#include <string.h>
#include "myfun.h"

#include <complex.h>









// extern lapack_int LAPACKE_dgeev( int matrix_order, char jobvl, char jobvr,
                          // lapack_int n, double* a, lapack_int lda, double* wr,
                          // double* wi, double* vl, lapack_int ldvl, double* vr,
                          // lapack_int ldvr );


void matrix_vector();



int main()
{	
int d=3;

//注意是否有1/2*G
double H[4]={0.5,0,0,0.5};
				
double h[2]={-3,-1};


//注意约束的写法ATx=b
//改为Ax=b 可以有效增减约束
 
double Ae[0];
double be[0];



double Ai[10]={ 0,-1,
				-1,-1,
				-1,0,
				1,0,
				0,1};
				
double bi[5]={-3,-4,-2,0,0};
double xk[2]={0,3};
				
				
				
erci(H,h,be,Ae,Bi,Ai,2,0,5,);





return 0;
}
 
 
 
 
 
 
 
