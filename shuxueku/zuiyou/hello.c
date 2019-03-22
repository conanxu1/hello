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
double H[9]={1,0,0,0,1,0,0,0,1};
				
double h[3]={0,0,0};


//注意约束的写法ATx=b
//改为Ax=b 可以有效增减约束
 
double Ae[1*3]={1,2,-1};
double be[1]={4};



double Ai[3]={2,2,2};
				
double bi[1]={-1};

qxt(H,h,be,Ae,bi,Ai,3,1,1);





return 0;
}
 
 
 
 
 
 
 
