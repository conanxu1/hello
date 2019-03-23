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

/*	
int d=3;

//注意是否有1/2*G

double H[4]={1,0,0,1};
				
double h[2]={0,0};


//注意约束的写法ATx=b
//改为Ax=b 可以有效增减约束
 
double Ae[1*2]={1,2};
double be[1]={4};



double Ai[2]={2,2};
				
double bi[1]={-1};

qxt(H,h,be,Ae,bi,Ai,2,1,1);


*/

int dim=3;

//注意是否有1/2*G

double H[9]={  1,0,0,
		0,1,0,
		0,0,160};

double zuoyong[3*3]={   1,2,1,
			-1,-2,1,
			0,0,1};
double tb[3]={0,0,0};
double tg[3]={1.6,3.2,80};


double *AM;

//xxwg(AM,zuoyong,3,3);

myqp(H,zuoyong,tg,tb,dim,3);


//shuchud(tg,dim,1);



return 0;
}
 
 
 
 
 
 
 
