#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stddef.h>
#include <malloc.h>
#include <string.h>
#include "myfun.h"

#include <complex.h>



 

int main()
{

double H[4]={5,-0.5,-0.5,1};

double h[2]={-2,-6};
double be[1]={-7};
double Ae[2]={-4,-5};

double bi[2]={2,2};

double Ai[4]={-1,0,1,0};


double *lam;
lam=cshi(1);
double x0[2]={0,0};




qxt(H,	h,be ,Ae,	bi,	Ai,2,	1,2,x0,lam);

shuchud(x0,2,1);
shuchud(lam,1,3);




/*
double L[4]={16.952368,-23.416053,
-23.416053,34.636579};
double Lz[4]={16.952368,-23.416053,
-23.416053,34.636579};

printf("\n\n--1-\n\n");

shuchud(L,2,2);

printf("\n\n-2--\n\n");
mychol(L,2,Lz);
printf("\n\n-3--\n\n");
shuchud(L,2,2);
shuchud(Lz,2,2);








double H[9]={5,-0.5,0,
	-0.5,1,0,
	0 ,0,100};

double h[3]={-1.853362,
-4.758316,
1674.134492,
};



double be[2]={0,0};
double Ae[6]={4.000000,5.000000,1.000000,
-1.000000,-8.000000,1.000000,};

double bi[0]={};

double Ai[0]={};


double *lam;
lam=cshi(1);
double x0[3]={0,0,0};




myqp(H,	Ae,h,	be,3,2);

shuchud(h,3,1);


 

*/






return 0;
}
 
 
 


/*


 

*/
 

