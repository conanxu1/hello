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
double be[2]={-7,-4};
double Ae[4]={-4,-5,-5,9};

double bi[2]={-2,-2};

double Ai[4]={-1,-8,1,-5};


double *lam;
lam=cshi(1);
double x0[2]={0,0};




qxt(H,	h,be ,Ae,	bi,	Ai,2,	2,2,x0,lam);

shuchud(x0,2,1);


 








return 0;
}
 
 
 
 

