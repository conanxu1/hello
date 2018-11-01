#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stddef.h>
#include <malloc.h>

#include "myfun.h"











void matrix_vector();



void main()
{
	
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

}




/*

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

	
	
	
	//
 
 */
 
 
 
 
 
 
 
 
 
 
 
 
 
