#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <cblas.h>
#include "myfun.h"
#include "ofun.h"


void matrix_vector();


float y[2] = {};


void main()
{
	
	double x[2]={0,1},y[2]={0,1};
	
	printf("\n\nhello%f\n",jifen(0,3,100,f1));
	//matrix_vector();
	
	
	
	
	
}






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
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 