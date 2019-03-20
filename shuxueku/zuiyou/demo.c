

//****** 傅里叶变换  */
double _Complex xxx[4]={1+0*I,2+0*I,-1+0*I,3+0*I},ffo[4]={1+0*I};


fft(xxx,ffo,4);
shuchuz(ffo,4,1);	

/* 复数运算 */

double complex a=1+8*I,b=1-8*I;	
printf("%1.1f,%1.1f\n\n\n",cimag(a+b),creal(a*b));

