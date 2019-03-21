

//****** 傅里叶变换  */
double _Complex xxx[4]={1+0*I,2+0*I,-1+0*I,3+0*I},ffo[4]={1+0*I};

//////////////////////////////////




fft(xxx,ffo,4);
shuchuz(ffo,4,1);	

/****** 复数运算 */

double complex a=1+8*I,b=1-8*I;	
printf("%1.1f,%1.1f\n\n\n",cimag(a+b),creal(a*b));
//////////////////////////////////



/********** 微分方程 */
double x0[3]={1,1,1};

double jieguo[(dim+1)*(n+1)];
//多一维时间 多一个初始值
//格式 n分量和时间

y=ode1(0,10,10000,lianxu,x0,3,jieguo);
//设定初始终止值 划分区间数 右端函数句柄  初始值 和维数
/////////////////////////////////////////




/*****     显示时间      */
time_t timer;//time_t就是long int 类型
struct tm *tblock;
timer = time(NULL);
tblock = localtime(&timer);
printf("Local time is: %s\n", asctime(tblock));
///////////////////////////




