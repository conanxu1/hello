#include<stdio.h>
#include<math.h>

double f(double x)
{return pow(x-0.2,2)+exp(x);}



int main()
{
	double ep=1e-8;
	double tau=(sqrt(5)-1)/2;
	double a=-10;
	double b=10;
	double alpha=0.5*(b-a);
	double al=a;
	double ar=b;

	while(b-a>ep)
	{
		al=a+(1-tau)*(b-a);
		ar=a+tau*(b-a);
		alpha=0.5*(a+b);
		if(f(al)<f(ar))
		{
			b=ar;
			ar=al;
			al=a+(1-tau)*(b-a);
		}
		else
		{
			a=al;
			al=ar;
			ar=a+tau*(b-a);
		}






		printf("计算中a=%lf,b=%lf\n",a,b);
	}
	alpha=0.5*(a+b);
	printf("结果%lf\n",alpha);







	return 0;
}

