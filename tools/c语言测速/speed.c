#include<time.h>
#include<stdio.h>


int main()
{


clock_t start,end;  
start = clock();  

int a=0;
for(int i=0;i<10000;i++)
{
	a=a+i;
}

end = clock();  
printf("time=%f\n",(double)(end-start)/CLOCKS_PER_SEC);


}