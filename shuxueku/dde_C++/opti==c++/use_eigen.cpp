#include "pro.h"

#ifdef USE_EIGEN
void double2mat(double *d_a,Rmn &m_a,int m,int n)
{
	for(int i=0;i<m;i++)
	for(int j=0;j<n;j++)
		m_a(i,j)=d_a[i*n+j];
}



void mat2double(Rmn &m_a,double *d_a,int m,int n)
{
	for(int i=0;i<m;i++)
	for(int j=0;j<n;j++)
		d_a[i*n+j]=m_a(i,j);
}


void double2vex(double *d_a,Rn &m_a,int m)
{
	for(int i=0;i<m;i++)
	m_a[i]=d_a[i];
}




void vec2double(Rn &m_a,double *d_a,int m)
{
	for(int i=0;i<m;i++)
	d_a[i]=m_a[i];
}



#endif
