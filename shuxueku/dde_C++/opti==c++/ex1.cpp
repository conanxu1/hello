#include "ex1.h"



Rn my_fxut(Rn x ,      Rn u     ,double        t )
{
		
	Rn f;
	f=Rninit(2);

	f[0]=x[1];
	f[1]=u[0];
	
	return f;
}

R my_opt_phi(Rn xt0 ,     double t0, Rn xtf ,     double tf )
{
	return 0;
}


R my_opt_L(Rn x ,      Rn u     ,double        t )
{
	return u[0]*u[0];
}

Rn pLpx(Rn x ,      Rn u     ,double        t,int dimx,int dimu )
{
	Rn p=Rninit(dimx);
	

	return p;
}

Rn pLpu(Rn x ,      Rn u     ,double        t ,int dimx,int dimu)
{
	Rn p=Rninit(dimu);
	for(int i=0;i<dimu;i++)
	{
		p[i]=2*u[i];
	}
	
	return p;
 
}





// // // // R myf(Rn x ,      Rn u     ,double        t )
// // // // {
	// // // // return x[0]*x[0];
// // // // }











Rn myut(double t)
{
	temut[0]=sin(t);
	
	return temut;
}
		