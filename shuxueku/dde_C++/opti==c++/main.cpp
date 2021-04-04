#ifndef _PRO_H_
	#include "pro.h"
#endif
 




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

Sol* pmysol;
int N=2;
int dimx=2;
int dimu=1;
	
	
	
	
// void fk(Rn a)
// {
	// a[0]=1;
// }	
	
	
int main()
{

		
		double tol=1e-8;
		
		Pro mypro(dimx,dimu,0,1,1,1);
		
		mypro.set_daef(my_fxut);
		mypro.set_opt_int(my_opt_L);
		mypro.set_opt_phi(my_opt_phi);
		mypro.set_opt_pLp(	 pLpx,	  pLpu);

	 
		
		//cout<<mypro.t0;
		
		
		Sol mysol(&mypro);
		mysol.set();
		mysol.display()		;
		mysol.sol(N);
		
	
		

		Rn th=Rninit(dimx*(N+1)+dimu*(N+1)+2);
		for(int i=0;i<(dimx*(N+1)+dimu*(N+1)+2);i++)
		{
			th[i]=i;	
		}
		

			
	mysol.hk_th(th,0);
	
  
	int	NN=3;
	ga_opti myga(dimx,NN,10);
	
	Rn lb =Rninit(dimx);
	Rn ub =Rninit(dimx);
	for(int i=0;i<dimx;i++)
	{
		lb[i]=-1;
		ub[i]=1;
	}

	myga.set_lbub(lb,ub);
	myga.initialize();
	myga.chooose();

		
		



	
	return 1;
}








