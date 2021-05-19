#ifndef _PRO_H_
#define _PRO_H_
	#include "pro.h"
#endif

#ifndef _EX1_H_
#define _EX1_H_
#include "ex1.h"
#endif
 

int main()
{

double ttA[2*2]={1,2,3,4};
double ttB[2*1]={-2,1};
double tx0[2]={1,1};


Rmn A=Rmninit(dimx,dimx);
Rmn B=Rmninit(dimx,dimu);
Rn x0=Rninit(dimx);

Rmn_copy(ttA,A,dimx,dimx);
Rmn_copy(ttB,B,dimx,dimu);
Rn_copy(tx0,x0,dimx);




 
  



Eu_Lode_Sol hSol(dimx,dimu,t0,tf,700000000);
hSol.set(A,B,myut,x0 );



hSol.sol( );


//test();







//<-------------------------------
		
		// // // double tol=1e-8;
		
		// // // Pro mypro(dimx,dimu,0,1,1,1);
		
		// // // mypro.set_daef(my_fxut);
		// // // mypro.set_opt_int(my_opt_L);
		// // // mypro.set_opt_phi(my_opt_phi);
		// // // mypro.set_opt_pLp(	 pLpx,	  pLpu);

	 
		
		// // // // // // // // // // //cout<<mypro.t0;
		
		
		// // // Sol mysol(&mypro);
		// // // mysol.set();
		// // // mysol.display()		;
		// // // mysol.sol(N);
		
	
		

		// // // Rn th=Rninit(dimx*(N+1)+dimu*(N+1)+2);
		// // // for(int i=0;i<(dimx*(N+1)+dimu*(N+1)+2);i++)
		// // // {
			// // // th[i]=i;	
		// // // }
		

			
	// // // mysol.hk_th(th,0);
	
  
	// // // int	NN=3;
	// // // ga_opti myga(dimx,NN,10);
	
	// // // Rn lb =Rninit(dimx);
	// // // Rn ub =Rninit(dimx);
	// // // for(int i=0;i<dimx;i++)
	// // // {
		// // // lb[i]=-1;
		// // // ub[i]=1;
	// // // }

	// // // myga.set_lbub(lb,ub);
	// // // myga.initialize();
	// // // myga.chooose();

//------------------------->



	
	return 1;
}



// void fk(Rn a)
// {
	// a[0]=1;
// }	
	



/*
MatrixXd* CC;
MatrixXd  AA(5,5);

double BB[4]={1,1,1,1};

cout<<"\n------>\n"<<sizeof(AA);



cout<<"\n======>\n"<<sizeof(BB);
cout<<"\n===CC=>\n"<<sizeof(CC);*/
