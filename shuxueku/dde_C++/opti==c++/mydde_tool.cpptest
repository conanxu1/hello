#include "ddetool.h"


using namespace std::chrono;


Rn myftx(double t,Rn x ,    Rn* xdelay         ,C_R fu)
{
	Rn dx;
	dx=(Rn)malloc( 1*sizeof(double ));
	dx[0]=x[0]  * xdelay[0][0] ;
	
	
	// dx[1]=-x[1]  -xdelay[0][2]  ; 
	
	// Rn d2;
	// d2=(T->*fx)(t);
	// cout<<d2[0];
	
	return dx;
}



Rn   myphi(double t0,double t   )
{
	Rn y;
	y=(Rn )malloc(1*sizeof(double));
	
	if( fabs(t-t0)<1e-10)
	{y[0]=1;
	}
	else{
		y[0]=0;
	}
	// y[1]=1;
	
	return y;
}


VectorXd   Vec_myphi(double t0,double t   )
{
	VectorXd y;
 
	if( fabs(t-t0)<1e-7)
	{
		y = VectorXd::Ones(2,1);
	}
	else{
		y=VectorXd::Zero(2,1);
	}
 
	
	return y;
}






Rn   myut(double t   )
{
	Rn y;
	y=(Rn )malloc(sizeof(double));
	
	y[0]=0;

	
	return y;
}


VectorXd   Vec_myut(double t   )
{
	VectorXd y;
	y=   MatrixXd::Zero(1,1);
	
	y(0,0)=sin(t);
	
	return y;
}




int main_haha()
{
	int dimx=2;
	int dimu=1;
	
	int tauN=1;
	
	
	RK4_ddesolver solver1(1,1);
	
	double tauk[1]={0.5};
	
	solver1.setTauN(tauN);
	
	solver1.setTauk(tauk);
	solver1.set_t0_tf(0,50);
	solver1.set_phi(myphi);
	solver1.set_ut(myut);
	
	solver1.set_eq(myftx);
	
	solver1.solve();
	
	
	
	
	
	solver1.display();
	solver1.display_sol();
	
	VectorXd x;
	VectorXd xdelay[tauN];
	// // // solver1.fAiut(0 , x ,   xdelay    , Vec_myut);
	
	Matrix<double, Dynamic, Dynamic> A0,A1,B0,B1;
	
	Matrix<double, Dynamic, Dynamic> tAi[tauN+1],tBi[tauN+1];
	

	
	
	A0=MatrixXd::Zero(dimx,dimx);
	A1=MatrixXd::Zero(dimx,dimx);

	B0=MatrixXd::Zero(dimx,dimu);
	B1=MatrixXd::Zero(dimx,dimu);
	
	 A0<<-2,0,
				 0,-1;
	 A1<<0,0,
				 0,0;

	B0<<	1,
				0;
	B1<<	0,
				0;
	
	
	
	tAi[0]=A0;
	tAi[1]=A1;
	
	tBi[0]=B0;
	tBi[1]=B1;
 
	
	solver1.set_Ai(tAi);
	solver1.set_Bi(tBi);
	
	solver1.set_vecphi(Vec_myphi);
	solver1.set_vecut(Vec_myut);
	solver1.display_AiBi();
	
	// auto t1=high_resolution_clock::now();

	solver1.solveAiBi();

	// auto t2=high_resolution_clock::now();	
	// auto duration = duration_cast<microseconds>(t2 - t1);
	// cout<<double(duration.count()) * microseconds::period::num / microseconds::period::den   << "秒" << endl;

	solver1.display_solVec();
	
	return 1;
	
	
	
}