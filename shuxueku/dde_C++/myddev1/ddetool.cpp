#include "ddetool.h"

void RK4_ddesolver :: display()
{
	cout<<"dimx:"<<this->dimx<<endl;
	cout<<"dimu:"<<this->dimu<<endl;


	cout<< "tauk"<<": ";
	for(int i=0; i<(this->tauN);i++)
	{
		cout<<(this->tauk)[i]<<", ";
	}
	cout <<endl;
	cout<<"taumax:"<<taumax<<endl;
	
	cout<< "[t0,tf]: "<<(this->t0)<<" , "<<(this->tf)<<endl;
	
	cout<< "phi"<<": "<<endl;
	
	
	
	
	
 }


void RK4_ddesolver :: display_AiBi()
{
	cout<<"Ai:"<<endl;
	for(int i=0; i<(this->tauN  + 1 );i++)
	{
		cout<<i<<"==="<<endl<<(this->Ai)[i]<<endl<<endl;
	}
	cout<<"Bi:"<<endl;
	for(int i=0; i<(this->tauN  + 1 );i++)
	{
		cout<<i<<"==="<<endl<<(this->Bi)[i]<<endl<<endl;
	}
	
}


void RK4_ddesolver :: display_sol()
{


	// for(int i=0;i<(this->stepN);i++)
	// {
	// cout<<setw(20)<<setprecision(10)<<(this->solx)[i][0]<<","<<endl;
	// }
	
	cout<<setw(20)<<setprecision(10)<<(this->solx)[this->stepN - 1  ][ 0  ]<<","<<endl;

}




void RK4_ddesolver :: display_solVec()
{
	// for(int i=0;i<(this->stepN);i++)
	// {
	// cout<<setw(20)<<setprecision(10)<<(this->solx)[i][0]<<","<<endl;
	// }
	
	cout<<setw(20)<<setprecision(10)<<(this->Vec_solx)[this->stepN - 1  ][ 0  ]<<","<<endl;

}





void RK4_ddesolver :: setTauN(int tauN)
{
this->tauN=tauN;
this->tauk=(Rn)malloc((this->tauN)*sizeof(double)) ;
}

void RK4_ddesolver :: setTauk(Rn tauk)
{
memcpy(this->tauk,tauk,(this->tauN)*sizeof(double));
this->taumax=tauk[this->tauN -1 ];
}

void RK4_ddesolver ::  set_t0_tf(double t0,double tf)
{
this->t0=t0;
this->tf=tf;


this->h0=(this->tf-this->t0)/(this->stepN   -  1);
}


void RK4_ddesolver ::  set_phi(Cphi mphi)
{
this->phi=mphi;
}

void RK4_ddesolver ::  set_ut(C_R mut)
{
this->ut=mut;
}


void RK4_ddesolver ::  set_vecphi(Vec_Cphi mphi)
{
this->vec_phi=mphi;
}

void RK4_ddesolver ::  set_vecut(Vec_C_R mut)
{
this->vec_u=mut;
}	 
 








void RK4_ddesolver ::  set_eq(functional   mFtx )
{
this->Ftx=mFtx;
}



void RK4_ddesolver :: set_Ai( Matrix<double, Dynamic, Dynamic> *tAi  )
{

this->Ai=tAi;
}

void RK4_ddesolver :: set_Bi( Matrix<double, Dynamic, Dynamic>* tBi  )
{
this->Bi=tBi;
}



void RK4_ddesolver ::  onestep(Rn xn1,Rn xn ,double hh,Rn ff   )
{

for(int i=0;i<(this->dimx);i++)
{
	xn1[i]=xn[i]+hh*ff[i];	
}


}




void RK4_ddesolver :: 	calx( C_R f,double t,Rn as        )
{
Rn tx;
tx=f(t);
memcpy( as,tx,(this->dimx)*sizeof(double)   );
delete tx;
}

void RK4_ddesolver :: 	calu( C_R f,double t,Rn as        )
{
Rn tx;
tx=f(t);
memcpy( as,tx,(this->dimu)*sizeof(double)   );
delete tx;
}



Rn* RK4_ddesolver ::  chi_t(double t)
{
	
	Rn* y;
	y=(Rn*)malloc(this->tauN  * sizeof(  Rn  ) );
	
	
	for(int i=0;i<tauN;i++)
	{
		
		
		//有点问题  t-tauk[i]
		if(      (t-tauk[i])   <= t0   )
		{
			y[i]=(this->phi)( this->t0, t - (this->tauk)[i]  );	
		}
		else
		{
			y[i]=(this->solx)[(int)round( ( t - this->t0   - (this->tauk)[i] )   /( this->h0  ) )];
		}
	}
	
	
	
	return y;
	
}




VectorXd* RK4_ddesolver ::  vec_chi_t(double t)
{
 
	
	VectorXd* y=new VectorXd[this->tauN];

	
	for(int i=0;i< (this->tauN);i++)
	{
		if(            ( t - (this->tauk)[i]  )          <= t0  )
		{
			y[i]=(this->vec_phi)( this->t0, t - (this->tauk)[i]  );	
		}
		else
		{
			y[i]=(this->Vec_solx )[(int)round( ( t - this->t0   - (this->tauk)[i] )   /( this->h0  ) )];
		}
	}
	
	
	// cout<<"+++"<<this->tauN<<";;;;["<<sizeof(y)<<"]";
	return y;
	
}






 
void RK4_ddesolver ::  solve( )
{
double h0;
h0=this->h0;
Rn tphi;
int Nmax;
Nmax=this->stepN;



(this->solt)=(Rn)malloc(Nmax*sizeof(double ));
(this->solx)=(Rn*)malloc(Nmax*sizeof(Rn));
(this->solx)[0]=(Rn)malloc(  (this->dimx)*sizeof(double));



tphi= (this->phi)(this->t0 , this->t0 );
memcpy( (this->solx)[0], tphi,(this->dimx)*sizeof(double)  );
delete tphi;

(this->solt)[0]=this->t0;

double tnow;
tnow=(this->t0);



Rn* tchi_t;
Rn tvec1;

int itr=0;
while(itr<Nmax)
{
	(this->solx)[itr+1]=(Rn)malloc(  (this->dimx)*sizeof(double));
	(this->solt)[itr+1]=tnow;


	tchi_t=chi_t(tnow);

	tvec1=(this->Ftx)(tnow , (this->solx)[itr]   ,   tchi_t  ,this->ut );
	onestep(  (this->solx)[itr+1],(this->solx)[itr] ,h0,tvec1                      );
	tnow=tnow+h0;
	itr++;


	delete tvec1;
	for(int i=0;i<tauN;i++)
	{
		delete tchi_t[i];
	}

}

}


void RK4_ddesolver ::  solveAiBi( )
{
	double h0;
	h0=this->h0;
	
	int Nmax;
	Nmax=this->stepN;
 
	cout<<endl<<Nmax<<endl;
	
	this->Vec_solx=new VectorXd[Nmax];
	this->solt=(Rn)malloc(Nmax*sizeof(double));
	 
	
	double tnow;
	tnow=this->t0;
	
	this->solt[0]=tnow;
	this->Vec_solx[0]=(this->vec_phi)(tnow,tnow);



	VectorXd *vec_delay;
	VectorXd tmp;
	
	int itr=1;
	while(itr<Nmax)
	{
		vec_delay=vec_chi_t(  tnow  );
		 
		
	 // tmp=fAiut( tnow , (this->Vec_solx)[itr-1] ,  vec_delay     , this->vec_u);
	 // cout<<tmp <<",,"<< itr/this->tf*h0 <<endl<<endl;
	 
	 
		(this->Vec_solx)[itr]=		(this->Vec_solx)[itr-1]+h0*fAiut( tnow , (this->Vec_solx)[itr-1] ,  vec_delay     , this->vec_u);
		
 
		tnow=tnow+h0;
		itr++;
 
	}
	
	
	


}





VectorXd RK4_ddesolver ::  fAiut(double t , VectorXd x ,   VectorXd *xdelay    , Vec_C_R fu)
{
	

	
	VectorXd  dxdt ;
	
	dxdt=(this->Ai)[0]*x+(this->Bi)[0]*fu(t);

	 
	 
	 for(int i=1;i<= (this->tauN);i++)
	 {
		 dxdt=dxdt+(this->Ai)[i]*xdelay[i-1]+(this->Bi)[i]*fu(t  - this->tauk[i-1]);
		
		
	}
	
 
	 

	return dxdt;
	
}



VectorXd RK4_ddesolver ::  f2(double t  ,VectorXd x ,VectorXd* xdelay, Vec_C_R fu)
{
	VectorXd  dxdt ;
	
	
	
	dxdt=x;
	
	
	
	return dxdt;
}




