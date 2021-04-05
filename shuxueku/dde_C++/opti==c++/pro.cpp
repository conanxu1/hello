#include "pro.h"








#ifdef noeigen
	Rn Rninit(int nnnn)
	{
		double *p;
		p=(double *)malloc(nnnn*sizeof(double));
		memset( p,nnnn,sizeof(double) );
		return p;
	}
	
	Rn Rmninit(int m,int n)
	{
		double *p;
		p=(double *)malloc((m*n)*sizeof(double));
		memset( p,m*n,sizeof(double) );
		return p;
	}
	
	

	void theta2x(Rn theta,Rn x,int k, int dimx    )
	{
		for(int i=0;i<dimx;i++)
		{
			x[i]=theta[ k*dimx+i  ];
		}
	}
	
	void theta2u(Rn theta,Rn u,int k, int dimx,int dimu,int N    )
	{
		for(int i=0;i<dimu;i++)
		{
			u[i]=theta[(N+1)*dimx+ k*dimu+i  ];
		}
	}

	void axbyc(double a,Rn x,double b,Rn y,Rn c,int hdim)
	{
		for(int i=0;i<hdim;i++)
		c[i]=a*x[i]+b*y[i];
	}


	

	


#endif
	
	
#ifndef noeigen	
	VectorXd Rninit(int nnnn)
	{
		
		VectorXd vec;
		vec=   MatrixXd::Zero(nnnn,1);
		return vec;
	}
	
	
	MatrixXd Rmninit(int m,int n)
	{
		MatrixXd mat(m,n);
		mat=MatrixXd::Zero(m,n);
		return mat;
	}
	
	
	
	void theta2x(Rn theta,Rn &x,int k, int dimx    )
	{
		for(int i=0;i<dimx;i++)
		{
			x[i]=theta[ k*dimx+i  ];
		}
	}
 
	void theta2u(Rn theta,Rn &u,int k, int dimx,int dimu,int N    )
	{
		for(int i=0;i<dimu;i++)
		{
			u[i]=theta[(N+1)*dimx+ k*dimu+i  ];
		}
	}
	
	void axbyc(double a,Rn x,double b,Rn y,Rn c,int hdim)
	{
		c=a*x+b*y;
	}
	
	
	
	void aMx_cC(double a,Rmn tM,Rn x   ,double c,Rmn C,int hdim)
	{
		C=a*tM*x+c*C;
		
	}
	
	
	
#endif







#ifdef noeigen	
void Rmn_copy(double* tm,Rmn tM,int m,int n)
{
	memcpy(tM, tm, (m*n)*sizeof(double));
}


void Rn_copy(double* tm,Rn tM,int n)
{
	memcpy(tM, tm, (n)*sizeof(double));
}
 
#endif


#ifndef noeigen	



void Rmn_copy(double* tm, Rmn (&tM),int m,int n)
{
	for(int i=0;i<m;i++)
	for(int j=0;j<n;j++)
		(tM)(i,j)=tm[ i*m+j ];
	
 
}


void Rn_copy(double* tm,Rn  (&ptM),int n)
{
	for(int i=0;i<n;i++)
		(tM)(i)=tm[i];
}

#endif




//============================================//

R testf(Rn x )
{
	return x[0]*x[0]+x[1]*x[1];
}






	void print_vec(Rn x,int n)
	{
		cout<<"[";
		for(int i=0;i<n-1;i++)
		{
			cout<<x[i]<<","<<endl;
		}
		cout<<x[n-1]<<"]"<<endl<<endl;
	}



void Pro  ::  set_daef(dae_f fxut )
{
	this->fxut=fxut;
}

void Pro :: set_opt_int(opt_int jj_opt )
{
	this->j_opt=jj_opt;	
}

void Pro :: set_opt_phi(opt_phi tphi_opt )
{
	this->phi_opt=tphi_opt;	
}


void Pro :: set_opt_pLp(	pfpx tpLpx,	pfpx tpLpu)
{
	this->pLpx=tpLpu;	
	this->pLpu=tpLpu;	
}





void Sol :: set()
{

	this->dimx=this->tPro->dimx;	
	this->dimu=this->tPro->dimu;	

	this->t0=this->tPro->t0;	
	this->tf=this->tPro->tf;	
	
	this->fxut=this->tPro->fxut;	
	this->opt_L=this->tPro->j_opt;	
	this->opt_Phi=this->tPro->phi_opt;	
	this->pLpx=this->tPro->pLpu;	
	this->pLpu=this->tPro->pLpu;	
	
	this->numB=this->tPro->numB;	
	this->numC=this->tPro->numC;	
	
	this->b_c=this->tPro->b_c;	
	this->p_c=this->tPro->p_c;	
	
	
	
}

// void Sol :: set2(fff f2)
// {
	// this->f=f2;

// }


void Sol :: display()
{
	
	cout<<"dimx:"<<(this->dimx)<<endl;
	cout<<"dimu:"<<(this->dimu)<<endl;
	
	cout<<"[t0,tf]:"<<this->t0<<","<<this->tf<<endl;

	
	// double x[1]={1};
	// double  u[1]={1};

	// cout<<this->opt_L(x,u,0)<<endl;
}




R Sol :: testf(Rn x)
{
	double val=0;
	for(int i=0; i<dimx;i++)
	{val=val+x[i]*x[i];}
	
	return val;
	
}

void Sol :: sol(int N)
{	
	this->N=N;
	cout<<"N:"<<this->N<<endl;
 

	this->h=(this->tf-this->t0)/this->N;
 
	this->tk=Rninit(N+1);

	for(int i=0;i<N+1;i++)
	{this->tk[i]=this->t0+i *(this->h) ;}

	int Ntheta=dimx*(N+1)+dimu*(N+1)+2;
	
	
}

R Sol:: objfun( Rn theta, Rn grad)
{

	double val=0;
	Rn tx=Rninit(dimx);
	
	
	Rn tu=Rninit(dimu);


	
	
	

	for(int i=0;i< this->N          ;i++)
	{
		theta2x( theta,tx,i, this->dimx    );
		theta2u(theta,tu,i, this->dimx,this->dimu,N    );
		val=val+(this->tk[i+1]-this->tk[i]   )*opt_L(tx,tu,this->tk[i]);
	}
	
	Rn xt0=Rninit(dimx);
 
	Rn xtf=Rninit(dimx);
 
		
	theta2x( theta,xt0,0, this->dimx    );
	theta2x( theta,xtf,N, this->dimx    );
	
	val=val+opt_Phi(xt0 ,     tk[0],  xtf ,    tk[N]);















	#ifdef noeigen
		free(xtf);
		free(xt0);
		free(tx);
		free(tu);
	#endif

return val;
}



Rn Sol::hk_th( Rn theta,int k)
{
	

	Rn x_k=Rninit(this->dimx);
	Rn u_k=Rninit(this->dimu);
	Rn x_k1=Rninit(this->dimx);
	Rn u_k1=Rninit(this->dimu);
	
	theta2x(theta,x_k, k, this->dimx    );
	theta2u(theta,u_k,  k,this->dimx,this->dimu,this->N    );
	theta2x(theta,x_k1, k+1, this->dimx    );
	theta2u(theta,u_k1,  k+1,this->dimx,this->dimu,this->N    );
	
	print_vec(x_k,dimx);
	print_vec(u_k,dimu);
	print_vec(x_k1,dimx);
	print_vec(u_k1,dimu);
 
	Rn Delta_k=Rninit(this->dimx);
	
	axbyc(1,x_k,1,x_k1,Delta_k,this->dimx);
 
	print_vec(Delta_k,this->dimx);

	Rn xc=Rninit(this->dimx);
	
	xc=x_c(this->fxut,x_k, u_k,0, x_k1, u_k1,this->h);

	print_vec(xc,this->dimx);
	cout<<this->h;
	return x_k;
}


		Rn Sol :: x_c(dae_f temf,Rn xk,Rn uk,double tk,Rn xk1,Rn uk1,double tk1)
		{
			Rn xc=Rninit(this->dimx);

			Rn tem1=Rninit(this->dimx);
			Rn tem2=Rninit(this->dimx);
			axbyc(0.5,xk,-0.5,xk1,tem1,this->dimx   );
			
			
			axbyc((this->h)/8,temf(xk,uk,tk ),-(this->h)/8,temf(xk1,uk1,tk1 ),tem2,this->dimx   );
			
 
			axbyc(1,tem1,1,tem2,xc  ,this->dimx     );
			return xc;
		}
	
	
	void ga_opti::gadisplay()
	{
		for(int i=0;i<this->mpop_size;i++)
		{
			cout<<endl;
			print_vec( this->group[i].state,this->dim );
			cout<<this->group[i].val<<"   val"<<endl ;		
			cout<<this->group[i].fitness<<"   fit"<<endl ;	
			cout<<this->group[i].prob<<"   prob"<<endl ;	
			cout<<this->group[i].sum_prob<<"   sprob"<<endl ;		
			cout<<this->group[i].cho_times<<"   cho"<<endl ;	
			
		}
	}
	
	
	
	
	
void ga_opti ::set_lbub(Rn llb,Rn uub)
{
	this->lb=llb;
	this->ub=uub;
	
}

void ga_opti ::initialize()
{

	
	
	
	individual* newgroup=new  individual[this->mpop_size];
	//individual* newgroup[this->mpop_size]; 不是全局的
	
	
	cout<<endl;cout<<endl;cout<<endl;cout<<endl;
	
	srand(time(NULL) );
	srand(rand() );
	
	for(int i=0;i<this->mpop_size;i++)
	{
		srand(rand() );
		newgroup[i].state=Rninit(this->dim);
		
		
		
		
		for(int j=0;j<this->dim;j++)
		{
			newgroup[i].state[j]=this->lb[j]+(this->ub[j]   -  this->lb[j] )*    (rand()%1000)/1000.0 ;			
		}


	}


	
	this->group=newgroup;
	



}

double ga_opti:: myfit(double tval)
{
	
	return atan(-tval)+3.5;
	
}


bool ga_cmp(individual &a, individual &b )
{
	return a.fitness>=b.fitness;
}


void ga_opti::comp_val_fitness_prob()
{
	this->sum_fit=0;

	for(int i=0;i<this->mpop_size;i++)
	{
		this->group[i].val=testf( this->group[i].state );
		this->group[i].fitness=myfit(this->group[i].val); 
		this->sum_fit=this->sum_fit+	this->group[i].fitness;
		
		
	}	
	
	
	double sum_prob=0;
	for(int i=0;i<this->mpop_size;i++)
	{
		this->group[i].prob=this->group[i].fitness /  this->sum_fit      ;
		sum_prob=sum_prob+this->group[i].prob;
		this->group[i].sum_prob=sum_prob ;
		this->group[i].cho_times=0;
	}	
	
}





bool cmp_mi(double a,double b)
{
	return a<=b;
}
Rn randmi(int NN)
{
		Rn mi=Rninit(NN);
		
		for(int i=0;i<NN;i++)
			mi[i]=rand()%100/100.0;
		
	
		
		#ifndef  noeigen
		sort(mi.data(),mi.data()+ NN,cmp_mi);
		#endif
		
		#ifdef  noeigen
		sort(mi,mi+NN,cmp_mi);
		#endif
		
		return mi;
}



void ga_opti:: chooose()
{
	
	double pr;
	comp_val_fitness_prob();
	Rn mi;
	
	
	
	for(int i=0;i<this->mpop_size;i++)
	{
		
		
		mi=randmi(this->mpop_size);
		print_vec(mi,  this->mpop_size );
		
		
		srand(rand());
		for(int j=0;j<this->mpop_size;j++)
		{
			pr=  rand()%100/100.0;
			cout<<pr<<endl;
			if( pr<=mi[j] )
			{
				++(this->group[j].cho_times);
				break;
				
			}
		}
		
		#ifdef noeigen
		delete(mi);
		#endif
	}









	

gadisplay();
 
	
	// sort(this->group,this->group+this->mpop_size,  ga_cmp );
	
// cout<<endl;
// cout<<endl;
// cout<<endl;
// cout<<endl;

// gadisplay();
 


	// gadisplay();
	
	
	
}

R L_zero(Rn x ,      Rn u     ,double        t )
{
	return 0;
}



void Eu_Lode_Sol::set(Rmn tA,Rmn tB,Rn_f tu ,Rn x0)
{
	this->tA=tA;
	this->tB=tB;
	this->tu=tu;
	this->x0=x0;
 
}



void Eu_Lode_Sol:: sol( ) 
{
Rn temxt=Rninit(dimx);
	
	double tnow;
	cout<<"N:"<<(this->lenN)<<endl;
	this->h=(this->tf-this->t0)/(this->lenN);
	
	
	(this->solt)=(R*)malloc(lenN*sizeof(double ));
	(this->solx)=(Rn*)malloc(lenN*sizeof(Rn));
	
	
	
	
	(this->solx)[0]=x0;
	(this->solt)[0]=this->t0;
	tnow=this->t0;
	
	int itr=1;
	while(itr<lenN)
	{
		
		(this->solx)[itr]=Rninit(dimx);
		
		
		aMx_cC(1,(this->tA),(this->solx)[itr-1] ,0,(this->solx)[itr]   ,dimx);
		aMx_cC(1,(this->tB),(this->tu)(tnow) ,1,(this->solx)[itr]   ,dimx);
		axbyc( 1,(this->solx)[itr-1],1,temxt,   (this->solx)[itr]     ,dimx );
		
		
		
		
		tnow=tnow+h;
		itr++;
	}
	
	cout<<dimx<<dimu<<(this->tA)<<endl;
	
	print_vec((this->solx)[itr-1],dimx);
	
	
}




