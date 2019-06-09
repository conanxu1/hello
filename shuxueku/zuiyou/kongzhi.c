#include <math.h>
#include <stdlib.h>
#include <string.h>



//J=phi(x0,t0,xf,tf)+int0 t    g(x,u,t)

//dx=f(x,u,t)

//phi1(x0,t0,xf,tf)=0

//C1(x,u,t)=0

//C2(x,u,t)<=0





double kzPHI(double **xftf)
{


return xftf[0][2];
 }

double kzg(double **xut)
{
	
return 0;
}


double* kzf(double **xut)
{
double *fxut;
fxut=(* double)malloc(3*sizeof(double));	
	
fxut[0]=xut[0][1];	
fxut[1]=xut[1][0];		
fxut[2]=xut[1][0]*xut[1][0]/2;	
	
	
	
	
return 	fxut;
}





double kzphii1(double **xftf)
{
	
	
return  xftf[0][0];

}

double kzphii2(double **xftf)
{
	
	
return  xftf[0][1];

}



double kzcik1(double **xut)
{
	
	
return xut[0][0]-1/9;	
}

double kzcik2(double **xut)
{
	
	
return -xut[0][0];	
}



// lyRnR *psii,
	
	
double* kzgPHI(double **xftf,int flag)
{
	double *pd;
	if(flag==1)
	{
		pd=(double *)malloc(3*sizeof(double));
		memset(pd,0,3*sizeof(double));
		pd[2]=1;
		return pd;
	}
	else if(flag==2)
	{
		pd=(double *)malloc(sizeof(double));
		memset(pd,0,sizeof(double));
		


		return pd;
		
	}
	
	
}


double* kzgg(double **xut,int flag)	
{
	double *pd;
	if(flag==1)
	{
		pd=(double *)malloc(4*sizeof(double));
		memset(pd,0,4*sizeof(double));
		
	
	
		return pd;
	}
	else if(flag==2)
	{
		pd=(double *)malloc(sizeof(double));
		memset(pd,0,sizeof(double));
		


		return pd;
		
	}	
	
	
	
	
}		
	

double* kzgf(double **xut,int flag)	
{
	double *pd;
	if(flag==1)
	{
		pd=(double *)malloc(3*4*sizeof(double));
		memset(pd,0,3*4*sizeof(double));
		
	pd[0*4+1]=1;
	
	pd[1*4+3]=1;
	
	pd[2*4+3]=xut[1][0];
	
	
	
	
	
	
	
	
	
		return pd;
	}
	else if(flag==2)
	{
		pd=(double *)malloc(3*sizeof(double));
		memset(pd,0,3*sizeof(double));
		


		return pd;
		
	}	
	
	
	
	
}		




// double* kzgcek1(double **xut,int flag)	
// {
	// double *pd;
	// int dim;
	
	// if(flag==1)
	// {dim=4;
		// pd=(double *)malloc(dim*sizeof(double));
		// memset(pd,0,dim*sizeof(double));
		
	// pd[0]=1;
	
	
	// }
	// else if(flag==2)
	// {dim=1;
		// pd=(double *)malloc(dim*sizeof(double));
		// memset(pd,0,dim*sizeof(double));
		


		
	// }	
	
	// return pd;
	

// }	


// double* kzgcek2(double **xut,int flag)	
// {
	// double *pd;
	// int dim;
	
	// if(flag==1)
	// {dim=4;
		// pd=(double *)malloc(dim*sizeof(double));
		// memset(pd,0,dim*sizeof(double));
		
	// pd[0]=-1;
	
	
	// }
	// else if(flag==2)
	// {dim=1;
		// pd=(double *)malloc(dim*sizeof(double));
		// memset(pd,0,dim*sizeof(double));
		


		
	// }	
	
	// return pd;
	

// }	




double* kzgphi1(double **xftf,int flag)	
{
	double *pd;
	int dim;
	
	if(flag==1)
	{dim=3;
		pd=(double *)malloc(dim*sizeof(double));
		memset(pd,0,dim*sizeof(double));
		
	pd[0]=1;
	
	
	}
	else if(flag==2)
	{dim=1;
		pd=(double *)malloc(dim*sizeof(double));
		memset(pd,0,dim*sizeof(double));
		


		
	}	
	
	return pd;
	

}	
	
double* kzgphi2(double **xftf,int flag)	
{
	double *pd;
	int dim;
	
	if(flag==1)
	{dim=3;
		pd=(double *)malloc(dim*sizeof(double));
		memset(pd,0,dim*sizeof(double));
		
	pd[1]=1;
	
	
	}
	else if(flag==2)
	{dim=1;
		pd=(double *)malloc(dim*sizeof(double));
		memset(pd,0,dim*sizeof(double));
		


		
	}	
	
	return pd;
	

}	

double* kzgcik1(double **xut,int flag)	
{
	double *pd;
	int dim;
	
	if(flag==1)
	{dim=4;
		pd=(double *)malloc(dim*sizeof(double));
		memset(pd,0,dim*sizeof(double));
		
	pd[0]=1;
	
	
	}
	else if(flag==2)
	{dim=1;
		pd=(double *)malloc(dim*sizeof(double));
		memset(pd,0,dim*sizeof(double));
		


		
	}	
	
	return pd;
	

}	





double* kzgcik2(double **xut,int flag)	
{
	double *pd;
	int dim;
	
	if(flag==1)
	{dim=4;
		pd=(double *)malloc(dim*sizeof(double));
		memset(pd,0,dim*sizeof(double));
		
	pd[0]=-1;
	
	
	}
	else if(flag==2)
	{dim=1;
		pd=(double *)malloc(dim*sizeof(double));
		memset(pd,0,dim*sizeof(double));
		


		
	}	
	
	return pd;
	

}	

	









// double* kzgphi2(double **xftf,int flag)	
// {
	// double *pd;
	// int dim;
	
	// if(flag==1)
	// {dim=3;
		// pd=(double *)malloc(dim*sizeof(double));
		// memset(pd,0,dim*sizeof(double));
		
	// pd[1]=1;
	
	
	// }
	// else if(flag==2)
	// {dim=1;
		// pd=(double *)malloc(dim*sizeof(double));
		// memset(pd,0,dim*sizeof(double));
		


		
	// }	
	
	// return pd;
	

// }	











