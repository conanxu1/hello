#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stddef.h>
#include <malloc.h>
#include <string.h>
#include "myfun.h"

#include <complex.h>



 

int main()
{

int i,j,k,geshu;
lyRnR *kkzcek,*kkzcik,*kkzphi,*kkzpsi;
piandao kkzgPHI,kkzgg,kkzgf,*kkzgcek,*kkzgphi,*kkzgcik,*kkzgpsi;
double *h,*be,*Ae,*bi,*Ai,*XUk,*tk,*x0,t0;
double *tauk,*wk,*Dki,*zuiyouX;
int zongdim;


//----------------------------------------------------------------
//各种约束的个数

int kkznumcek=0,kkznumcik=2,kkznumphii=2, kkznumpsii=0;

kkzcek=malloc(kkznumcek*sizeof(lyRnR));
kkzcik=malloc(kkznumcik*sizeof(lyRnR));
kkzphi=malloc(kkznumphii*sizeof(lyRnR));
kkzpsi=malloc(kkznumpsii*sizeof(lyRnR));
kkzgcek=malloc(kkznumcek*sizeof(piandao));
kkzgphi=malloc(kkznumphii*sizeof(piandao));
kkzgcik=malloc(kkznumcik*sizeof(piandao));
kkzgpsi=malloc(kkznumpsii*sizeof(piandao));



//----------------------------------------------------------------
//问题的维数 几阶勒让德


int dimx=3,dimu=1;
int Ntau=3;
x0=cshi(dimx);

tauk=cshi(Ntau);
tk=cshi(Ntau);
wk=cshi(Ntau);
Dki=cshi(Ntau*Ntau);


h=cshi(Ntau*(dimx+dimu)+1);

double *H;
H=cshi((Ntau*(dimx+dimu)+1)*(Ntau*(dimx+dimu)+1));
for(i=0;i<(Ntau*(dimx+dimu)+1);i++)
	H[i*(Ntau*(dimx+dimu)+1)+i]=1;

//----------------------------------------------------------------
//初始化一些函数

kkzgPHI=kzgPHI;

kkzcik[0]=kzcik1;
kkzcik[1]=kzcik2;

kkzphi[0]=kzphii1;
kkzphi[1]=kzphii2;


kkzgcik[0]=kzgcik1;
kkzgcik[1]=kzgcik2;


kkzgphi[0]=kzgphi1;
kkzgphi[1]=kzgphi2;






t0=0;
x0[0]=0;
x0[1]=1;
x0[2]=0;


 

lgd(tauk,Ntau-1);
lg_AD(tauk, wk,Dki,Ntau-1);

	
shuchud(Dki,Ntau,Ntau);
printf("ee	\n");


zuiyouX=cshi((dimx+dimu)*Ntau+1);



geshu=(dimx+dimu)*Ntau+1;
h=cshi(geshu);


 


geshu=kkznumcek*Ntau+kkznumphii+dimx*Ntau;
be=cshi(geshu);


geshu=(kkznumcek*Ntau+kkznumphii+dimx*Ntau)*((dimx+dimu)*Ntau+1);
Ae=cshi(geshu);

geshu=kkznumcik*Ntau+kkznumpsii;
bi=cshi(geshu);

geshu=(kkznumcik*Ntau+kkznumpsii)*((dimx+dimu)*Ntau+1);
Ai=cshi(geshu);
printf("%d",geshu);


geshu=(dimx+dimu)*Ntau+1;
XUk=cshi(geshu);


//-------------------------------

//首末时间不能一样


for(i=1;i<=Ntau*(dimx+dimu);i++)
XUk[i-1]=i;
 

XUk[(dimx+dimu)*Ntau]=1;
//t0  neq tf
 printf("XUk\n");
 shuchud(XUk,Ntau*(dimx+dimu)+1,1);
 
 


xishu(kzPHI, kzg,kzf,kkzcek,kkznumcek,kkzphi,kkznumphii,kkzcik,
kkznumcik,kkzpsi,kkznumpsii,kkzgPHI,kzgg,kzgf,kkzgcek,kkzgphi,kkzgcik,kkzgpsi,
dimx,dimu,Ntau,tauk,wk, Dki,zuiyouX ,h,be,Ae, bi, Ai,XUk, tk,x0,t0);

zongdim=Ntau*(dimx+dimu)+1;

qxt(H,	h,be ,Ae,	bi,	Ai,zongdim,	kkznumcek*Ntau+dimx*Ntau+kkznumphii,	kkznumcik*Ntau+kkznumphii	);
















/*



 
printf("\n\n...............\n");
 
geshu=kkznumcek*Ntau+kkznumphii+dimx*Ntau;
shuchud(be,geshu,1);





printf("\n\n...............\n");
 
shuchud(Ae,(kkznumcek*Ntau+kkznumphii+dimx*Ntau),(dimx+dimu)*Ntau);






printf("\n\n...............\n");
geshu=kkznumcik*Ntau+kkznumpsii;
 

shuchud(be,geshu,1);



printf("\n\n...............\n");
shuchud(Ae,(kkznumcik*Ntau+kkznumpsii),((dimx+dimu)*Ntau+1));


 
 
 
printf("\n\n...............\n");
geshu=(dimx+dimu)*Ntau+1;
XUk=(double *)malloc(((dimx+dimu)*Ntau+1)*sizeof(double)); shuchud(XUk,geshu,1);



*/



















return 0;
}
 
 
 
 

