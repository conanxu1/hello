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

 
int i,j,k;

lyRnR *kkzcek,*kkzcik,*kkzphi,*kkzpsi;
int kkznumcek=0,kkznumcik=2,kkznumphii=2, kkznumpsii=0;
piandao kkzgPHI,kkzgg,kkzgf,*kkzgcek,*kkzgphi,*kkzgcik,*kkzgpsi;


kkzgPHI=kzgPHI;





int geshu;

int dimx=3,dimu=1;
int Ntau=3;
double *tauk,*wk,*Dki,*zuiyouX;



kkzcek=malloc(kkznumcek*sizeof(lyRnR));
kkzcik=malloc(kkznumcik*sizeof(lyRnR));
kkzphi=malloc(kkznumphii*sizeof(lyRnR));
kkzpsi=malloc(kkznumpsii*sizeof(lyRnR));
kkzgcek=malloc(kkznumcek*sizeof(piandao));
kkzgphi=malloc(kkznumphii*sizeof(piandao));
kkzgcik=malloc(kkznumcik*sizeof(piandao));
kkzgpsi=malloc(kkznumpsii*sizeof(piandao));


kkzcik[0]=kzcik1;
kkzcik[1]=kzcik2;

kkzphi[0]=kzphii1;
kkzphi[1]=kzphii2;


kkzgcik[0]=kzgcik1;
kkzgcik[1]=kzgcik2;


kkzgphi[0]=kzgphi1;
kkzgphi[1]=kzgphi2;







double *h,*be,*Ae,*bi,*Ai,*XUk,*tk,*x0,t0;




h=cshi(Ntau*(dimx+dimu)+1);

x0=(double *)malloc((dimx)*sizeof(double));
x0[0]=0;
x0[1]=1;
x0[2]=0;




double *ioi,*ioi1;
ioi=cshi(1);
ioi1=cshi(1);
 

free(ioi);
ioi=NULL;
ioi=cshi(1);

free(ioi1);
ioi1=NULL;
ioi1=cshi(dimx);

 
 cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,1,1,dimx, 1,x0, dimx,ioi1,1, 1,ioi, 1);


 
	








t0=0;






tauk=(double *)malloc((Ntau)*sizeof(double));
tk=(double *)malloc((Ntau)*sizeof(double));



wk=(double *)malloc((Ntau)*sizeof(double));
Dki=(double *)malloc((Ntau*Ntau)*sizeof(double));
lgd(tauk,Ntau-1);


lg_AD(tauk, wk,Dki,Ntau-1);

	


zuiyouX=(double *)malloc(((dimx+dimu)*Ntau+1)*sizeof(double));



geshu=(dimx+dimu)*Ntau+1;
h=cshi(geshu);


shuchud(h,geshu,1);




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


for(i=1;i<=Ntau;i++)
XUk[(dimx+dimu)*i-1]=1;
 

XUk[(dimx+dimu)*Ntau]=1;
//t0  neq tf
 


xishu(kzPHI, kzg,kzf,kkzcek,kkznumcek,kkzphi,kkznumphii,kkzcik,
kkznumcik,kkzpsi,kkznumpsii,kkzgPHI,kzgg,kzgf,kkzgcek,kkzgphi,kkzgcik,kkzgpsi,
dimx,dimu,Ntau,tauk,wk, Dki,zuiyouX ,h,be,Ae, bi, Ai,XUk, tk,x0,t0);






printf("\n\ntauk...............\n");
shuchud(tauk,Ntau,1);

printf("\n\ntauk...............\n");
shuchud(tauk,Ntau,1);










printf("\n\nbe...............\n");
geshu=(kkznumcek*Ntau+kkznumphii+dimx*Ntau) ;
shuchud(be,geshu,1);
  

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
 
 
 
 

