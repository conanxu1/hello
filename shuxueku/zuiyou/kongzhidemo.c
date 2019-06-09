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


lyRnR *kkzcek,*kkzcik,*kkzphi,*kkzpsi;
int kkznumcek=0,kkznumcik=2,kkznumphii=2, kkznumpsii=0;
piandao kkzgPHI,kkzgg,kkzgf,*kkzgcek, *kkzgphi,*kkzgcik,*kkzgpsi,

int geshu;

int dimx=3,dimu=1;
int Ntau=3
double *tauk,*wk,*Dki,*zuiyouX;



kkzcek=malloc(kznumcek*sizeof(lyRnR));
kkzcik=malloc(kznumcik*sizeof(lyRnR));
kkzphi=malloc(kznumphii*sizeof(lyRnR));
kkzpsi=malloc(kznumpsii*sizeof(lyRnR));
kkzgcek=malloc(kznumcek*sizeof(piandao));
kkzgphi=malloc(kznumphii*sizeof(piandao));
kkzgcik=malloc(kznumcik*sizeof(piandao));
kkzgpsi=malloc(kznumpsii*sizeof(piandao));


kkzcik[0]=kzcik1;
kkzcik[1]=kzcik2;

kkzphi[0]=kzphi1;
kkzphi[1]=kzphi2;












double *h,*be,*Ae*bi,*Ai,*XUk,tk;
 



tauk=(double *)malloc((Ntau)*sizeof(double));
tk=(double *)malloc((Ntau)*sizeof(double));



wk=(double *)malloc((Ntau)*sizeof(double));
Dki=(double *)malloc((Ntau*Ntau)*sizeof(double));
lgd(tauk,Ntau-1);
lg_AD(tauk, Ak,Dki,Ntau-1);




zuiyouX=(double *)malloc(((dimx+dimu)*Ntau+1)*sizeof(double));



geshu=(dimx+dimu)*Ntau+1;
h=(double *)malloc(geshu*sizeof(double));

geshu=kznumcek*Ntau+kznumphii+dimx*Ntau;
be=(double *)malloc(geshu*sizeof(double));


geshu=(kznumcek*Ntau+kznumphii+dimx*Ntau)*((dimx+dimu)*Ntau+1);
Ae=(double *)malloc(geshu*sizeof(double));

geshu=kznumcik*Ntau+kznumpsii;
bi=(double *)malloc(geshu*sizeof(double));

geshu=(kznumcik*Ntau+kznumpsii)*((dimx+dimu)*Ntau+1);
Ai=(double *)malloc(geshu*sizeof(double));



geshu=(dimx+dimu)*Ntau+1;
XUk=(double *)malloc(((dimx+dimu)*Ntau+1)*sizeof(double));
 


xishu(kzPHI, kzg,kzf,kkzcek,kkznumcek,kkzphii,kkznumphii,kkzcik,
kkznumcik,kkzpsii,kkznumpsii,kkzgPHI,kkzgg,kkzgf,kkzgcek,kkzgphi,kkzgcik,kkzgpsi,
dimx,dimu,Ntau,tauk,wk, Dki,zuiyouX ,h,be,Ae, bi, Ai,XUk, tk)















return 0;
}
 
 
 
 

