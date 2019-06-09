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

printf("\n\n...............\n");
 

lyRnR *kkzcek,*kkzcik,*kkzphi,*kkzpsi;
int kkznumcek=0,kkznumcik=2,kkznumphii=2, kkznumpsii=0;
piandao kkzgPHI,kkzgg,kkzgf,*kkzgcek, *kkzgphi,*kkzgcik,*kkzgpsi;

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

kkzphi[0]=kzphi1;
kkzphi[1]=kzphi2;












double *h,*be,*Ae*bi,*Ai,*XUk,tk,*x0,t0;






x0=(double *)malloc((dimx)*sizeof(double));
x0[0]=0;
x0[1]=1;
x0[2]=0;



t0=0;






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
dimx,dimu,Ntau,tauk,wk, Dki,zuiyouX ,h,be,Ae, bi, Ai,XUk, tk,x0,t0)




printf("\n\n...............\n");
geshu=(dimx+dimu)*Ntau+1;
shuchud(h,geshu,1)
 
 
printf("\n\n...............\n");
 
geshu=kznumcek*Ntau+kznumphii+dimx*Ntau;
shuchud(be,geshu,1);





printf("\n\n...............\n");
 
shuchud(Ae,(kznumcek*Ntau+kznumphii+dimx*Ntau),((dimx+dimu)*Ntau);






printf("\n\n...............\n");
geshu=kznumcik*Ntau+kznumpsii;
 

shuchud(be,geshu,1);



printf("\n\n...............\n");
shuchud(Ae,(kznumcik*Ntau+kznumpsii),((dimx+dimu)*Ntau+1));


 
 
 
printf("\n\n...............\n");
geshu=(dimx+dimu)*Ntau+1;
XUk=(double *)malloc(((dimx+dimu)*Ntau+1)*sizeof(double)); shuchud(XUk,geshu,1);























return 0;
}
 
 
 
 

