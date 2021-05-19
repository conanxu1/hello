#include "pro.h"


Sol* pmysol;
int N=2;
int dimx=2;
int dimu=1;
double t0=0;
double tf=2;
Rn temut=Rninit(dimu);

Rn my_fxut(Rn x ,      Rn u     ,double        t );
R my_opt_phi(Rn xt0 ,     double t0, Rn xtf ,     double tf );
R my_opt_L(Rn x ,      Rn u     ,double        t );

Rn pLpx(Rn x ,      Rn u     ,double        t,int dimx,int dimu );
Rn pLpu(Rn x ,      Rn u     ,double        t ,int dimx,int dimu);
Rn myut(double t);




