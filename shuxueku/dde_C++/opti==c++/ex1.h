#ifndef _PRO_H_
#define _PRO_H_
#include "pro.h"
#endif

#ifndef _EX1_H_
#define _EX1_H_
#include "ex1.h"
#endif


extern int N;
extern int dimx;
extern int dimu;
extern double t0;
extern double tf;






Rn my_fxut(Rn x ,      Rn u     ,double        t );
R my_opt_phi(Rn xt0 ,     double t0, Rn xtf ,     double tf );
R my_opt_L(Rn x ,      Rn u     ,double        t );

Rn pLpx(Rn x ,      Rn u     ,double        t,int dimx,int dimu );
Rn pLpu(Rn x ,      Rn u     ,double        t ,int dimx,int dimu);
Rn myut(double t);




