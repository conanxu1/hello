/*
 * main.c
 *
 *  Created on: Oct 9, 2018
 *      Author: lgh
 */


#include <stdio.h>
#include <iostream>

#include <math.h>
#include "nlopt/include/nlopt.h"
using namespace std;

#define INF (1.0/0.0)
int i=0;

//目标函数；
double utility(unsigned n, const double *x, double *grad, void *data)
{
    if(grad){
      grad[0]=2*x[0];
      grad[1]=2*x[1];
    }
    printf("迭代次数 i= %d, x[0]=%f, x[1]= %f，f(x1,x2)=%f\n",i++,x[0],x[1],x[0]*x[0]+x[1]*x[1]+8);
  return ( x[0]*x[0]+x[1]*x[1]+8 );
}

//等式限制条件；
double constraint(unsigned n, const double *x, double *grad, void *data)
{
    if(grad){
          grad[0]= -1.0;
          grad[1]= -2*x[1];
    }

  return (-x[0]-x[1]*x[1]+2);
}

//不等式限制条件；
double inconstraint(unsigned n, const double *x, double *grad, void *data)
{
    if(grad){
     grad[0]= -2*x[0];
     grad[1]= 1.0;
    }
  return (-x[0]*x[0]+x[1]);
}


int main(int argc, char const *argv[])
{
  double tol=1e-8;
  double lb[2]={0,0};          //x1、x2的下边界；
  double ub[2]={INF,INF};
  double x[2]={1, 1};      //给x1、x2赋予初始值；
  double f_max;

  nlopt_opt opter=nlopt_create( NLOPT_LD_SLSQP, 2);

  //设置自变量下限；
  nlopt_set_lower_bounds(opter, lb);

  // 目标函数；
  nlopt_set_min_objective(opter, utility, NULL);

  // 不等式约束；
  nlopt_add_inequality_constraint(opter, inconstraint, NULL, tol);

  // 等式约束；
  nlopt_add_equality_constraint(opter, constraint, NULL, tol);

  // 停止时需要的条件；
  nlopt_set_xtol_rel(opter, tol);

  // 开始优化；
  nlopt_result result=nlopt_optimize(opter, x, &f_max);
cout<<"k";
  if (result)
  {
    printf("目标函数最大值=%g, x=(%g,%g)\n", f_max, x[0], x[1]);
  }

  //free
  nlopt_destroy(opter);
  return 0;
}