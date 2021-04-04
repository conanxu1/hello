/*
https://nlopt.readthedocs.io/en/latest/NLopt_Reference/

*/

#include <stdio.h>
#include <math.h>
#include "nlopt.h"
#define INF (1.0/0.0)

double utility(unsigned n, const double *x, double *grad, void *data){
  grad[0]=1.0/x[0];
  grad[1]=1.0/x[1];
  printf("%f, %f, %f ", x[0],x[1],log(x[0])+log(x[1]));
  return log(x[0])+log(x[1]);
}

double constraint(unsigned n, const double *x, double *grad, void *data){
  double *p=(double *)data;
  grad[0]=*p;
  grad[1]=*(p+1);
  printf("Constraint: %f\n", x[0]*(*p)+x[1]*(*(p+1))-5);
  return x[0]*(*p)+x[1]*(*(p+1))-5;
}

double inconstraint(unsigned n, const double *x, double *grad, void *data){
  grad[0]=1;
  grad[1]=-1;
  return x[0]-x[1];
}

int main(int argc, char const *argv[]) {
  double p[2]={1,2};
  double tol=1e-8;
  double lb[2]={0,0};
  double ub[2]={INF,INF};
  double x[2]={1,1};
  double f_max=-INF;
  // set up optimizer
  //初始化
  nlopt_opt opter=nlopt_create(NLOPT_LD_SLSQP, 2);
  
  // lower and upper bound
  //上下界
  nlopt_set_lower_bounds(opter, lb);
  nlopt_set_upper_bounds(opter, ub);
  
  // objective function
  //目标函数
  nlopt_set_max_objective(opter, utility, NULL);
  
  // equality constraint
  //等式约束
  nlopt_add_equality_constraint(opter, constraint, p, tol);
  
  // inequality constraint
  //不等式约束
  nlopt_add_inequality_constraint(opter, inconstraint, NULL, tol);
  
  // stopping criterion
  //停止准则
  nlopt_set_xtol_rel(opter, tol);
  nlopt_set_ftol_abs(opter, tol);
  nlopt_set_force_stop(opter, tol);
  
  
  // optimize
  nlopt_result result=nlopt_optimize(opter, x, &f_max);
  if (result)
    printf("Maximum utility=%f, x=(%f,%f)\n", f_max, x[0], x[1]);

  // free
  //释放 
 nlopt_destroy(opter);
  return 0;
}
