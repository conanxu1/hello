/*
https://nlopt.readthedocs.io/en/latest/NLopt_Reference/

*/

#include <stdio.h>
#include <math.h>
#include "nlopt.h"
#define INF (1.0/0.0)

double utility(unsigned n, const double *x, double *grad, void *data){
  grad[0]= 0.5*x[0];
  grad[1]= 0.5*x[1];
  
  
  
  printf("%f, %f\n", x[0],x[1]);
  return pow(x[0],2) + pow(x[1],2);
}






int main(int argc, char const *argv[]) {

  double tol=1e-8;
  double lb[2]={-10,-10};
  double ub[2]={10,10};
  double x[2]={1,3};
  double f_max=0;
  // set up optimizer
  //初始化
  nlopt_opt opter=nlopt_create(NLOPT_LD_SLSQP, 2);
  
  // lower and upper bound
  //上下界
  nlopt_set_lower_bounds(opter, lb);
  nlopt_set_upper_bounds(opter, ub);
  
  // objective function
  //目标函数
  nlopt_set_min_objective(opter, utility, NULL);
  



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
