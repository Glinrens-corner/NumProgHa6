#include "iteration.h"

void iterate_cg(double eps, addeval_func addeval, add_func add,
                scal_func scal, clear_func clear,
                dotprod_func dotprod, norm_func norm, void *data,
                void *r, void *p, void *a, void *b,
                void *x)
{
  clear(p);
  clear(r);
  add(1.0, b, r);
  addeval(-1.0, data, x, r);
  add(1.0, r,p);
  
  int niteration = 0;
  double error = norm(p);
  while (error > eps && niteration < 2000){
    clear(a);
    addeval(1.0,data, p,a);
    double beta = dotprod(p,a);
    double alpha = dotprod(p,r)/beta;
    add(alpha, p,x );
    add(-alpha,a,r);
    double gamma = dotprod(a,r)/beta;
    
    scal(-gamma,p);
    add(1.0, r, p);
    
    error = norm(p);

    /* if(niteration %10 == 0){ */
    /*   printf("iteration %i ", niteration); */
    /*   printf("error: %f \n", error); */
    /* } */
    niteration++;
  }
  


}


