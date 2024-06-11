#include "iteration.h"

void iterate_cg(double eps, addeval_func addeval, add_func add,
                scal_func scal, clear_func clear,
                dotprod_func dotprod, norm_func norm, void *data,
                void *r, void *p, void *a, void *b,
                void *x)
{
  clear(p);
  printf("nrm(b):%f\n", norm(b));
  add(1.0, b, p);
  addeval(-1.0, data, x, p);
  int niteration = 0;
  double error = norm(p);
  printf("error&eps: %f &%f\n", error,eps);
  while (error > eps && niteration < 100){
    printf("iteration %i\n", niteration);
    clear(a);
    addeval(1.0,data, p,a);
    double alpha = error*error/dotprod(p,a);
    add(alpha, p,x );
    add(-alpha,a,p);

    error = norm(p);
    niteration++;
  }
  


}


void iterate_gradient(double eps, addeval_func addeval, add_func add,
                scal_func scal, clear_func clear,
                dotprod_func dotprod, norm_func norm, void *data,
                void *r, void *p, void *a, void *b,
                void *x)
{
  clear(p);
  printf("nrm(b):%f\n", norm(b));
  add(1.0, b, p);
  addeval(-1.0, data, x, p);
  int niteration = 0;
  double error = norm(p);
  printf("error&eps: %f &%f\n", error);
  while (error > eps && niteration < 100){
    printf("iteration %i ", niteration);
    clear(a);
    addeval(1.0,data, p,a);
    double alpha = error*error/dotprod(p,a);
    add(alpha, p,x );
    add(-alpha,a,p);

    error = norm(p);
    niteration++;
    printf("error: %f \n", error);
  }
  


}
