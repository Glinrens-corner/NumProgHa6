#ifndef ITERATION_H
#define ITERATION_H

#include <stdio.h>

typedef void (*addeval_func)(double a, void *data, void *x, void *y); // y <-- y + a A x, data stores information for evaluation of A
typedef void (*add_func)(double a, void *x, void *y);                 // y <-- y + ax
typedef void (*scal_func)(double a, void *x);                         // x <-- ax
typedef void (*clear_func)(void *x);                                  // x <-- 0
typedef double (*dotprod_func)(void *x, void *y);                     // Compute <x,y>
typedef double (*norm_func)(void *x);                                 // Compute ||x||

/* General conjugate gradient method using callback functions.
 * The accuracy eps is used for a termination criterion.
 * Confer lecture notes for further explanations. 
 * 'data' is being passed to the callback 'addeval'.*/
void iterate_cg(double eps, addeval_func addeval, add_func add,
                scal_func scal, clear_func clear,
                dotprod_func dotprod, norm_func norm, void *data,
                void *r, void *p, void *a, void *b, void *x);


#endif
