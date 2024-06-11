#include "gridfunction.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "basic.h"
#include "iteration.h"

/**
 * @brief Set reasonable initial values to a grid function.
 *
 * @param x x-coordinate
 * @param y y-coordinate
 * @return evaluates a continous function f at f(x,y).
 */
static double initial(double x, double y)
{
  double r, z;
  r = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
  if (r < 0.01)
  {
    z = 200.0 * r - 2.0;
  }
  else
  {
    z = 0.0;
  }
  return z;
}

/**
 * @brief Set reasonable initial values to a grid function.
 *
 * @param x x-coordinate
 * @param y y-coordinate
 * @return evaluates a continous function f at f(x,y).
 */
static double linear(double x, double y)
{
  double z = x + y;

  return z;
}

/**
 * @brief Set reasonable boundary values to a grid function.
 *
 * @param x x-coordinate
 * @param y y-coordinate
 * @return evaluates a continous function f at f(x,y).
 */
static double zero(double x, double y)
{
  double z = 0.0;

  return z;
}

int main(void)
{

  pstopwatch sw;
  pgridfunction u, v;
  pgridfunction r, p, a, x;
  int n;
  double eps_solve, t;
  double norm;

  ////////////////////////////////////////////////////////////////////////////////
  // Setup parameters
  ////////////////////////////////////////////////////////////////////////////////

  n = 500;            // Number of inner gridpoints per dimension
  eps_solve = 1.0e-6; // accuracy for solving linear systems with CG-method

  sw = new_stopwatch();

  ////////////////////////////////////////////////////////////////////////////////
  // Setup gridfunctions
  ////////////////////////////////////////////////////////////////////////////////

  printf("Discretizing [0,1] x [0,1] by a grid with %d x %d (%d) points and mesh width %.1e\n", n, n, n * n, 1.0 / (n + 1.0));
  printf("Initializing grid functions\n");
  start_stopwatch(sw);
  u = new_gridfunction(n);
  v = new_gridfunction(n);
  a = new_gridfunction(n);
  r = new_gridfunction(n);
  p = new_gridfunction(n);
  x = new_gridfunction(n);
  clear_gridfunction(v);
  clear_gridfunction(a);
  clear_gridfunction(r);
  clear_gridfunction(p);
  clear_gridfunction(x);
  t = stop_stopwatch(sw);
  printf("  %.2f ms\n", t * 1.0e3);

  ////////////////////////////////////////////////////////////////////////////////
  // testing addevalLaplace
  ////////////////////////////////////////////////////////////////////////////////

  printf("Testing discrete Lapalace operator:\n");
  start_stopwatch(sw);
  init_gridfunction(&linear, u);
  boundary_gridfunction(&linear, u);
  addevalLaplace_gridfunction(1.0, NULL, u, v);
  t = stop_stopwatch(sw);
  printf("  %.2f ms\n", t * 1.0e3);
  printf("  ||u|| = %.5e\n", norm_gridfunction(u) / (n * n));
  norm = norm_gridfunction(v) / (n * n);
  printf("  ||v|| = %s%.5e" NORMAL "\n", (norm < 1.0e-12) ? BGREEN : ((norm < 1.0e-8) ? BYELLOW : BRED), norm); // Should be near zero if laplacian works correctly.

  ////////////////////////////////////////////////////////////////////////////////
  // Solving linear system with cg-method
  ////////////////////////////////////////////////////////////////////////////////

  printf("Solve for initial value using cg-method:\n");
  start_stopwatch(sw);
  init_gridfunction(&initial, u);
  boundary_gridfunction(&zero, u);
  clear_gridfunction(v);
  addevalLaplace_gridfunction(1.0, NULL, u, v);
  
  iterate_cg(eps_solve, (addeval_func)addevalLaplace_gridfunction, (add_func)add_gridfunction,
             (scal_func)scal_gridfunction, (clear_func)clear_gridfunction, (dotprod_func)dotprod_gridfunction,
             (norm_func)norm_gridfunction, NULL, (void *)r, (void *)p, (void *)a, (void *)v, (void *)x);
  t = stop_stopwatch(sw);
  printf("  %.2f ms\n", t * 1.0e3);
  add_gridfunction(-1.0, u, x);
  norm = norm_gridfunction(x) / (n * n) / norm_gridfunction(u);
  printf("  ||u - \\tilde u|| / ||u|| = %s%.5e"NORMAL"\n", (norm < 1.0e-15) ? BGREEN : ((norm < 1.0e-11) ? BYELLOW : BRED), norm);

  ////////////////////////////////////////////////////////////////////////////////
  // Cleanup
  ////////////////////////////////////////////////////////////////////////////////

  printf("Cleaning up\n");
  start_stopwatch(sw);
  del_gridfunction(u);
  del_gridfunction(v);
  del_gridfunction(r);
  del_gridfunction(a);
  del_gridfunction(p);
  del_gridfunction(x);

  t = stop_stopwatch(sw);
  printf("  %.2f ms\n", t * 1.0e3);
  del_stopwatch(sw);

  return EXIT_SUCCESS;
}
