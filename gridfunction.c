#include "gridfunction.h"
#include "miniblas.h"
#include "iteration.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/***************************
Constructors and destructors
***************************/

pgridfunction new_gridfunction(int n)
{
  pgridfunction u;

  u = (pgridfunction)calloc(1, sizeof(gridfunction));
  u->n = n;
  u->h = 1.0 / (n + 1.0);
  u->val = (double *)calloc(1, sizeof(double) * (n + 2) * (n + 2));

  return u;
}

void del_gridfunction(pgridfunction u)
{
  if (u)
  {
    free(u->val);
    free(u);
  }
}

void copy_gridfunction(pcgridfunction u, pgridfunction v)
{
  int i, j, n, incy;
  double *uval, *vval;

  n = u->n;
  uval = u->val;
  vval = v->val;

  assert(n == v->n);

  for (j = 1; j <= n; j++)
  {
    incy = j * (n + 2);
    for (i = 1; i <= n; i++)
    {
      vval[i + incy] = uval[i + incy];
    }
  }
}

void init_gridfunction(function f, pgridfunction u)
{
  int n, i, j, incy;
  double h, *val;

  n = u->n;
  val = u->val;
  h = u->h;

  for (j = 1; j <= n; j++)
  {
    incy = j * (n + 2);
    for (i = 1; i <= n; i++)
    {
      val[i + incy] = f(i * h, j * h);
    }
  }
}

void boundary_gridfunction(function f, pgridfunction u)
{
  int n, i, j, incy;
  double h, *val;

  n = u->n;
  val = u->val;
  h = u->h;

  for (j = 1; j <= n; j++)
  {
    incy = j * (n + 2);
    i = 0;
    val[i + incy] = f(i * h, j * h);
    i = n + 1;
    val[i + incy] = f(i * h, j * h);
  }
  for (i = 0; i <= n + 1; i++)
  {
    j = 0;
    incy = j * (n + 2);
    val[i + incy] = f(i * h, j * h);
    j = n + 1;
    incy = j * (n + 2);
    val[i + incy] = f(i * h, j * h);
  }
}

/**********************
In- and output
**********************/

void saveToFile_gridfunction(pcgridfunction u, char *filename)
{
  FILE *fp;
  int i, j, n;
  double h, *uval;

  n = u->n;
  uval = u->val;
  h = u->h;

  fp = fopen(filename, "w");
  for (i = 0; i <= n + 1; i++)
  {
    for (j = 0; j <= n + 1; j++)
    {
      fprintf(fp, "%lf  %lf  %lf\n", i * h, j * h, uval[i + j * (n + 2)]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

void output_gridfunction(pcgridfunction u, char *iteration, int frame)
{
  char filename[255];

  sprintf(filename, "%s/wave2d_%s_%05d.dat", iteration, iteration, frame);

  saveToFile_gridfunction(u, filename);
}

/****************************
Basic linear algebra routines
****************************/

void add_gridfunction(double a, pcgridfunction u, pgridfunction v)
{
  int n, j, incy;
  double *uval, *vval;

  n = u->n;
  uval = u->val;
  vval = v->val;

  assert(n == v->n);

  for (j = 1; j <= n; j++)
  {
    incy = j * (n + 2);
    axpy(n, a, uval + incy + 1, 1, vval + incy + 1, 1);
  }
}

double dotprod_gridfunction(pcgridfunction u, pcgridfunction v)
{
  int n, j, incy;
  double prod, *uval, *vval;

  n = u->n;
  uval = u->val;
  vval = v->val;
  assert(n == v->n);

  prod = 0.0;
  for (j = 1; j <= n; j++)
  {
    incy = j * (n + 2);
    prod += dot(n, uval + incy + 1, 1, vval + incy + 1, 1);
  }

  return prod;
}

double norm_gridfunction(pcgridfunction u)
{
  return sqrt(dotprod_gridfunction(u, u));
}

void scal_gridfunction(double a, pgridfunction u)
{
  int n, j, incy;
  double *uval;

  n = u->n;
  uval = u->val;

  for (j = 1; j <= n; j++)
  {
    incy = j * (n + 2);
    scal(n, a, uval + incy + 1, 1);
  }
}

void clear_gridfunction(pgridfunction u)
{
  scal_gridfunction(0.0, u);
}

void addevalLaplace_gridfunction(double a, void *data, pcgridfunction u,
                                 pgridfunction v)
{
  int n;
  double *uval, *vval, h /*, b*/;

  n = u->n;
  h = u->h;
  vval = v->val;
  uval = u->val;

  (void)data;

  assert(n == v->n && h == v->h);

  double factor = a/(h*h);


  int ld = n+2;
  
  for (int icol = 1; icol < n+1 ; icol++){
    int copy_start_offset =1+icol*(n+2);
    axpy(n, 4.0*factor, uval+ copy_start_offset, 1, vval +copy_start_offset,1);
    axpy(n, -1.0*factor, uval+copy_start_offset+ld,1, vval+copy_start_offset,1 );
    axpy(n, -1.0*factor, uval+copy_start_offset-ld,1, vval+copy_start_offset,1 );
    axpy(n, -1.0*factor, uval+copy_start_offset-1,1, vval+copy_start_offset,1 );
    axpy(n, -1.0*factor, uval+copy_start_offset+1,1, vval+copy_start_offset,1 );
    

  }
}

/*****************************
 * Callback functions
 *****************************/

void(add_gridfunc)(double a, void *x, void *y)
{
  add_gridfunction(a, (pcgridfunction)x, (pgridfunction)y);
}

void(scal_gridfunc)(double a, void *x)
{
  scal_gridfunction(a, (pgridfunction)x);
}

void(clear_gridfunc)(void *x)
{
  clear_gridfunction((pgridfunction)x);
}

double(dotprod_gridfunc)(void *x, void *y)
{
  return dotprod_gridfunction((pcgridfunction)x, (pcgridfunction)y);
}

double(norm_gridfunc)(void *x)
{
  return norm_gridfunction((pcgridfunction)x);
}
