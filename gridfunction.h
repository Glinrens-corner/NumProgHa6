#ifndef GRIDFUNCTION_H
#define GRIDFUNCTION_H

/*************************
Definitions
*************************/

/* Define gridfunction object. */
typedef struct _gridfunction gridfunction;
/* Define pointer to gridfunction object. */
typedef gridfunction *pgridfunction;
/* Define pointer to constant gridfunction object. */
typedef const gridfunction *pcgridfunction;

struct _gridfunction
{

  /* Number of inner grid points in each direction. In total there are
     (n+2)*(n+2) many points including the boundary, but only n*n inner points.*/
  int n;

  /* Mesh width. */
  double h;

  /* Function values for all grid points, arranged in column-major order.
   * For example, (u(0,0),u(1,0),u(2,0),u(0,1),u(1,1),u(2,1),u(0,2),u(1,2),u(2,2)). */
  double *val;
};

/* Define function object for real functions on a two-dimensional space. */
typedef double (*function)(double, double);

/***************************
Constructors and destructors
***************************/

/* Allocate memory for a new gridfunction object and set up parameters.
 * n: Number of inner grid points in each direction.
 * returns: Address of new gridfunction object.
 * notes: Values will by default be set to zero. */
pgridfunction new_gridfunction(int n);

/* Delete a gridfunction object.
 * u: Gridfunction object to be deleted. */
void del_gridfunction(pgridfunction u);

/* Copy a grid function.
 * u: Gridfunction object to be copied.
 * v: Gridfunction to be updated.
 * notes: v will be overwritten by u.
 *        Only the values of the inner grid points are taken into account.
 *        Memory for the copy is supposed to already be allocated. */
void copy_gridfunction(pcgridfunction u, pgridfunction v);

/* Initialize a gridfunction object u with the function values of the
 * inner points.
 * f: Continuous counterpart of the grid function.
 * u: Gridfunction object to be initialized.
 * notes: u will be modified in row-major order. */
void init_gridfunction(function f, pgridfunction u);

/* Initialize a gridfunction object u with the function values of the
 * boundary points only.
 * f: Continuous counterpart of the grid function.
 * u: Gridfunction object to be initialized. */
void boundary_gridfunction(function f, pgridfunction u);

/****************************
In- and output
****************************/

/* Save the a grid function into a file.
 * u: Respective gridfunction object.
 * filename: Name of the file where the grid function shall be saved.
 * notes: An output file of the given name will be produced.
          It will contain only the x- and y-coordinates of all inner
          and boundary gridpoints and corresponding values of the grid
          function with no additional data. */
void saveToFile_gridfunction(pcgridfunction u, char *filename);

/* Write a gridfunction to file "y/wave2d_y_xxxxx.dat",
 * where 'xxxxx' is given by the current frame and 'y' by the iteration used.
 * u: Respective gridfunction object.
 * iteration: Name of the iteration in use.
 * frame: Current frame. */
void output_gridfunction(pcgridfunction u, char *iteration, int frame);

/****************************
Basic linear algebra routines
****************************/

/* Compute a linear combination of two grid functions by scaling one and
 * adding it to the other one.
 * a: Scaling factor for u.
 * u: Gridfunction object to be added.
 * v: Gridfunction object to be updated.
 * notes: Only the values of the inner grid points are taken into account.*/
void add_gridfunction(double a, pcgridfunction u, pgridfunction v);

/* Evaluate the negative discrete Laplacian on a grid function,
 * scale the result by a certain factor and add it to another grid function.
 * a: Scaling factor.
 * u: Gridfunction object that the Laplacian acts on.
 * v: Resulting gridfunction object.
 * notes: The boundary values will not be modified.
 *        Memory for the result is supposed to already be allocated. */
void addevalLaplace_gridfunction(double a, void *data, pcgridfunction u, pgridfunction v);

/* Scale all inner values of a grid function by a. */
void scal_gridfunction(double a, pgridfunction u);

/* Set all inner values of a grid function to zero.*/
void clear_gridfunction(pgridfunction u);

/* Compute the scalar product of two gridfunctions. */
double dotprod_gridfunction(pcgridfunction u, pcgridfunction v);

/* Compute the norm of a gridfunction. */
double norm_gridfunction(pcgridfunction u);

/*****************************
 * Callback functions
 *****************************/

void(add_gridfunc)(double a, void *x, void *y);
void(scal_gridfunc)(double a, void *x);
void(clear_gridfunc)(void *x);
double(dotprod_gridfunc)(void *x, void *y);
double(norm_gridfunc)(void *x);
#endif