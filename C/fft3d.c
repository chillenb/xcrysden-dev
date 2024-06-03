#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

/*#include "shankland.h"*/
#include "fft3d.h"
#include "memory.h"
#include "tensor.h"

/*
  BEWARE: ***src is a tensor3f of general grid (despite it contain
  periodic function), hence we need to the following pipeline:

  1. transform general-grid (tensor3f) --> periodic-grid (double*) 
  2. interpolate
  3. transform periodic-grid (double*) to general-grid (tensor3f)
*/

float*** general_grid_fft_interpolator_tensor3f(int n[3], int degree[3], float ***src)
{
  int i, j, k, in, ijn;

  int n1 = n[0]-1;
  int n2 = n[1]-1;
  int n3 = n[2]-1;

  int n1new = n1*degree[0];
  int n2new = n2*degree[1];
  int n3new = n3*degree[2];

  float ***result;
  double *out, *func;

  MALLOC_TENSOR3(float, result, n1new+1, n2new+1, n3new+1); 
  
  func = (double*)malloc(sizeof(double) * n1*n2*n3 );

  for(i=0; i<n1; i++) 
    {
      in = i*n2;
      for(j=0; j<n2; j++) 
	{
	  ijn = (in+j)*n3;

	  for(k=0; k<n3; k++)
	    func[ijn+k] = src[i][j][k];
	}
    }

  out = fft3d_interpolator(n1, n2, n3, n1new, n2new, n3new, func);

  /* take care of redundant points */
  for(i=0; i<n1new+1; i++) 
    {
      in = (i%n1new)*n2new;
      for(j=0; j<n2new+1; j++) 
	{	  
	  ijn = (in+(j%n2new))*n3new;
	  
	  for(k=0; k<n3new+1; k++)
	    result[i][j][k] = out[ijn+(k%n3new)];
	}
    }

  free(out);
  free(func);

  return result;
}


double*
fft3d_interpolator(int n1, int n2, int n3, int n1new, int n2new, int n3new, double *func)
{
  double *result = (double*) calloc(n1new*n2new*n3new, sizeof(double));
  
  fft3d_interpolate(n1, n2, n3, n1new, n2new, n3new, func, result);
  
  return result;
}

void
fft3d_interpolate(int n1, int n2, int n3, int n1new, int n2new, int n3new, double *func, double *result)
{
  double inv_omega;
  int i, dim3, dim3new;
  
  fftw_complex *aux, *fftw2;
  fftw_plan    pfor, pbak;

  inv_omega = 1.0 / (double) (n1*n2*n3);

  dim3    = n3/2+1;
  dim3new = n3new/2+1;

  aux    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1*n2*dim3);
  fftw2  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1new*n2new*dim3new);

  pfor = fftw_plan_dft_r2c_3d(n1, n2, n3, func, aux, FFTW_ESTIMATE);
  pbak = fftw_plan_dft_c2r_3d(n1new, n2new, n3new, fftw2, result, FFTW_ESTIMATE);

  /*
    FFTW_FORWARD
  */
 
  fftw_execute(pfor);

  for(i=0; i<n1*n2*dim3; i++) COMPLEX_SCALE(aux[i], inv_omega);

  fft3d_enlarge_complex_r2c(n1, n2, n3, aux, n1new, n2new, n3new, fftw2);

  /*
    FFTW_BACKWARD
  */

  fftw_execute(pbak);

  /*
    destroy FFTW plans
  */
  fftw_destroy_plan(pfor);
  fftw_destroy_plan(pbak);

  /*
    deallocate memmory
  */
  
  fftw_free((void*) fftw2);
  fftw_free((void*) aux);
}


