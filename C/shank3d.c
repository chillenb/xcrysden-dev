#ifdef XC_WITH_SHANKLAND

#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <meschach/matrix2.h>

#include "shankland.h"
#include "memory.h"
#include "tensor.h"

float*** general_grid_shankland_interpolator_tensor3f(int n[3], int degree[3], float ***src);
double *shankland3d_interpolator(int n1, int n2, int n3, int n1new, int n2new, int n3new, double *func, double (*roughnessFunc)(double m2, double sigma), double log10_damp);
void shankland3d_interpolate(int n1, int n2, int n3, int n1new, int n2new, int n3new, double *func, double (*roughnessFunc)(double m2, double sigma), double log10_damp, double *result);

/*
  BEWARE: ***src is a tensor3f of general grid (despite it contain
  periodic function), hence we need to the following pipeline:

  1. transform general-grid (tensor3f) --> periodic-grid (double*) 
  2. interpolate
  3. transform periodic-grid (double*) to general-grid (tensor3f)
*/

float*** general_grid_shankland_interpolator_tensor3f(int n[3], int degree[3], float ***src)
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

  fprintf(stderr, "\n *** executing general_grid_shankland_interpolator_tensor3f ***\n");
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

  out = shankland3d_interpolator(n1, n2, n3, n1new, n2new, n3new, func, gauss, 6.0);

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
shankland3d_interpolator(int n1, int n2, int n3, int n1new, int n2new, int n3new, double *func, double (*roughnessFunc)(double m2, double sigma), double log10_damp)
{
  double *result  = (double*) calloc(n1new*n2new*n3new, sizeof(double));

  shankland3d_interpolate(n1, n2, n3, n1new, n2new, n3new, func, roughnessFunc, log10_damp, result);

  return result;
}


void
shankland3d_interpolate(int n1, int n2, int n3, int n1new, int n2new, int n3new, double *func, double (*roughnessFunc)(double m2, double sigma), double log10_damp, double *result)
{
  double sigma, f_aver, inv_omega, inv_rough_i, inv_rough_ij, fac;
  int n1min, n1max, n2min, n2max, n3min, n3max, i, j, k, i_ind, j_ind, k_ind, i0_ind, j0_ind, dim3;
  
  double       *lambda;
  fftw_complex *fftw2;
  fftw_plan    pfor, pbak;

  lambda = (double*) calloc(n1*n2*n3, sizeof(double));
  fftw2  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n1new*n2new*n3new);

  pfor = fftw_plan_dft_r2c_3d(n1new, n2new, n3new, result, fftw2, FFTW_ESTIMATE);
  pbak = fftw_plan_dft_c2r_3d(n1new, n2new, n3new, fftw2, result, FFTW_ESTIMATE);

  /*
    estimate sigma
  */ 
  if ( n1==n1new && n2==n2new && n3==n3new ) {
    sigma = 0.0;
  } else {
    sigma = estimate_sigma(log10_damp, 2*MAX(n1,MAX(n2,n3)), roughnessFunc);  
  }

#ifdef DEBUG
  printf("sigma check: log10_damp = %f, 1.0/roughnessFunc(%e*mmax^2) = %e\n", sigma, 
	 1.0/roughnessFunc(MAX(n1,MAX(n2,n3))*MAX(n1,MAX(n2,n3)),sigma));
#endif


  /*
    calc lambda matrix
  */
    
  f_aver = shank3d_lambda(n1, n2, n3, n1new, n2new, n3new, sigma, func, lambda, roughnessFunc);

  /* 
     lambda ==> result:
     ----------------
     fill the zeroes into result so as to match the [n1new][n2new][n3new] dimensions
  */

  for(i=0; i<n1; i++)
    {
      i0_ind = i*n2*n3;
      i_ind  = i*(n1new/n1)*n2new*n3new;

      for (j=0; j<n2; j++)
  	{
	  j0_ind = j*n3;
	  j_ind = j*(n2new/n2)*n3new;
	  
	  for (k=0; k<n3; k++)
	    {
	      /* _result[i*(n1new/n1)][j*(n2new/n2)][k*(n3new/n3)] := _lambda[i][j][k]; */
	      result[i_ind + j_ind + k*(n3new/n3)] = lambda[i0_ind + j0_ind + k];
	    }
	}
    }

  /*
    FFTW_FORWARD
  */
 
  fftw_execute(pfor);

  dim3 = n3new/2 + 1;

  /*
    adopt/calculate am's
  */
  
  n1min = FT_MIN(n1new);
  n1max = FT_MAX(n1new);
  
  n2min = FT_MIN(n2new);
  n2max = FT_MAX(n2new);
  
  n3min = FT_MIN(n3new);
  n3max = FT_MAX(n3new);
  
  inv_omega = 1.0 / (double) (n1*n2*n3);

  for(i=n1min; i<=n1max; i++)
    {
      i_ind  = FT_INDEX(i,n1new)*n2new*dim3;      
      inv_rough_i = inv_omega / roughnessFunc((double)(i*i), sigma);

      for(j=n2min; j<=n2max; j++)
	{
	  j_ind = FT_INDEX(j,n2new)*dim3;
	  inv_rough_ij = inv_rough_i / roughnessFunc((double)(j*j), sigma);
	  
	  for (k=0; k<=n3max; k++) 
	    {
	      fac = inv_rough_ij / roughnessFunc((double)(k*k), sigma);

	      fftw2[i_ind + j_ind + k][0] *= fac;
	      fftw2[i_ind + j_ind + k][1] *= fac;
	    }
	}
    }
  COMPLEX_SET(fftw2[0], 0.0, 0.0);

  /*
    FFTW_BACKWARD
  */

  fftw_execute(pbak);

 /* 
     add the f_aver to result
  */
  for (i=0; i<n1new*n2new*n3new; i++) result[i] += f_aver;

  /*
    destroy FFTW plans
  */
  fftw_destroy_plan(pfor);
  fftw_destroy_plan(pbak);

  /*
    deallocate memmory
  */
  
  fftw_free((void*) fftw2);
  free((void*) lambda);  
}
#else
typedef int make_iso_compilers_happy;
#endif
