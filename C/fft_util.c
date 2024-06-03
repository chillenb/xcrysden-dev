#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

/*#include "shankland.h"*/
#include "fft3d.h"

void
fft2d_enlarge_complex_r2c(size_t n1,    size_t n2,    fftw_complex *in, 
			  size_t n1new, size_t n2new, fftw_complex *out) 
{
  int n1min, n1max, n2max, dim2, dim2new, i, j, i_ind, ii_ind;

  if (in == out) {
    fprintf(stderr, "from fft2d_enlarge_complex_r2c: inplace enlargement is not supported, aborting ...");
    abort();
  }

  dim2    = n2/2+1;
  dim2new = n2new/2+1;
  
  n1min = FT_MIN(n1);
  n1max = FT_MAX(n1);
  n2max = FT_MAX(n2);
  
  /* clear the *out */
  for(i=0; i<n1new*dim2new; i++) COMPLEX_SET(out[i], 0.0, 0.0);

  for(i=n1min; i<=n1max; i++)
    {
      ii_ind = FT_INDEX(i,n1new)*dim2new;
      i_ind  = FT_INDEX(i,n1)*dim2;
      
      for(j=0; j<=n2max; j++)
      	{
      	  COMPLEX_COPY ( out[ii_ind + j], in[i_ind + j] );
      	}
    }

  /* account for double counting in 2nd dim if n2 is even */

  if ( n2%2 == 0 )
    {
      for(i=n1min; i<=n1max; i++)
	{
	  ii_ind = FT_INDEX(i,n1new)*dim2new;
	  COMPLEX_SCALE(out[ii_ind + n2max], 0.5);
	}
    }
}

void
fft3d_enlarge_complex_r2c(size_t n1,    size_t n2,    size_t n3,    fftw_complex *in, 
			  size_t n1new, size_t n2new, size_t n3new, fftw_complex *out) 
{
  int n1min, n1max, n2min, n2max, n3max, dim3, dim3new, i, j, k, i_ind, ii_ind, j_ind, jj_ind;

  if (in == out) {
    fprintf(stderr, "from fft3d_enlarge_complex_r2c: inplace enlargement is not supported, aborting ...");
    abort();
  }

  dim3    = n3/2+1;
  dim3new = n3new/2+1;
  
  n1min = FT_MIN(n1);
  n1max = FT_MAX(n1);
  n2min = FT_MIN(n2);
  n2max = FT_MAX(n2);
  n3max = FT_MAX(n3);
  
  /* clear the *out */
  for(i=0; i<n1new*n2new*dim3new; i++) COMPLEX_SET(out[i], 0.0, 0.0);

  for(i=n1min; i<=n1max; i++)
    {
      ii_ind = FT_INDEX(i,n1new)*n2new;
      i_ind  = FT_INDEX(i,n1)*n2;
      
      for(j=n2min; j<=n2max; j++)
	{
	  jj_ind = (ii_ind + FT_INDEX(j,n2new)) * dim3new;
	  j_ind  = (i_ind + FT_INDEX(j,n2)) * dim3;

	  for(k=0; k<=n3max; k++)
	    {
	      COMPLEX_COPY ( out[jj_ind + k], in[j_ind + k] );
	    }
	}
    }

  /* account for double counting in 3rd dim if n3 is even */

  if ( n3%2 == 0 )
    {
      for(i=n1min; i<=n1max; i++)
	{
	  ii_ind = FT_INDEX(i,n1new)*n2new;
      
	  for(j=n2min; j<=n2max; j++)
	    {
	      jj_ind = (ii_ind + FT_INDEX(j,n2new)) * dim3new;
	      COMPLEX_SCALE(out[jj_ind + n3max], 0.5);
	    }
	}
    }
}

