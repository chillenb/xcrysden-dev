#ifdef XC_WITH_SHANKLAND

#include <stdio.h>
#include <math.h>
#include <meschach/matrix2.h>
#include <fftw3.h>

#include "shankland.h"

extern const double two_pi;

/*
  NOTICE: m2 stands for m*m
*/

double one_func(double m2, double sigma)
{
  return 1.0;
}

double first_order(double m2, double sigma)
{
  return (1.0 + sigma*m2);
}

double first_sqrt_order(double m2, double sigma)
{
  m2=sigma*m2;
  return (1.0 + m2*(1.0 + 0.01*sqrt(m2)));
}

double second_order(double m2, double sigma)
{
  m2=sigma*m2;
  return (1.0 + m2*(1.0 + 0.001*m2));
}

double pickett(double m2, double sigma)
{
  double xm;

  m2 = sigma*m2;
  xm = (1.0 - m2);
  return (xm*xm + m2*m2*m2);
}

double gauss(double m2, double sigma)
{
  return (exp(sigma*m2));
}

double exp_abs(double m2, double sigma)
{
  return (exp(sigma*sqrt(m2)));
}


double estimate_sigma(double power, int m, double (*roughnessFunc)(double m2, double sigma))
{
  const int max=20;
  int i, ind;
  double sigma, d, diff, max_num = pow(10.0, power), min=10e16;

  m = FT_MAX(m);
  sigma = 1.0;
  
  if ( roughnessFunc(m*m, sigma) > max_num) 
    {
      d = 0.9;
      for(i=0; i<2000; i++) 
	{
	  if ( roughnessFunc(m*m, sigma) < max_num )
	    {
	      return ( sigma*1.05263 );
	    }
	  sigma *= d;	  
	}
    }
  else /* roughnessFunc(m*m, sigma) < max_num */
    {
      d = 1.1;
      for(i=0; i<2000; i++) 
	{
	  if ( roughnessFunc(m*m,sigma) > max_num )
	    {
	      return ( sigma*.95238 );
	    }
	  sigma *= d;	  
	}
    }

  /* if we arrived here, sigma was not found */

  ind = -max;
  for (i=-max; i<=max; i++)
    {
      sigma = pow(10.0,i);
      diff = ABS( 1.0/max_num - 1.0/roughnessFunc(m*m, sigma) );
     
      if ( diff < min ) {
	ind = i;
	min = diff;
      }
    }
  
  return ( pow(10.0,ind) );
}


void
print_double_vector(FILE *stream, const char* name, size_t m, double *mat)
{
  int i;
  fprintf(stream,"%s == \n", name);
  for (i=0; i<m; i++) 
    fprintf(stream,"| %21.14e |\n", mat[i]);
  fprintf(stream,"\n");
}

#ifdef DEBUG
void
print_complex_vector(FILE *stream, const char* name, size_t m, fftw_complex *mat)
{
  int i;
  fprintf(stream,"%s == \n", name);
  for (i=0; i<m; i++)
    fprintf(stream,"| (%21.13e,%21.13e) |\n", mat[i][0],mat[i][1]);
  fprintf(stream,"\n");
}
#endif


void
print_MAT_matrix(FILE *stream, const char* name, MAT *mat)
{
  int i, j;
  fprintf(stream,"%s == \n", name);
  for (i=0; i<mat->m; i++) {
    fprintf(stream,"| ");
    for (j=0; j<mat->n; j++)
      fprintf(stream,"%21.13e ", mat->me[i][j]);
    fprintf(stream,"|\n");
  }
}

void
print_double_matrix(FILE *stream, const char* name, size_t m, size_t n, double *mat)
{
  int i, j;
  fprintf(stream,"%s == \n", name);
  for (i=0; i<m; i++) {
    fprintf(stream,"| ");
    for (j=0; j<n; j++)
      fprintf(stream,"%21.13e ", mat[i*n+j]);
    fprintf(stream,"|\n");
  }
}

#ifdef DEBUG
void
print_complex_matrix(FILE *stream, const char* name, size_t m, size_t n, fftw_complex *mat)
{
  int i, j;
  fprintf(stream,"%s == \n", name);
  for (i=0; i<m; i++) {
    fprintf(stream,"| ");
    for (j=0; j<n; j++)
      fprintf(stream,"(%21.13e,%21.13e) ", mat[i*n+j][0],mat[i*n+j][1]);
    fprintf(stream,"|\n");
  }
}
#endif

#else
typedef int make_iso_compilers_happy;
#endif
