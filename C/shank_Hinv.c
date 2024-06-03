#ifdef XC_WITH_SHANKLAND

#include <math.h>
#include <stdlib.h>
#include <meschach/matrix2.h>

#include "shankland.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif
const double two_pi = 2.0*M_PI;

/*
  construct the H^-1 Shankland matrix

*/

MAT*
shank_HinvObj(size_t Npoints, 
	      size_t Nwaves, 
	      double sigma, 
	      double (*roughnessFunc)(double m2, double sigma))
{
  MAT *Hinv = m_get(Npoints,Npoints);
  shank_Hinv(Npoints, Nwaves, sigma, roughnessFunc, Hinv);
  return Hinv;
}

void
shank_Hinv(size_t Npoints,
	   size_t Nwaves, 
	   double sigma, 
	   double (*roughnessFunc)(double m2, double sigma), 
	   MAT *Hinv)
{
  size_t i, j;
  int m, m_min, m_max;
  double inv_Npoints, inv_Nwaves, tpi_bn, tpi_ibn, scale;

  double *Hi_j = (double*) malloc(sizeof(double)*Npoints);
  MAT *H       = m_get(Npoints,Npoints);

  /* new experimentation */
  /* if (Nwaves < 100 && sigma != 0.0) { */
  /*   Nwaves = Npoints * (100/Npoints + 1); */
  /* } */
  /* END */

  inv_Npoints = 1.0 / (double) Npoints;
  inv_Nwaves  = 1.0 / (double) Nwaves;
  tpi_bn      = two_pi * inv_Npoints;
  scale       = (double)Npoints/(double)Nwaves;

  m_min = FT_MIN(Nwaves);
  m_max = FT_MAX(Nwaves);

#ifdef DEBUG
  printf("m=%d, m_min=%d, m_max=%d\n", Nwaves, m_min, m_max);
  printf("rougness factor from: %e to %e\n",
	 1.0/roughnessFunc(1.0,sigma), 
	 1.0/roughnessFunc(m_max*m_max,sigma));
#endif

  for (i=0; i<Npoints; i++)
    {
      Hi_j[i] = 0.0;
      tpi_ibn = tpi_bn * (double) i;
      
      for (m=m_min; m<=m_max; m++)
	{	  
	  /* if (m==0 && sigma == 0.0 ) continue;  */
	  Hi_j[i] += cos( tpi_ibn*m ) / roughnessFunc(m*m,sigma);
	} 
      Hi_j[i] *=  inv_Npoints;
    }  

  for (i=0; i<Npoints; i++) 
    {
      H->me[i][i] = Hi_j[0];
      
      for (j=i+1; j<Npoints; j++) 
	H->me[i][j] = H->me[j][i] = Hi_j[j-i];	
    }


#ifdef DEBUG
  print_MAT_matrix(stdout,"H", H);
#endif

  m_inverse(H, Hinv);

#ifdef DEBUG
  print_MAT_matrix(stdout,"Hinv", Hinv);
#endif

  m_free(H);
  free((void*)Hi_j);
}


double*
shank1d_lambdaObj(int n1, int n1new, double sigma, double* func, double (*roughnessFunc)(double m2, double sigma), double *f_aver)
{
  double* lambda = (double*) malloc(sizeof(double) * n1);
  *f_aver = shank1d_lambda(n1, n1new, sigma, func, lambda, roughnessFunc);
  return lambda;
} 

double
shank1d_lambda(int n1, int n1new, double sigma, double* func, double* lam, double (*roughnessFunc)(double m2, double sigma))
{
  int i1, j1;
  MAT *Hinv1;
  double f_aver = 0.0;

#ifdef DEBUG
  double up = 0.0;
  double down = 0.0;
  double sum_lam = 0.0;
#endif

  Hinv1 = shank_HinvObj(n1, n1new, sigma, roughnessFunc);

  /* calculate f_aver */

  for (i1=0; i1<n1; i1++)
    {
      f_aver += func[i1];
    }
  f_aver /= (double) n1;

#ifdef DEBUG    
  printf("f_aver = %f\n", f_aver);
#endif

  /*
    calculate lambda matrix
  */

  for (i1=0; i1<n1; i1++) 
    {
      lam[i1] = 0.0;	
	
      for (j1=0; j1<n1; j1++) 
	{
	  lam[i1] += Hinv1->me[i1][j1] * (func[j1] - f_aver);
#ifdef DEBUG	      
	  up   += Hinv1->me[i1][j1] * func[j1];
	  down += Hinv1->me[i1][j1];
#endif		
	}
    } 

#ifdef DEBUG
  printf("f_aver - f_mean = %e - %e = %e\n", f_aver, up/down, f_aver - up/down);
  print_double_vector(stdout,"lambda",n1, lam);

  for (i1=0; i1<n1; i1++)
    sum_lam += lam[i1];

  printf("Sum_i lambda_i = %f\n", sum_lam);  
#endif
  
  m_free((void*)Hinv1);

  return f_aver;
}


double*
shank2d_lambdaObj(int n1, int n2, int n1new, int n2new, double sigma, double* func, double (*roughnessFunc)(double m2, double sigma), double *f_aver)
{
  double* lambda = (double*) malloc(sizeof(double) * n1*n2);
  *f_aver = shank2d_lambda(n1, n2, n1new, n2new, sigma, func, lambda, roughnessFunc);
  return lambda;
} 

double
shank2d_lambda(int n1, int n2, int n1new, int n2new, double sigma, double* func, double* lam, double (*roughnessFunc)(double m2, double sigma))
{
  int i1,i2, j1,j2, indL, j1n;
  MAT *Hinv1;
  MAT *Hinv2;
  double f_aver = 0.0;

#ifdef DEBUG
  double up = 0.0;
  double down = 0.0;
  double sum_lam = 0.0;
#endif

  Hinv1 = shank_HinvObj(n1, n1new, sigma, roughnessFunc);
  if (n1 == n2) 
    Hinv2 = Hinv1;
  else
    Hinv2 = shank_HinvObj(n2, n2new, sigma, roughnessFunc);

  /* calculate f_aver */

  for (i1=0; i1<n1; i1++)
    {
      j1n = i1*n2;
      for (i2=0; i2<n2; i2++) f_aver += func[j1n+i2];
    }
  f_aver /= (double) (n1*n2);

#ifdef DEBUG    
  printf("f_aver = %f\n", f_aver);
#endif

  /*
    calculate lambda matrix
  */

  for (i1=0; i1<n1; i1++) 
    for (i2=0; i2<n2; i2++) 
      {
	indL = i1*n2 + i2;
	lam[indL] = 0.0;	
	
	for (j1=0; j1<n1; j1++) 
	  {
	    j1n = j1*n2;
	    
	    for (j2=0; j2<n2; j2++) {
	      lam[indL] += Hinv1->me[i1][j1] * Hinv2->me[i2][j2] * (func[j1n+j2] - f_aver);
#ifdef DEBUG	      
	      up += Hinv1->me[i1][j1] * Hinv2->me[i2][j2] * func[j1n+j2];
	      down += Hinv1->me[i1][j1] * Hinv2->me[i2][j2];
#endif		
	    }
	  } 
      }

#ifdef DEBUG
  printf("f_aver - f_mean = %e - %e = %e\n", f_aver, up/down, f_aver - up/down);
  print_double_matrix(stdout,"lambda",n1, n2, lam);

  for (i1=0; i1<n1*n2; i1++)
    sum_lam += lam[i1];

  printf("Sum_i lambda_i = %f\n", sum_lam);  
#endif
  
  if (Hinv2 != Hinv1) 
    m_free((void*)Hinv2);
  
  m_free((void*)Hinv1);

  return f_aver;
}


double*
shank3d_lambdaObj(int n1, int n2, int n3, int n1new, int n2new, int n3new, double sigma, double* func, double (*roughnessFunc)(double m2, double sigma), double *f_aver)
{
  double* lambda = (double*) malloc(sizeof(double) * n1*n2*n3);
  *f_aver = shank3d_lambda(n1, n2, n3, n1new, n2new, n3new, sigma, func, lambda, roughnessFunc);
  return lambda;
} 

double
shank3d_lambda(int n1, int n2, int n3, int n1new, int n2new, int n3new, double sigma, double* func, double* lam, double (*roughnessFunc)(double m2, double sigma))
{
  int i1,i2,i3, j1,j2,j3, indL, i1n, i12n, j1n, j12n;
  MAT *Hinv1;
  MAT *Hinv2;
  MAT *Hinv3;
  double H12, f_aver = 0.0;

#ifdef DEBUG
  double up = 0.0;
  double down = 0.0;
  double sum_lam = 0.0;
#endif

  Hinv1 = shank_HinvObj(n1, n1new, sigma, roughnessFunc);

  if (n2 == n1) Hinv2 = Hinv1;
  else          Hinv2 = shank_HinvObj(n2, n2new, sigma, roughnessFunc);
  
  if (n3 == n1)      Hinv3 = Hinv1;
  else if (n3 == n2) Hinv3 = Hinv2;
  else               Hinv3 = shank_HinvObj(n3, n3new, sigma, roughnessFunc);

  /* calculate f_aver */

  for (i1=0; i1<n1; i1++)
    {
      j1n = i1*n2;
      
      for (i2=0; i2<n2; i2++) 
	{
	  j12n = (j1n + i2)*n3;
	  
	  for (i3=0; i3<n3; i3++) f_aver += func[j12n+i3];
	}
    }
  f_aver /= (double) (n1*n2*n3);

#ifdef DEBUG    
  printf("f_aver = %f\n", f_aver);
#endif

  /*
    calculate lambda matrix
  */

  for (i1=0; i1<n1; i1++) 
    {
      i1n = i1*n2;
      
      for (i2=0; i2<n2; i2++) 
	{
	  i12n = (i1n + i2)*n3;
	  
	  for (i3=0; i3<n3; i3++) 
	    {
	      indL = i12n + i3;
	      lam[indL] = 0.0;	
	    
	      for (j1=0; j1<n1; j1++) 
		{
		  j1n = j1*n2;
		  
		  for (j2=0; j2<n2; j2++) 
		    {
		      j12n = (j1n + j2)*n3;
		      H12 = Hinv1->me[i1][j1] * Hinv2->me[i2][j2];

		      for (j3=0; j3<n3; j3++) 			
			{
			  lam[indL] += H12 * Hinv3->me[i3][j3] * (func[j12n+j3] - f_aver);
#ifdef DEBUG	      
			  up += H12 * Hinv3->me[i3][j3] * func[j12n+j3];
			  down += H12 * Hinv3->me[i3][j3];
#endif		
			}
		    }
		}
	    }
	}
    }
  
#ifdef DEBUG
  printf("f_aver - f_mean = %e - %e = %e\n", f_aver, up/down, f_aver - up/down);
  print_double_vector(stdout,"lambda",n1*n2*n3, lam); 

  for (i1=0; i1<n1*n2*n3; i1++)
    sum_lam += lam[i1];

  printf("Sum_i lambda_i = %f\n", sum_lam);  
#endif
  
  if (Hinv3 != Hinv2 && Hinv3 != Hinv1) 
    m_free((void*)Hinv3);

  if (Hinv2 != Hinv1) 
    m_free((void*)Hinv2);
  
  m_free((void*)Hinv1);
  
  return f_aver;
}

#else
typedef int make_iso_compilers_happy;
#endif
