/*
 
 *****************************************************************************
 * Author:                                                                   *
 * ------                                                                    *
 *  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  *
 *  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    *
 *  Jozef Stefan Institute                          Fax: x 386 1 477 3811    *
 *  Jamova 39, SI-1000 Ljubljana                                             *
 *  SLOVENIA                                                                 *
 *                                                                           *
 * Source: $XCRYSDEN_TOPDIR/C/linear.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <stdio.h>
#include <stdlib.h>
#include "tensor.h"
#include "memory.h"

/* linear.c */
float ***cryGeneralGridRegPeriodInterpolator_f_Linear3D(int n[3], int degree[3], float ***src);

float ***cryGeneralGridRegPeriodInterpolator_f_Linear3D(int n[3],      /*!< number of input points of general grid */
							int degree[3], /*!< degree of interp   */
							float ***src)  /*!< INPUT general grid */
{
  register int i0, i1, j0, j1, k0, k1, ii, jj;
  int n0, n1, n2, nn[3];
  float ay0, ay1, az0, az1, fr[3];
  float ***dst;

  nn[0] = n[0]-1;
  nn[1] = n[1]-1;
  nn[2] = n[2]-1;

  n0 =  nn[0] * degree[0];
  n1 =  nn[1] * degree[1];
  n2 =  nn[2] * degree[2];
  MALLOC_TENSOR3(float, dst, n0+1, n1+1, n2+1);  

  /* 3D case */
  for (i0=0; i0<nn[0]; i0++)    
    for  (i1=0; i1<degree[0]; i1++) 
      {
	fr[0] = (float)i1 / (float)degree[0];	
	ii    = i0*degree[0] + i1;

	for (j0=0; j0<nn[1]; j0++)
	  for  (j1=0; j1<degree[1]; j1++) 
	    {
	      fr[1] = (float)j1 / (float)degree[1];	
	      jj    = j0*degree[1] +j1;

	      for (k0=0; k0<nn[2]; k0++)
		for  (k1=0; k1<degree[2]; k1++) 
		  {
		    fr[2] = (float)k1 / (float)degree[2];	
	      	      
		    /* Z0 */
		    ay0 = (1.0-fr[0])*src[i0][j0][k0]   + fr[0]*src[i0+1][j0][k0];
		    ay1 = (1.0-fr[0])*src[i0][j0+1][k0] + fr[0]*src[i0+1][j0+1][k0];
		    /**/
		    az0 = (1.0-fr[1])*ay0 + fr[1]*ay1;

		    /* Z1 */
		    ay0 = (1.0-fr[0])*src[i0][j0][k0+1]   + fr[0]*src[i0+1][j0][k0+1];
		    ay1 = (1.0-fr[0])*src[i0][j0+1][k0+1] + fr[0]*src[i0+1][j0+1][k0+1];
		    /**/
		    az1 = (1.0-fr[1])*ay0 + fr[1]*ay1;

		    dst[ii][jj][k0*degree[2] + k1] = (1.0-fr[2])*az0 + fr[2]*az1;		    
		  }
	    }
      }

  /* assing the plane-border points: xy, xz, yz */
  for (i0=0; i0<n0; i0++)
    for (j0=0; j0<n1; j0++)
      dst[i0][j0][n2] = dst[i0][j0][0];
  
  for (i0=0; i0<n0; i0++)
    for (k0=0; k0<n2; k0++)
      dst[i0][n1][k0] = dst[i0][0][k0];
  
  for (j0=0; j0<n1; j0++)
    for (k0=0; k0<n2; k0++)
      dst[n0][j0][k0] = dst[0][j0][k0];
  
  /* assing the line-border points x, y, z*/
  for (i0=0; i0<n0; i0++)
    dst[i0][n1][n2] = dst[i0][0][0];
  
  for (j0=0; j0<n1; j0++)
    dst[n0][j0][n2] = dst[0][j0][0];
  
  for (k0=0; k0<n2; k0++)
    dst[n0][n1][k0] = dst[0][0][k0];
  
  /* the last point */
  dst[n0][n1][n2] = dst[0][0][0];

  return dst;
}



/* 
   BEWARE: this interpolation is NOT YET IMPLEMENTED !!!

   Below function is merely a copy of cryGeneralGridRegPeriodInterpolator_f_Linear3D !!!
   
*/
float ***cryGeneralGridRegInterpolator_f_Linear3D(int n[3],      /*!< number of input points of general grid */
						  int degree[3], /*!< degree of interp   */
						  float ***src)  /*!< INPUT general grid */
{
  register int i0, i1, j0, j1, k0, k1, ii, jj;
  int n0, n1, n2, nn[3];
  float ay0, ay1, az0, az1, fr[3];
  float ***dst;

  nn[0] = n[0]-1;
  nn[1] = n[1]-1;
  nn[2] = n[2]-1;

  n0 =  nn[0] * degree[0];
  n1 =  nn[1] * degree[1];
  n2 =  nn[2] * degree[2];
  MALLOC_TENSOR3(float, dst, n0+1, n1+1, n2+1);  

  /* 3D case */
  for (i0=0; i0<nn[0]; i0++)    
    for  (i1=0; i1<degree[0]; i1++) 
      {
	fr[0] = (float)i1 / (float)degree[0];	
	ii    = i0*degree[0] + i1;

	for (j0=0; j0<nn[1]; j0++)
	  for  (j1=0; j1<degree[1]; j1++) 
	    {
	      fr[1] = (float)j1 / (float)degree[1];	
	      jj    = j0*degree[1] +j1;

	      for (k0=0; k0<nn[2]; k0++)
		for  (k1=0; k1<degree[2]; k1++) 
		  {
		    fr[2] = (float)k1 / (float)degree[2];	
	      	      
		    /* Z0 */
		    ay0 = (1.0-fr[0])*src[i0][j0][k0]   + fr[0]*src[i0+1][j0][k0];
		    ay1 = (1.0-fr[0])*src[i0][j0+1][k0] + fr[0]*src[i0+1][j0+1][k0];
		    /**/
		    az0 = (1.0-fr[1])*ay0 + fr[1]*ay1;

		    /* Z1 */
		    ay0 = (1.0-fr[0])*src[i0][j0][k0+1]   + fr[0]*src[i0+1][j0][k0+1];
		    ay1 = (1.0-fr[0])*src[i0][j0+1][k0+1] + fr[0]*src[i0+1][j0+1][k0+1];
		    /**/
		    az1 = (1.0-fr[1])*ay0 + fr[1]*ay1;

		    dst[ii][jj][k0*degree[2] + k1] = (1.0-fr[2])*az0 + fr[2]*az1;		    
		  }
	    }
      }

  /* assing the plane-border points: xy, xz, yz */
  for (i0=0; i0<n0; i0++)
    for (j0=0; j0<n1; j0++)
      dst[i0][j0][n2] = dst[i0][j0][0];
  
  for (i0=0; i0<n0; i0++)
    for (k0=0; k0<n2; k0++)
      dst[i0][n1][k0] = dst[i0][0][k0];
  
  for (j0=0; j0<n1; j0++)
    for (k0=0; k0<n2; k0++)
      dst[n0][j0][k0] = dst[0][j0][k0];
  
  /* assing the line-border points x, y, z*/
  for (i0=0; i0<n0; i0++)
    dst[i0][n1][n2] = dst[i0][0][0];
  
  for (j0=0; j0<n1; j0++)
    dst[n0][j0][n2] = dst[0][j0][0];
  
  for (k0=0; k0<n2; k0++)
    dst[n0][n1][k0] = dst[0][0][k0];
  
  /* the last point */
  dst[n0][n1][n2] = dst[0][0][0];

  return dst;
}


/*   /\* 1D case *\/ */
/*   for (i0=0; i0<n[0]; i0++) */
/*     { */
/*       for (i1=0; i1<degree[0]; i1++)  */
/* 	{ */
/* 	  fr = (float)i1 / (float)degree[0];	 */
/* 	  dst[i0*degree[0] + i1] = (1.0-fr)*src[i0] + fr*src[i0+1]; */
/* 	} */
/*     } */

/*   /\* 2D case *\/ */
/*   for (i0=0; i0<n[0]; i0++)     */
/*     for  (i1=0; i1<degree[0]; i1++)  */
/*       { */
/* 	fr[0] = (float)i1 / (float)degree[0];	 */
	
/* 	for (j0=0; j0<n[1]; j0++) */
/* 	  for  (j1=0; j1<degree[1]; j1++)  */
/* 	    { */
/* 	      fr[1] = (float)j1 / (float)degree[1];	 */
	      
/* 	      ay0 = (1.0-fr[0])*src[i0][j0]   + fr[0]*src[i0+1][j0]; */
/* 	      ay1 = (1.0-fr[0])*src[i0][j0+1] + fr[0]*src[i0+1][j0+1]; */
	      
/* 	      dst[i0*degree[0] + i1][j0*degree[1] +j1] = (1.0-fr[1])*ay0 + fr[1]*ay1; */
/* 	    } */
/*       } */

