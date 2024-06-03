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
 * Source: $XCRYSDEN_TOPDIR/C/splineInt.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "memory.h"

void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void xcRegularSplineInt2(float gridX[], float gridY[], float **Fn,
			 int nX, int nY, float ndegree, float **newFn);
void xcRegSplineInt(float xa[], float ya[], float y2a[], 
		    int i, float x, float *y);
void xcRegularSplineInt3(float gridX[], float gridY[], float gridZ[], 
			 float ***Fn, int nX, int nY, int nZ,
			 float ndegree, float ***newFn);
void spline3(float x[], int n, float ***y, int j, int k, 
	     float yp1, float ypn, float y2[]);
void xcRegSpline3Int(float xa[], float ***ya, int i, int  j, int k, \
		     float y2a[], float x, float *y);


void spline(float x[], float y[], int n, float yp0, float ypnn, float y2[])
{
  /*  
      MEMORIZE: nn = n-1;
  */
  int i,k;
  int nn = n-1;
  float p,qn,sig,un,*u;
  
  u = xcMallocVectorf( nn );
  if (yp0 > 0.99e30) /* The lower boundary condition is set either to be "natural" */
    y2[0] = u[0] =0.0;
  else { /* or else to have a specified first derivative. */
    y2[0] = -0.5;
    u[0]  = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp0);
  }
  for (i=1; i<nn; i++) { 
    /* 
       This is the decomposition loop of the tridiagonal algorithm. 
       y2 and u are used for temporary storage of the decomposed factors.
    */
    sig   = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
    p     = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0) / p;
    u[i]  = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
    u[i]  = (6.0*u[i] / (x[i+1]-x[i-1]) - sig*u[i-1]) / p;
  }
  if (ypnn > 0.99e30) 
    /* The upper boundary condition is set either to be "natural" */
    qn = un =0.0;
  else { /* or else to have a specified first derivative. */
    qn = 0.5;
    un = (3.0 / (x[nn]-x[nn-1])) * (ypnn-(y[nn]-y[nn-1]) / (x[nn]-x[nn-1]));
  }
  y2[nn] = (un-qn*u[nn-1]) / (qn*y2[nn-1]+1.0);
  for (k=nn-1;k>=0;k--)
    /* This is the backsubstitution loop of the tridiagonal algorithm */
    y2[k] = y2[k]*y2[k+1] + u[k];
  xcFree_Vectorf( u );
}


void 
splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a)
{
  /*
    MEMORIZE: nn,mm=n-1,m-1
    -----------------------
  */
  int j;
  
  for (j=0;j<m;j++)
    /* 
       printf("splie2: j=%d\n",j);
       fflush( stdout );
    */
    /* Values 1e30 signal a natural spline */
    spline( x2a, ya[j], n, 1.0e30, 1.0e30, y2a[j] );
}


void 
xcRegularSplineInt2(float gridX[], float gridY[], float **Fn,
		    int nX, int nY, float ndegree, float **newFn)
{
  /* gridX[] ........ x-mesh (can be purely arbitrary)
     gridY[] ........ y-mesh
     **Fn    ........ tabulated function values f(x,y) --> Fn[x][y]
     nX      ........ X-dimension
     nY      ........ Y-dimension
     ndegree ........ degree of interpolation (how many new segments is made 
                      out of old segment; 
		      if oldsegmant = [0,1] && ndegree = 2 --> 
		      seg1 = [0,0.5], seg2 = [0.5,1]
     **newFn ........ new Values goes there
   */
  int ii, jj, nd, i, j, k, l;
  float dx, dy, x, y;
  float **Der2 = xcMallocMatrixf ( nX, nY );
  float *value = xcMallocVectorf ( nX );
  float *der2  = xcMallocVectorf ( nX ); 
  int nXX = nX - 1;
  int nYY = nY - 1;

  splie2( gridX, gridY, Fn, nX, nY, Der2 );

  jj = 0;
  nd = ndegree;
  dy = (gridY[1] - gridY[0]) / (float) nd;
  dx = (gridX[1] - gridX[0]) / (float) nd;
  for (j=0; j<nY; j++) {
    if (j == nYY) nd = 1; /* this is for the last time; */
    for (l=0; l<nd; l++) {
      if (l == 0) {
	for(i=0; i<nX; i++)
	  value[i] = Fn[i][j];
      } else {
	y = gridY[j] + l*dy;
	for(i=0; i<nX; i++)
	  xcRegSplineInt(gridY, Fn[i], Der2[i], j, y, &value[i]);
      }
      
      /* now assing the last value */
      spline(gridX, value, nX, 1.e30, 1.e30, der2);
      
      ii = 0;
      for(i=0; i<nXX; i++) {
	newFn[ii++][jj] = value[i];
	x = gridX[i];
	for (k=1; k<ndegree; k++) {
	  x += dx;
	  xcRegSplineInt(gridX, value, der2, i, x, &newFn[ii++][jj] );
	}
      }
      newFn[ii][jj++] = value[i];
    }
  }

  xcFree_Vectorf ( value );
  xcFree_Vectorf ( der2 );
  xcFree_Matrixf ( Der2 );
}


void 
xcRegSplineInt(float xa[], float ya[], float y2a[], int i, float x, float *y)
{
  int klo,khi;
  float h,b,a;

  klo = i;
  khi = i + 1;
  h = xa[khi]-xa[klo];
  /* The xa's must be distinct. */
  if (h == 0.0) xcError("Bad xa input to routine xcRegSplineInt"); 
  a  = (xa[khi]-x) / h;
  b  = (x-xa[klo]) / h; /* Cubic spline polynomial is now evaluated */
  *y = a*ya[klo] + b*ya[khi] + ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h)/6.0;
}


/* 3D --- 3D --- 3D --- 3D */
void 
xcRegularSplineInt3(float gridX[], float gridY[], float gridZ[], 
		    float ***Fn, int nX, int nY, int nZ,
		    float ndegree, float ***newFn)
{
  /* gridX[] ........ x-mesh (can be purely arbitrary)
     gridY[] ........ y-mesh
     **Fn    ........ tabulated function values f(x,y) --> Fn[x][y]
     nX      ........ X-dimension
     nY      ........ Y-dimension
     ndegree ........ degree of interpolation (how many new segments is made 
                      out of old segment; 
		      if oldsegmant = [0,1] && ndegree = 2 --> 
		      seg1 = [0,0.5], seg2 = [0.5,1]
     **newFn ........ new Values goes there
   */
  int ii, jj, kk, nd1, nd2, i, j, k, l, ij, ik;
  int nny = ndegree * (nY-1) + 1;
  int nnz = ndegree * (nZ-1) + 1;
  int nXX = nX - 1;
  int nYY = nY - 1;
  int nZZ = nZ - 1;
  float dx, x;  
  float *der2   = xcMallocVectorf ( nX ); 
  float ***newf = xcMallocTensor3f ( nX, nny, nnz );
  
  for(i=0; i<nX; i++) {
    printf("i=%d\n",i);
    fflush(stdout);
    xcRegularSplineInt2(gridY, gridZ, Fn[i], nY, nZ, ndegree, newf[i]);
  }

  nd1 = ndegree;
  kk=0;
  dx = (gridX[1] - gridX[0]) / (float) ndegree;
  for(k=0; k<nZ; k++) {
    if (k==nZZ) nd1 = 1;
    for(ik=0; ik<nd1; ik++) {
      jj = 0;
      nd2 = ndegree;
      for(j=0; j<nY; j++) {
	if (j==nYY) nd2 = 1;
	for(ij=0; ij<nd2; ij++) {
	  spline3(gridX, nX, newf, jj, kk, 1.e30, 1.e30, der2);
	  ii = 0;
	  for(i=0; i<nXX; i++) {	   
	    newFn[ii++][jj][kk] = newf[i][jj][kk];
	    x = gridX[i];
	    for(l=1; l<ndegree; l++) {
	      x += dx;
	      xcRegSpline3Int( gridX, newf, i, jj, kk, der2, x, 
			       &newFn[ii++][jj][kk] );
	    }
	  }
	  newFn[ii][jj][kk] = newf[i][jj][kk];
	  jj++;
	} /* ij */
      } /* j */
      kk++;
    } /* ik */
  } /* k */

  xcFree_Vectorf( der2 );
  xcFree_Tensor3f( newf );
}


void spline3(float x[], int n, float ***y, int j, int k, 
	     float yp0, float ypnn, float y2[])
{
  /*  
      MEMORIZE: nn=n-1
      ----------------
  */
  int i, l;
  int nn = n - 1;
  float p, qn, sig, un, *u;
  u = xcMallocVectorf( nn );
  if (yp0 > 0.99e30) /* The lower boundary condition is set either to be "natural" */
    y2[0] = u[0] = 0.0;
  else { /* or else to have a specified first derivative. */
    y2[0] = -0.5;
    u[0]  = (3.0 / (x[1]-x[0])) * ((y[1][j][k]-y[0][j][k]) / (x[1]-x[0])-yp0);
  }
  for (i=1; i<nn; i++) { 
    /* 
       This is the decomposition loop of the tridiagonal algorithm. 
       y2 and u are used for temporary storage of the decomposed factors.
    */
    sig   = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
    p     = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0) / p;
    u[i]  = (y[i+1][j][k]-y[i][j][k]) / (x[i+1]-x[i]) - (y[i][j][k]-y[i-1][j][k]) / (x[i]-x[i-1]);
    u[i]  = (6.0*u[i] / (x[i+1]-x[i-1]) - sig*u[i-1]) / p;
  }
  if (ypnn > 0.99e30) 
    /* The upper boundary condition is set either to be "natural" */
    qn = un =0.0;
  else { /* or else to have a specified first derivative. */
    qn = 0.5;
    un = (3.0 / (x[nn]-x[nn-1])) * (ypnn-(y[nn][j][k]-y[nn-1][j][k]) / (x[nn]-x[nn-1]));
  }
  y2[nn] = (un-qn*u[nn-1]) / (qn*y2[nn-1]+1.0);
  for (l=nn-1; l>=0; l--) 
    /* This is the backsubstitution loop of the tridiagonal algorithm */
    y2[l] = y2[l]*y2[l+1] + u[l];
  xcFree_Vectorf( u );
}



void 
xcRegSpline3Int(float xa[], float ***ya, int i, int  j, int k, \
		float y2a[], float x, float *y)
{
  int klo,khi;
  float h,b,a;

  klo = i;
  khi = i + 1;
  h = xa[khi]-xa[klo];
  /* The xa's must be distinct. */
  if (h == 0.0) xcError("Bad xa input to routine splint"); 
  a  = (xa[khi]-x) / h;
  b  = (x-xa[klo]) / h; /* Cubic spline polynomial is now evaluated */
  *y = a*ya[klo][j][k] + b*ya[khi][j][k] + 
    ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h)/6.0;
}
