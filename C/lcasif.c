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
 * Source: $XCRYSDEN_TOPDIR/C/lcasif.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

/* ------------------------------------------------------------------------
 *  
 * float version of LCASI interpolation
 *
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include "tensor.h"
#include "memory.h"

static double invL[3][3] = {
  { 0.5,-1.0, 0.5},
  {-0.5, 0.0, 0.5},
  { 0.0, 1.0, 0.0}
};
static double invR[3][3] = {
  { 0.5,-1.0, 0.5},
  {-1.5, 2.0,-0.5},
  { 1.0, 0.0, 0.0}
};

#define M33_MULT_V3(dst, m, v ) \
do { \
  (dst)[0] = (m)[0][0] * (v)[0]  +  (m)[0][1] * (v)[1]  +  (m)[0][2] * (v)[2];\
  (dst)[1] = (m)[1][0] * (v)[0]  +  (m)[1][1] * (v)[1]  +  (m)[1][2] * (v)[2];\
  (dst)[2] = (m)[2][0] * (v)[0]  +  (m)[2][1] * (v)[1]  +  (m)[2][2] * (v)[2];\
} while(0)


/* lcasi.c */
extern void cryError(int status, char *function, char *message);

/* lcasif.c */
extern float *cryRegPeriodInterpolator_f_LCASI1D(int N, int DEGREE, float *src);
extern float **cryRegPeriodInterpolator_f_LCASI2D(int n[2], int degree[2], float **src);
extern float ***cryRegPeriodInterpolator_f_LCASI3D(int n[3], int degree[3], float ***src);
extern void cryRegPeriodInterpolate_f_LCASI1D(int N, int DEGREE, float *src, float *dst);
extern void cryRegPeriodInterpolate_f_LCASI2D(int N1, int N2, int DEGREE1, int DEGREE2, float **src, float **dst);
extern void cryRegPeriodInterpolate_f_LCASI3D(int N1, int N2, int N3, int DEGREE1, int DEGREE2, int DEGREE3, float ***src, float ***dst);



/*!
 *-----------------------------------------------------------------------------
 *
 * cryRegPeriodInterpolator_f_LCASI1D -- \brief 1D LCASI interpolator
 *
 * \description

 *    This function is an object driver for 2D LCASI interpolation on
 *    regular periodic grid. The output grid can only be of
 *    integer-multiple size compared to input grid. The integer degree
 *    interpolation is specified via "DEGREE" input parameter.
 *
 * \return
 *    Interpolated 1D regular periodic datagrid.
 *
 * \sideeffect
 *    None.
 *
 *-----------------------------------------------------------------------------
 */

float *cryRegPeriodInterpolator_f_LCASI1D(int N, /* number of input points*/
					   int DEGREE, /*degree of interp   */
					   float *src)/*INPUT func. values */
{
  int n_out = N * DEGREE;
  float *dst;

  dst = (float*) malloc(sizeof(float) * n_out);
  if (!dst) {
    fprintf( stderr, 
	     "\n***\nMemory Error: can't allocate %ld bytes\n", 
	     sizeof(float) * n_out );
    abort();
  }

  cryRegPeriodInterpolate_f_LCASI1D(N, DEGREE, src, dst);
  return dst;
}


/*!
 *-----------------------------------------------------------------------------
 *
 * cryRegPeriodInterpolator_f_LCASI2D -- \brief 2D LCASI interpolator
 *
 * \description

 *    This function is an object driver for 2D LCASI interpolation on
 *    regular periodic grid. The output grid can only be of
 *    integer-multiple size compared to input grid. The integer degree
 *    interpolation is specified via "degree" input parameter.
 *
 * \return
 *    Interpolated 2D regular periodic datagrid.
 *
 * \sideeffect
 *    None.
 *
 *-----------------------------------------------------------------------------
 */

float **cryRegPeriodInterpolator_f_LCASI2D(int n[2], /*!< number of input points*/
					    int degree[2], /*!< degree of interp */
					    float **src) /*!< INPUT func. values*/
{
  int n_out[2];
  float **dst;

  n_out[0] =  n[0] * degree[0];
  n_out[1] =  n[1] * degree[1];
  
  MALLOC_TENSOR2(float, dst, n_out[0], n_out[1]);  
  cryRegPeriodInterpolate_f_LCASI2D(n[0], n[1], 
				    degree[0], degree[1], src, dst);
  return dst;
}



/*!
 *-----------------------------------------------------------------------------
 *
 * cryRegPeriodInterpolator_f_LCASI3D -- \brief 3D LCASI interpolator
 *
 * \description

 *    This function is an object driver for 3D LCASI interpolation on
 *    regular periodic grid. The output grid can only be of
 *    integer-multiple size compared to input grid. The integer degree
 *    interpolation is specified via "degree" input parameter.
 *
 * \return
 *    Interpolated 3D regular periodic datagrid.
 *
 * \sideeffect
 *    None.
 *
 *-----------------------------------------------------------------------------
 */

float ***cryRegPeriodInterpolator_f_LCASI3D(int n[3],/*!<number of input points*/
					     int degree[3],/*!< degree of interp   */
					     float ***src)/*!< INPUT func. values */
{
  int n_out[3];
  float ***dst;
  
  n_out[0] =  n[0] * degree[0];
  n_out[1] =  n[1] * degree[1];
  n_out[2] =  n[2] * degree[2];

  MALLOC_TENSOR3(float, dst, n_out[0], n_out[1], n_out[2]);  
  cryRegPeriodInterpolate_f_LCASI3D(n[0], n[1], n[2],
				    degree[0], degree[1], degree[2], src, dst);
  return dst;
}


/*
  a version of periodic cryRegPeriodInterpolator_f_LCASI3D that operates on
  general grids !!!
*/
float ***cryGeneralGridRegPeriodInterpolator_f_LCASI3D(int n[3],      /*!< number of input points of general grid */
						       int degree[3], /*!< degree of interp   */
						       float ***src)  /*!< INPUT general grid */
{
  register int i0, i1, i2, ii0, ii1, ii2, n0, n1, n2;
  float ***aux_src, ***aux, ***dst;
  
  /* 
     copy general "src" to periodic "src_aux" grid 
  */
  MALLOC_TENSOR3(float, aux_src, n[0]-1, n[1]-1, n[2]-1);
  for (i0=0; i0<n[0]-1; i0++) 
    for (i1=0; i1<n[1]-1; i1++) 
      for (i2=0; i2<n[2]-1; i2++) 
	aux_src[i0][i1][i2] = src[i0][i1][i2];
  
  /*
    interpolate
  */
  n0 =  (n[0]-1) * degree[0];
  n1 =  (n[1]-1) * degree[1];
  n2 =  (n[2]-1) * degree[2];
  MALLOC_TENSOR3(float, aux, n0, n1, n2);  
  MALLOC_TENSOR3(float, dst, n0+1, n1+1, n2+1);  
  cryRegPeriodInterpolate_f_LCASI3D(n[0]-1, n[1]-1, n[2]-1, 
				    degree[0], degree[1], degree[2], 
				    src, aux);

  /* 
     copy interpolated periodic grid to a general grid 
  */
  for (i0=0; i0<n0+1; i0++) 
    {
      ii0 = i0 % n0;
      for (i1=0; i1<n1+1; i1++) 
	{
	  ii1 = i1 % n1;
	  for (i2=0; i2<n2+1; i2++) 
	    {
	      ii2 = i2 % n2;
	      dst[i0][i1][i2] = aux[ii0][ii1][ii2];
	    }
	}
    }
  xcFree_Tensor3f (aux_src);
  xcFree_Tensor3f (aux);
  return dst;
}



/*!
 *-----------------------------------------------------------------------------
 *
 * cryRegPeriodInterpolate_f_LCASI1D -- \brief 1D LCASI interpolation routine
 *
 * \description

 *    This function performs 2D LCASI interpolation on regular
 *    periodic grid. The output grid can only be of integer-multiple
 *    size compared to input grid. The integer degree interpolation is
 *    specified via "degree" input parameter.
 *
 * \return
 *    Nothing.
 *
 * \sideeffect
 *    On exit "dst" contains interpolated 1D regular periodic datagrid.
 *
 *-----------------------------------------------------------------------------
 */

void cryRegPeriodInterpolate_f_LCASI1D(int N, /* number of input points*/
				       int DEGREE, /*degree of interp   */
				       float *src,/*INPUT func. values */
				       float *dst)/*OUTPUT func. values*/
{
  /* *dst must be preallocated !!! */

  int    i, j, ii;
  float x, yL[3], yR[3], abcL[3], abcR[3], wL, wR;

  if ( dst == NULL ) 
    cryError(1, "cryRegPeriodInterpolate_f_LCASI1D",
	     "dst pointer is NULL");
  
  /* calculate 1st Left (a,b,c) coefficients */
  yL[0] = src[N-1];
  yL[1] = src[0];
  yL[2] = src[1];
  M33_MULT_V3 (abcL, invL, yL);

  for (i=0; i<N; i++) {
    ii = i * DEGREE;

    /* calculate Right (a,b,c) coefficients */
    for (j=0; j<3; j++) {
      if (i+j < N) yR[j] = src[i+j];
      else yR[j] = src[i+j - N];
    }
    M33_MULT_V3(abcR, invR, yR);
    
    dst[ii] = src[i];
    for (j=1; j<DEGREE; j++) {
      x = (float) j / (float) DEGREE;
      wR = x;
      wL = 1.0 - x;
      dst[ii+j] =
	wL * (x*(abcL[0]*x + abcL[1]) + abcL[2]) +
	wR * (x*(abcR[0]*x + abcR[1]) + abcR[2]);
    }

    /* calculate Left (a,b,c) coefficients */
    for (j=0; j<3; j++) {
      if (i+j < N) yL[j] = src[i+j];
      else yL[j] = src[i+j - N];
    }
    M33_MULT_V3(abcL, invL, yL);    
  }
}


/*!
 *-----------------------------------------------------------------------------
 *
 * cryRegPeriodInterpolate_f_LCASI2D -- \brief 2D LCASI interpolation routine
 *
 * \description

 *    This function performs 2D LCASI interpolation on regular
 *    periodic grid. The output grid can only be of integer-multiple
 *    size compared to input grid. The integer degree interpolation is
 *    specified via "degree" input parameter.
 *
 * \return
 *    Nothing.
 *
 * \sideeffect
 *    On exit "dst" contains interpolated 2D regular periodic datagrid.
 *
 *-----------------------------------------------------------------------------
 */

void cryRegPeriodInterpolate_f_LCASI2D(int N1, 
				       int N2,
				       int DEGREE1,
				       int DEGREE2,
				       float **src,
				       float **dst) {
  int ID1 = DEGREE1;
  int ID2 = DEGREE2;
  int IN1 = N1*ID1;
  int IN2 = N2*ID2;
  int i1, ii1, ii2, j;
  float x, yL[3], yR[3], abcL[3], abcR[3], wL, wR;

  /* **dst must be preallocated !!! */

  if ( dst == NULL ) 
    cryError(1, "cryRegPeriodInterpolate_f_LCASI2D",
	     "dst pointer is NULL");


  /* 
   * doing the 1st dimension (means: intepolating the 2nd)
   */
  for (i1=0; i1<N1; i1++)
    cryRegPeriodInterpolate_f_LCASI1D (N2,ID2,src[i1],dst[i1*ID1]);  
    
  /*
   * doing the 2nd dimension (means: intepolating the 1st)
   */
  for (ii2=0; ii2<IN2; ii2++) {
    /* calculate 1st Left (a,b,c) coefficients */
    yL[0] = dst[IN1-ID1][ii2];
    yL[1] = dst[0      ][ii2];
    yL[2] = dst[ID1    ][ii2];
    M33_MULT_V3 (abcL, invL, yL);
    
    /* intepolating the 1st dimension */
    for (i1=0; i1<N1; i1++) {
      ii1 = i1*ID1;

      /* calculate Right (a,b,c) coefficients */
      for (j=0; j<3; j++) {
	if (i1+j < N1) yR[j] = dst[ii1+j*ID1][ii2];
	else yR[j] = dst[ii1+j*ID1 - IN1][ii2];
      }
      M33_MULT_V3(abcR, invR, yR);
    
      for (j=1; j<ID1; j++) {
	x = (float) j / (float) ID1;
	wR = x;
	wL = 1.0 - x;
	dst[ii1+j][ii2] =
	  wL * (x*(abcL[0]*x + abcL[1]) + abcL[2]) +
	  wR * (x*(abcR[0]*x + abcR[1]) + abcR[2]);
      }

      /* calculate Left (a,b,c) coefficients */
      for (j=0; j<3; j++) {
	if (i1+j < N1) yL[j] = dst[ii1+j*ID1][ii2];
	else yL[j] = dst[ii1+j*ID1 - IN1][ii2];
      }
      M33_MULT_V3(abcL, invL, yL);    
    }
  }
}



/*!
 *-----------------------------------------------------------------------------
 *
 * cryRegPeriodInterpolate_f_LCASI3D -- \brief 3D LCASI interpolation routine
 *
 * \description

 *    This function performs 3D LCASI interpolation on regular
 *    periodic grid. The output grid can only be of integer-multiple
 *    size compared to input grid. The integer degree interpolation is
 *    specified via "degree" input parameter.
 *
 * \return
 *    Nothing.
 *
 * \sideeffect
 *    On exit "dst" contains interpolated 3D regular periodic datagrid.
 *
 *-----------------------------------------------------------------------------
 */

void cryRegPeriodInterpolate_f_LCASI3D(int N1, 
				       int N2,
				       int N3,
				       int DEGREE1,
				       int DEGREE2,
				       int DEGREE3,
				       float ***src,
				       float ***dst) {
  register int ID1 = DEGREE1;
  register int ID2 = DEGREE2;
  register int ID3 = DEGREE3;
  register int IN1 = N1*ID1;
  register int IN2 = N2*ID2;
  register int IN3 = N3*ID3;
  register int i1, ii1, ii2, ii3, j;
  float x, yL[3], yR[3], abcL[3], abcR[3], wL, wR;

  /* **dst must be preallocated !!! */

  if ( dst == NULL ) 
    cryError(1, "cryRegPeriodInterpolate_f_LCASI3D",
	     "dst pointer is NULL");

  /* 
   * doing the 1st dimension (means: intepolating the 2nd & 3rd)
   */
  for (i1=0; i1<N1; i1++)
    cryRegPeriodInterpolate_f_LCASI2D (N2, N3, ID2, ID3,
				       src[i1], dst[i1*ID1]);  
    
  /*
   * doing the 2nd & 3th dimension (means: intepolating the 1st)
   */
  for (ii2=0; ii2<IN2; ii2++) {
    for (ii3=0; ii3<IN3; ii3++) {
      /* calculate 1st Left (a,b,c) coefficients */
      yL[0] = dst[IN1-ID1][ii2][ii3];
      yL[1] = dst[0      ][ii2][ii3];
      yL[2] = dst[ID1    ][ii2][ii3];
      M33_MULT_V3 (abcL, invL, yL);

      /* intepolating the 1st dimension */
      for (i1=0; i1<N1; i1++) {
	ii1 = i1*ID1;

	/* calculate Right (a,b,c) coefficients */
	for (j=0; j<3; j++) {
	  if (i1+j < N1) yR[j] = dst[ii1+j*ID1][ii2][ii3];
	  else yR[j] = dst[ii1+j*ID1 - IN1][ii2][ii3];
	}
	M33_MULT_V3(abcR, invR, yR);
    
	for (j=1; j<ID1; j++) {
	  x = (float) j / (float) ID1;
	  wR = x;
	  wL = 1.0 - x;
	  dst[ii1+j][ii2][ii3] =
	    wL * (x*(abcL[0]*x + abcL[1]) + abcL[2]) +
	    wR * (x*(abcR[0]*x + abcR[1]) + abcR[2]);
	}

	/* calculate Left (a,b,c) coefficients */
	for (j=0; j<3; j++) {
	  if (i1+j < N1) yL[j] = dst[ii1+j*ID1][ii2][ii3];
	  else yL[j] = dst[ii1+j*ID1 - IN1][ii2][ii3];
	}
	M33_MULT_V3(abcL, invL, yL);    
      }
    }
  }
}


