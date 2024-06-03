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
 * Source: $XCRYSDEN_TOPDIR/C/lcasi.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/



/* ------------------------------------------------------------------------
 *  
 * unsigned-char version of LCASI interpolation
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


/* lcasi.c */

extern void cryError(int status, char *function, char *message);
extern unsigned char **cryRegPeriodInterpolator_uc_LCASI2D(int n[2], int degree[2], unsigned char **src);
extern void cryRegPeriodInterpolate_uc_LCASI1D(int N, int DEGREE, unsigned char *src, unsigned char *dst);
extern void cryRegPeriodInterpolate_uc_LCASI2D(int N1, int N2, int DEGREE1, int DEGREE2, unsigned char **src, unsigned char **dst);


/* V = MATRIX .mult. VECTOR */
#define M33_MULT_V3(dst, m, v ) \
do { \
  (dst)[0] = (m)[0][0] * (v)[0]  +  (m)[0][1] * (v)[1]  +  (m)[0][2] * (v)[2];\
  (dst)[1] = (m)[1][0] * (v)[0]  +  (m)[1][1] * (v)[1]  +  (m)[1][2] * (v)[2];\
  (dst)[2] = (m)[2][0] * (v)[0]  +  (m)[2][1] * (v)[1]  +  (m)[2][2] * (v)[2];\
} while(0)

void cryError(int status, char *function, char *message) {
  fprintf(stderr,"\n ==============================================================================\n");
  fprintf(stderr, "  ERROR: an error has occurred in function %s\n", function);
  fprintf(stderr, "         %s\n", message);
  fprintf(stderr," ==============================================================================\n\n");
  exit(status);
}



/*!
 *-----------------------------------------------------------------------------
 *
 * cryRegPeriodInterpolator_uc_LCASI2D -- \brief 2D LCASI interpolator
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

unsigned char **cryRegPeriodInterpolator_uc_LCASI2D(int n[2], /*!< number of input points*/
						    int degree[2], /*!< degree of interp */
						    unsigned char **src) /*!< INPUT func. values*/
{
  int n_out[2];
  unsigned char **dst;

  n_out[0] = n[0] * degree[0];
  n_out[1] = n[1] * degree[1];
  
  MALLOC_TENSOR2(unsigned char, dst, n_out[0], n_out[1]);  
  cryRegPeriodInterpolate_uc_LCASI2D(n[0], n[1], 
				    degree[0], degree[1], src, dst);
  return dst;
}


/*!
 *-----------------------------------------------------------------------------
 *
 * cryRegPeriodInterpolate_uc_LCASI1D -- \brief 1D LCASI interpolation routine
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

void cryRegPeriodInterpolate_uc_LCASI1D(int N, /* number of input points*/
					int DEGREE, /*degree of interp   */
					unsigned char *src,/*INPUT func. values */
					unsigned char *dst)/*OUTPUT func. values*/
{
  /* *dst must be preallocated !!! */

  int    i, j, ii, point;
  float x, yL[3], yR[3], abcL[3], abcR[3], wL, wR;

  if ( dst == NULL ) 
    cryError(1, "cryRegPeriodInterpolate_uc_LCASI1D",
	     "dst pointer is NULL");
  
  /* calculate 1st Left (a,b,c) coefficients */
  yL[0] = (float) src[N-1];
  yL[1] = (float) src[0];
  yL[2] = (float) src[1];
  M33_MULT_V3 (abcL, invL, yL);

  for (i=0; i<N; i++) {
    ii = i * DEGREE;

    /* calculate Right (a,b,c) coefficients */
    for (j=0; j<3; j++) {
      if (i+j < N) yR[j] = (float) src[i+j];
      else yR[j] = (float) src[i+j - N];
    }
    M33_MULT_V3(abcR, invR, yR);
    
    dst[ii] = src[i];
    for (j=1; j<DEGREE; j++) {
      x = (float) j / (float) DEGREE;
      wR = x;
      wL = 1.0 - x;
      point = (int)
	wL * (x*(abcL[0]*x + abcL[1]) + abcL[2]) +
	wR * (x*(abcR[0]*x + abcR[1]) + abcR[2]);

      if (point > 255) dst[ii+j] = 255;
      else if (point < 0 ) dst[ii+j] = 0;
      else dst[ii+j] = (unsigned char) point;
    }

    /* calculate Left (a,b,c) coefficients */
    for (j=0; j<3; j++) {
      if (i+j < N) yL[j] = (float) src[i+j];
      else yL[j] = (float) src[i+j - N];
    }
    M33_MULT_V3(abcL, invL, yL);    
  }
}


/*!
 *-----------------------------------------------------------------------------
 *
 * cryRegPeriodInterpolate_uc_LCASI2D -- \brief 2D LCASI interpolation routine
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

void cryRegPeriodInterpolate_uc_LCASI2D(int N1, 
					int N2,
					int DEGREE1,
					int DEGREE2,
					unsigned char **src,
					unsigned char **dst) {
  int ID1 = DEGREE1;
  int ID2 = DEGREE2;
  int IN1 = N1*ID1;
  int IN2 = N2*ID2;
  int i1, ii1, ii2, j, point;
  float x, yL[3], yR[3], abcL[3], abcR[3], wL, wR;

  /* **dst must be preallocated !!! */

  if ( dst == NULL ) 
    cryError(1, "cryRegPeriodInterpolate_uc_LCASI2D",
	     "dst pointer is NULL");


  /* 
   * doing the 1st dimension (means: intepolating the 2nd)
   */
  for (i1=0; i1<N1; i1++)
    cryRegPeriodInterpolate_uc_LCASI1D (N2,ID2,src[i1],dst[i1*ID1]);  
    
  /*
   * doing the 2nd dimension (means: intepolating the 1st)
   */
  for (ii2=0; ii2<IN2; ii2++) {
    /* calculate 1st Left (a,b,c) coefficients */
    yL[0] = (float) dst[IN1-ID1][ii2];
    yL[1] = (float) dst[0      ][ii2];
    yL[2] = (float) dst[ID1    ][ii2];
    M33_MULT_V3 (abcL, invL, yL);
    
    /* intepolating the 1st dimension */
    for (i1=0; i1<N1; i1++) {
      ii1 = i1*ID1;

      /* calculate Right (a,b,c) coefficients */
      for (j=0; j<3; j++) {
	if (i1+j < N1) yR[j] = (float) dst[ii1+j*ID1][ii2];
	else yR[j] = (float) dst[ii1+j*ID1 - IN1][ii2];
      }
      M33_MULT_V3(abcR, invR, yR);
      
      for (j=1; j<ID1; j++) {
	x = (float) j / (float) ID1;
	wR = x;
	wL = 1.0 - x;
	
	point = (int)
	  wL * (x*(abcL[0]*x + abcL[1]) + abcL[2]) +
	  wR * (x*(abcR[0]*x + abcR[1]) + abcR[2]);

	if (point > 255) dst[ii1+j][ii2] = 255;
	else if (point < 0 ) dst[ii1+j][ii2] = 0;
	else dst[ii1+j][ii2] = (unsigned char) point;
      }
      
      /* calculate Left (a,b,c) coefficients */
      for (j=0; j<3; j++) {
	if (i1+j < N1) yL[j] = (float) dst[ii1+j*ID1][ii2];
	else yL[j] = (float) dst[ii1+j*ID1 - IN1][ii2];
      }
      M33_MULT_V3(abcL, invL, yL);    
    }
  }
}
