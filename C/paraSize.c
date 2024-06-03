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
 * Source: $XCRYSDEN_TOPDIR/C/paraSize.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <math.h>

/* --- auxulproc.c ---*/
extern void VecSum2f(float vec1[], float vec2[], int i, int j, float res[]);
extern float distfv(float *v);

double DetermineParapipedSized(double *vec1, double *vec2, 
			       double *vec3, double *orig);

float DetermineParapipedSize(float *vec1, float *vec2, float *vec3, 
			     float *orig);
static void GetParaVec( float *ed, float *vec0, float *vec1, float *vec2,
			 float *res );
static void GetPara2DVec( float *ed, float *vec0, float *vec1, float *res );
double papipedMaxDiagSize(double v0[3], double v1[3], double v2[3]);

double 
DetermineParapipedSized(double *vec1, double *vec2, double *vec3, double *orig)
{
  int i;
  float vecf1[3], vecf2[3], vecf3[3], origf[3];
  
  for (i=0; i<3; i++) {
    vecf1[i] = (float) vec1[i];
    vecf2[i] = (float) vec2[i];
    vecf3[i] = (float) vec3[i];
    origf[i] = (float) orig[i];
  }
  return (double) DetermineParapipedSize(vecf1, vecf2, vecf3, origf);
}


/* get the longest distance of eight parapiped's corners from a
   given origin orig[3] */
float 
DetermineParapipedSize(float *vec1, float *vec2, float *vec3, float *orig)
{
  register int i;
  float max, d;
  float v[3], v1[3];
  float ed[8][3] = {
    {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {1.0, 0.0, 0.0},
    {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}
  };
  
  max = 0.0;
  for( i=0; i<8; i++) {
    GetParaVec( ed[i], vec1, vec2, vec3, v1 );
    VecSum2f( v1, orig, 1, -1, v );
    d = distfv( v );
    if ( max < d ) max = d;
  }
  return max;
}

static void
GetParaVec( float *ed, float *vec0, float *vec1, float *vec2, float *res )
{
  res[0] = ed[0]*vec0[0] + ed[1]*vec1[0] + ed[2]*vec2[0]; 
  res[1] = ed[0]*vec0[1] + ed[1]*vec1[1] + ed[2]*vec2[1]; 
  res[2] = ed[0]*vec0[2] + ed[1]*vec1[2] + ed[2]*vec2[2]; 
}


static void
GetPara2DVec( float *ed, float *vec0, float *vec1, float *res )
{
  res[0] = ed[0]*vec0[0] + ed[1]*vec1[0];
  res[1] = ed[0]*vec0[1] + ed[1]*vec1[1];
  res[2] = ed[0]*vec0[2] + ed[1]*vec1[2];
}


/* get the longest distance of four parallelogram's corners from a
   given origin orig[3] */
float 
DetermineParalleSize(float *vec1, float *vec2, float *orig)
{
  register int i;
  float max, d;
  float v[3], v1[3];
  float ed[4][2] = {
    {0.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {1.0, 0.0} 
  };
  
  max = 0.0;
  for( i=0; i<4; i++) {
    GetPara2DVec( ed[i], vec1, vec2, v1 );
    VecSum2f( v1, orig, 1, -1, v );
    d = distfv( v );
    if ( max < d ) max = d;
  }
  return max;
}


double papipedMaxDiagSize(double v0[3], double v1[3], double v2[3])
{
  register int i;
  double diag_d[4][3], sq_diag[4], max;

  for (i=0; i<3; i++)
    {
      /* 1st diagonal: (0,0,0)-->(1,1,1) */
      diag_d[0][i] = v0[i] + v1[i] + v2[i];
         
      /* 2nd diagonal: (1,1,0)-->(0,0,1) */
      diag_d[1][i] = v2[i] - v0[i] - v1[i]; 
         
      /* 3rd diagonal: (1,0,0)-->(0,1,1) */
      diag_d[2][i] = v1[i] + v2[i] - v0[i];
         
      /* 4th diagonal: (0,1,0)-->(1,0,1) */
      diag_d[3][i] = v0[i] + v2[i] - v1[i];
    }

  for (i=0; i<4; i++) {
    sq_diag[i] =
      diag_d[i][0]*diag_d[i][0] +
      diag_d[i][1]*diag_d[i][1] +
      diag_d[i][2]*diag_d[i][2];
  }

  max = sq_diag[0];    

  if (sq_diag[1] > max) max = sq_diag[1];
  if (sq_diag[2] > max) max = sq_diag[2];
  if (sq_diag[3] > max) max = sq_diag[3];

  return sqrt(max);
}
