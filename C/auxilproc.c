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
 * Source: $XCRYSDEN_TOPDIR/C/auxilproc.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

/* FUNCTION PROTOTYPES */
double dist6(double x1, double x2, double y1, double y2, double z1, double z2);
double dist3(double x, double y, double z);
double distdv(double *vec);
int normalizepv(double *x, double *y, double *z);
int iround(double x);
void Rotate(double *x, double *y, double *z, double cosfi, double sinfi);
void xcMat44Copyd( double mat1[][4], double mat2[][4], int n, int m );

/* float version of functions */
int normalizepvf(float *x, float *y, float *z);
int normalizepvfv(float *vec);
void RevertVectorfv(float *vec);
float dist3f(float x, float y, float z);
float distfv(float *v);
int iroundf(float x);
void xcDeleteBitFlags(unsigned int *var, unsigned int flags);
void MatVecMult33f(float mat[3][3], float vec[3], float new[3]);
void MatMult33f(float a[3][3], float b[3][3], float new[3][3]);
void MatCopy33f(float a[3][3], float b[3][3]);
void MatToZero33f(float a[3][3]);
int IsSamePointvf(float *a, float *b, float limit);
int IsEqualf(float a, float b, float limit);
int WriteArgv(char *argv[], int *argc, char *fmt, ... );
void VecSum2f(float vec1[], float vec2[], int i, int j, float res[]);
void VecSum3f(float vec1[], float vec2[], float vec3[], 
	      int i, int j, int k, float res[]);
void VecProduct3f(float v0[], float v1[], float res[]);
float CramerRule2x2Detf(float d[2][3], int i);
void xcRotate3fv( float fi, float *v, float *m, float *p );

double det3x3d(const double *v0, const double *v1, const double *v2);
float det3x3f(const float *v0, const float *v1, const float *v2);

/* DISTANCE calculate distance between two points in space */
double 
dist6(double x1, double x2, double y1, double y2, double z1, double z2)
{
  double xij2, yij2, zij2, dis;
  /*printf("In DIST6%s\n","=DIST6");*/
  xij2 = (x2 - x1) * (x2 - x1);
  yij2 = (y2 - y1) * (y2 - y1); 
  zij2 = (z2 - z1) * (z2 - z1);
  dis = sqrt(xij2 + yij2 + zij2);
  
  return(dis);
}


/* dist3 = lenght of vector (x,y,z) */
double 
dist3(double x, double y, double z)
{
  double dis;
  dis = sqrt(x * x + y * y + z * z);
  /*printf("dis=%f\n",dis);*/
  return(dis);
}
    

/* distdv = lenght of 3D vector vec */
double 
distdv(double *vec)
{
  return( sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]) );
}
    

/* normalize a vector; pointers must be submited
   when function is called */
int 
normalizepv(double *x, double *y, double *z)
{
  double dis;

  dis = dist3(*x, *y, *z);
  if ( dis < 1.e-12 ) {
    fprintf(stderr,"WARNING: zero length vector\n");
    return 0;
  }
  *x /= dis;
  *y /= dis;
  *z /= dis;

  return 1;
}


/* custom ROUND FUNCTION */
int
iround(double x)
{
  int n, m, mm;
  
  n = (int) x;
  m = 10 * n;
  mm = (int) ( 10 * x );
  if ( mm - m  < 5 ) return n;
  else return n + 1;
}


/* ROTATION PROC; this is rotation around Z axe */
void
Rotate(double *x, double *y, double *z, double cosfi, double sinfi)
{
  double x1, y1;
  
  x1 =  *x * cosfi + -*y * sinfi;
  y1 =  *x * sinfi + *y * cosfi;
  *x = x1;
  *y = y1;
}



/*****************************************************************************/
/* function copy mat2[][4] to mat1[][4]                                      */
void
xcMat44Copyd( double mat1[][4], double mat2[][4], int n, int m ) 
{
  int i, j;
  
  for ( i=0; i<n; i++)
    for (j=0; j<m; j++)
      mat1[i][j] = mat2[i][j];
}


/*****************************************************************************/
/* here are float versions of functions                                      */
int 
normalizepvf(float *x, float *y, float *z)
{
  float dis;

  dis = dist3f(*x, *y, *z);
  /* printf("DISTANCE:: distance = %f",dis);
     fflush(stdout); */
  if ( dis < 1.e-7 ) {
    fprintf(stderr,"WARNING: zero length vector\n");
    return 0;
  }
  *x /= dis;
  *y /= dis;
  *z /= dis;

  return 1;
}


int 
normalizepvfv(float *vec)
{
  float dis;

  dis = dist3f(vec[0], vec[1], vec[2]);
  /*  printf("DISTANCE:: distance = %f",dis);
  fflush(stdout); */
  if ( dis < 1.e-7 ) {
    fprintf(stderr,"WARNING: zero length vector\n");
    return 0;
  }
  vec[0] /= dis;
  vec[1] /= dis;
  vec[2] /= dis;

  return 1;
}


void
RevertVectorfv(float *vec)
{
  vec[0] *= -1.0;
  vec[1] *= -1.0;
  vec[2] *= -1.0;
}


float 
dist3f(float x, float y, float z)
{
  float dis;
  dis = sqrt(x * x + y * y + z * z);
  return(dis);
}


float 
distfv(float *v)
{
  float dis;
  dis = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  return(dis);
}


/* custom FLOAT-ROUND FUNCTION */
int
iroundf(float x)
{
  int n, m, mm;
  
  n = (int) x;
  m = 10 * n;
  mm = (int) ( 10 * x );
  if ( mm - m  < 5 ) return n;
  else return n + 1;
}


/*****************************************************************************/
/* this function delete flags-bits from *var                                 */
void
xcDeleteBitFlags(unsigned int *var, unsigned int flags)
{
  unsigned int i, res;

  res = *var & flags;
  for (i=1; i<=res; i*=2) 
    if ( res & i ) *var ^= i; 
}


void 
MatVecMult33f(float mat[3][3], float vec[3], float new[3])
{
  new[0] = mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2];
  new[1] = mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2];
  new[2] = mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2];
}


void
MatMult33f(float a[3][3], float b[3][3], float new[3][3])
{
  int i, j, k;
  
  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      for(k=0; k<3; k++)
	new[i][j] = new[i][j] + a[i][k] * b[k][j];
}


void
MatCopy33f(float a[3][3], float b[3][3])
{ 
  int i,j;

  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      b[i][j] = a[i][j];
}


void
MatToZero33f(float a[3][3])
{
  int i, j;

  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      a[i][j] = 0.0;
}

    
int
IsSamePointvf(float *a, float *b, float limit)
{
  if ( IsEqualf(a[0], b[0], limit) && \
       IsEqualf(a[1], b[1], limit) && \
       IsEqualf(a[2], b[2], limit) ) {
    return 1;
  } else {
    return 0;
  }
}


int 
IsEqualf(float a, float b, float limit)
{
  /*double fabs();*/

  if ( fabs(a-b) < limit ) {
    return 1;
  } else {
    return 0;
  }
}


int
WriteArgv(char *argv[], int *argc, char *fmt, ... )
{
  va_list ap;
  int    d;
  double f;
  
  va_start(ap, fmt);
  while (*fmt)
    switch(*fmt++) {
    case 's':           /* string */
      strcpy(argv[(*argc)++], va_arg(ap, char *));
      break;
    case 'd':           /* int */
      d = va_arg(ap, int);
      sprintf(argv[(*argc)++],"%d", d);
      break;
    case 'f':           /* float */
      f = (double) va_arg(ap, double);
      sprintf(argv[(*argc)++],"%f", f);
      break;
    }
  va_end(ap);
  return 1;
}


void VecSum2f(float vec1[], float vec2[], int i, int j, float res[])
{
  res[0] = (float)i*vec1[0] + (float)j*vec2[0];
  res[1] = (float)i*vec1[1] + (float)j*vec2[1];
  res[2] = (float)i*vec1[2] + (float)j*vec2[2];
}

void VecSum3f(float vec1[], float vec2[], float vec3[], 
	      int i, int j, int k, float res[])
{
  res[0] = (float)i*vec1[0] + (float)j*vec2[0] + (float)k*vec3[0];
  res[1] = (float)i*vec1[1] + (float)j*vec2[1] + (float)k*vec3[1];
  res[2] = (float)i*vec1[2] + (float)j*vec2[2] + (float)k*vec3[2];
}

void VecProduct3f(float v0[], float v1[], float res[])
{
  res[0] = v0[1] * v1[2] - v1[1] * v0[2];
  res[1] = v0[2] * v1[0] - v1[2] * v0[0];
  res[2] = v0[0] * v1[1] - v1[0] * v0[1];
}


/* 
   solutions of 2x2 linear equantions problem;
   returns determinants D or Dx or Dy; example: x = Dx/D
*/
float 
CramerRule2x2Detf(float d[2][3], int i)
{
  float det;
  if ( i== 0 ) {
    det = d[0][0]*d[1][1] - d[1][0]*d[0][1];
  } else if ( i == 1 ) {
    det = d[0][2]*d[1][1] - d[1][2]*d[0][1];
  } else {
    det = d[0][0]*d[1][2] - d[1][0]*d[0][2];
  }
  return det;
}


/*
  rotation of point *p for fi around vector *v and point *m
*/
void
xcRotate3fv( float fi, float *v, float *m, float *p )
{
  /* fi in radians */
  register int i, j;
  float u[3], p1[3], cosfi, sinfi, uij, uij1; 
  float S[3][3] = { {0.0, 0.0, 0.0},
		    {0.0, 0.0, 0.0},
		    {0.0, 0.0, 0.0} };
  float M[3][3] = { {0.0, 0.0, 0.0},
		    {0.0, 0.0, 0.0},
		    {0.0, 0.0, 0.0} };
  float I[3][3] = { {1.0, 0.0, 0.0},
		    {0.0, 1.0, 0.0},
		    {0.0, 0.0, 1.0} };

  cosfi = cos( fi );
  sinfi = sin( fi );

  for (i=0; i<3; i++) {
    u[i]  = v[i];
    p[i] -= m[i];
    p1[i] = p[i];
  }

  normalizepvfv( u );

  S[0][1] = -u[2];
  S[0][2] =  u[1];
  S[1][0] =  u[2];
  S[1][2] = -u[0];
  S[2][0] = -u[1];
  S[2][1] =  u[0];
  
  for (i=0; i<3; i++)
    for (j=0; j<3; j++) {
      uij  = u[i]*u[j];
      uij1 = I[i][j] - uij;
      M[i][j] = uij + cosfi*uij1 + sinfi*S[i][j];
    }

  MatVecMult33f( M, p1, p );
  for (i=0; i<3; i++)
    p[i] += m[i];
}


double 
det3x3d(const double *v0, const double *v1, const double *v2)
{
  return (
	  v0[0]*v1[1]*v2[2] + v0[1]*v1[2]*v2[0] + v0[2]*v1[0]*v2[1] -
	  v2[0]*v1[1]*v0[2] - v2[1]*v1[2]*v0[0] - v2[2]*v1[0]*v0[1]
	  );
}


float 
det3x3f(const float *v0, const float *v1, const float *v2)
{
  return (
	  v0[0]*v1[1]*v2[2] + v0[1]*v1[2]*v2[0] + v0[2]*v1[0]*v2[1] -
	  v2[0]*v1[1]*v0[2] - v2[1]*v1[2]*v0[0] - v2[2]*v1[0]*v0[1]
	  );
}
