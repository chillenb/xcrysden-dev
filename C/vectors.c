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
 * Source: $XCRYSDEN_TOPDIR/C/vectors.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#define  XC_CPP_NO_STDIO    /* look in struct.h file */
#include <math.h>
#include "struct.h"

/* --- function prototypes --- */
void vecMultMat(double major[4][4], double new[4][4], 
		double old[4][4], double vec[16]); 
void vecMatToVec(double major[4][4], double vec[16]);
void vecVecToMat(double vec[16], double major[4][4]);
void vecOldMat(double major[4][4], double old[4][4]);
void vecRotXYMat(double new[4][4], double cosfiX, double sinfiX,
		 double cosfiY, double sinfiY);
static void vecMatMultXY( GLdouble new[][4], GLdouble rotx[][4], GLdouble roty[][4] );
void vecRotX( GLdouble new[][4], GLdouble cosfiX, GLdouble sinfiX );
void vecRotY( GLdouble new[][4], GLdouble cosfiY, GLdouble sinfiY );
void vecRotZ( GLdouble new[][4], GLdouble cosfiZ, GLdouble sinfiZ );
void VecRotateX(double cosfi, double sinfi);
void VecRotateY(double cosfi, double sinfi);
void VecRotateZ(double cosfi, double sinfi);
void VecRotateXY(double cosfiX, double sinfiX, double cosfiY, double sinfiY);
void VecRotTo_XY(void);
void VecRotTo_XZ(void);
void VecRotTo_YZ(void);
static void VecProduct3df(double v0[], double v1[], float res[]);
void xcRotatefv(float fi, float u[]);
void VecRotTo_AB(void);
void VecRotTo_AC(void);
void VecRotTo_BC(void);


/* --- auxilproc.c --- */
extern int normalizepvfv(float *vec);
extern float distfv(float *v);
extern void VecProduct3f(float v0[], float v1[], float res[]);

/* --- 3D.c --- */
extern void LoadNull(GLdouble matrix[][4]);

/* ============================================== */
/* multiplication of matrix in COLOUMN MAJOR mode */
/* major[Kcol,Irow] = SUM(j){ new[Jcol,Irow] * old[Kcol,Jrow]
 *
 * also assign a vector xvec[] from major[4][4] in COLOUMN MAJOR mode
 */
void
vecMultMat(double major[4][4], double new[4][4], double old[4][4], double vec[16])
{
  int ivec = 0;
  int i,j,k;

  for(k=0; k<4; k++)
    for(i=0; i<4; i++) 
      {
	for(j=0; j<4;j++) 
	  major[k][i] = major[k][i] + new[j][i] * old[k][j];
	vec[ivec++] = major[k][i];
      }
}
 
 
void 
vecMatToVec(double major[4][4], double vec[16])
{
  int ivec = 0;
  int i,j;

  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      vec[ivec++] = major[i][j];
}


void 
vecVecToMat(double vec[16], double major[4][4] )
{
  int ivec = 0;
  int i,j;
  
  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      major[i][j] = vec[ivec++];    
}
  

void
vecOldMat(double major[4][4], double old[4][4])
{
  int i,j; 
  
  for(i=0; i<4; i++) 
    for(j=0; j<4;j++) {
      old[i][j] = major[i][j];
      major[i][j] = 0.0;
    }
}


/* In C 2D arrays are in row-major mode, but I'm dealing with coloumn major
 * mode; so whan initializing a matrix I must take care of above 
 * transposition
 */
void
vecRotXYMat( double new[4][4], double cosfiX, double sinfiX,
	     double cosfiY, double sinfiY)
{  
  /* 2.2 is just to mark where a cos() or sin() will go */
  GLdouble rotx[4][4] = { 
    { 1.0, 0.0, 0.0, 0.0 },
    { 0.0, 2.2, 2.2, 0.0 },
    { 0.0, 2.2, 2.2, 0.0 },
    { 0.0, 0.0, 0.0, 1.0 }
  };
  GLdouble roty[4][4] = {
    { 2.2, 0.0, 2.2, 0.0 },
    { 0.0, 1.0, 0.0, 0.0 },
    { 2.2, 0.0, 2.2, 0.0 },
    { 0.0, 0.0, 0.0, 1.0 }
  };

  rotx[1][1] = cosfiX;
  rotx[2][1] = -sinfiX;
  rotx[1][2] = sinfiX;
  rotx[2][2] = cosfiX;

  roty[0][0] = cosfiY;
  roty[2][0] = sinfiY;
  roty[0][2] = -sinfiY;
  roty[2][2] = cosfiY;
  /* first rotate around X, than around Y -> that means NEW=Y*X */
  vecMatMultXY( new, roty, rotx );
}


/* assign NewMat = rotx * roty */
static void
vecMatMultXY( GLdouble new[][4], GLdouble rotx[][4], GLdouble roty[][4] )
{
  int i,j,k;

  LoadNull( new );
  for(k=0; k<4; k++)
    for(i=0; i<4; i++)       
	for(j=0; j<4;j++) 
	  new[k][i] = new[k][i] + rotx[j][i] * roty[k][j];
}


void
vecRotX( GLdouble new[][4], GLdouble cosfiX, GLdouble sinfiX )
{
  new[0][0] = 1.0;
  new[1][0] = 0.0;
  new[2][0] = 0.0;
  new[3][0] = 0.0;

  new[0][1] = 0.0;
  new[1][1] = cosfiX;
  new[2][1] = -sinfiX;
  new[3][1] = 0.0;

  new[0][2] = 0.0;
  new[1][2] = sinfiX;
  new[2][2] = cosfiX;
  new[3][2] = 0.0;

  new[0][3] = 0.0;
  new[1][3] = 0.0;
  new[2][3] = 0.0;
  new[3][3] = 1.0;
}

void
vecRotY( GLdouble new[][4], GLdouble cosfiY, GLdouble sinfiY )
{
  new[0][0] = cosfiY;
  new[1][0] = 0.0;
  new[2][0] = sinfiY;
  new[3][0] = 0.0;

  new[0][1] = 0.0;
  new[1][1] = 1.0;
  new[2][1] = 0.0;
  new[3][1] = 0.0;

  new[0][2] = -sinfiY;
  new[1][2] = 0.0;
  new[2][2] = cosfiY;
  new[3][2] = 0.0;

  new[0][3] = 0.0;
  new[1][3] = 0.0;
  new[2][3] = 0.0;
  new[3][3] = 1.0;
}


void
vecRotZ( GLdouble new[][4], GLdouble cosfiZ, GLdouble sinfiZ )
{
  new[0][0] = cosfiZ;
  new[1][0] = -sinfiZ;
  new[2][0] = 0.0;
  new[3][0] = 0.0;

  new[0][1] = sinfiZ;
  new[1][1] = cosfiZ;
  new[2][1] = 0.0;
  new[3][1] = 0.0;

  new[0][2] = 0.0;
  new[1][2] = 0.0;
  new[2][2] = 1.0;
  new[3][2] = 0.0;

  new[0][3] = 0.0;
  new[1][3] = 0.0;
  new[2][3] = 0.0;
  new[3][3] = 1.0;
}


void 
VecRotateX(double cosfi, double sinfi)
{
  vecRotX( vec.crdnew, cosfi, sinfi);
  vecOldMat( vec.crdmajor, vec.crdold );
  vecMultMat( vec.crdmajor, vec.crdnew, vec.crdold, vec.crdvec ); 
}


void 
VecRotateY(double cosfi, double sinfi)
{
  vecRotY( vec.crdnew, cosfi, sinfi);
  vecOldMat( vec.crdmajor, vec.crdold );
  vecMultMat( vec.crdmajor, vec.crdnew, vec.crdold, vec.crdvec ); 
}


void 
VecRotateZ(double cosfi, double sinfi)
{
  vecRotZ( vec.crdnew, cosfi, sinfi);
  vecOldMat( vec.crdmajor, vec.crdold );
  vecMultMat( vec.crdmajor, vec.crdnew, vec.crdold, vec.crdvec ); 
}


void 
VecRotateXY(double cosfiX, double sinfiX, double cosfiY, double sinfiY)
{
  vecRotXYMat( vec.crdnew, cosfiX, sinfiX, cosfiY, sinfiY);
  vecOldMat( vec.crdmajor, vec.crdold );
  vecMultMat( vec.crdmajor, vec.crdnew, vec.crdold, vec.crdvec ); 
}


/* orientate to XY plane */
void
VecRotTo_XY(void) {
  register int i, k, ivec;
  GLdouble v[4][4] = {
    {1.0, 0.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 0.0, 1.0}
  };

  ivec = 0;
  for(k=0; k<4; k++)
    for(i=0; i<4; i++) 
      {
	vec.crdmajor[k][i] = v[k][i];
	vec.crdvec[ivec++] = vec.crdmajor[k][i];
      }
}

/* orientate to XZ plane (90 degree rotation around X axis) */
void
VecRotTo_XZ(void) {
  register int i, k, ivec;
  GLdouble v[4][4] = {
    {1.0, 0.0, 0.0, 0.0},
    {0.0, 0.0,-1.0, 0.0},
    {0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 1.0}
  };

  ivec = 0;
  for(k=0; k<4; k++)
    for(i=0; i<4; i++) 
      {
	vec.crdmajor[k][i] = v[k][i];
	vec.crdvec[ivec++] = vec.crdmajor[k][i];
      }
}

/* orientate to YZ plane (-90 degree rotation around Y axis) */
void
VecRotTo_YZ(void) {
  register int i, k, ivec;
  GLdouble v[4][4] = {
    { 0.0, 0.0, 1.0, 0.0},
    { 0.0, 1.0, 0.0, 0.0},
    {-1.0, 0.0, 0.0, 0.0},
    { 0.0, 0.0, 0.0, 1.0}
  };

  ivec = 0;
  for(k=0; k<4; k++)
    for(i=0; i<4; i++) 
      {
	vec.crdmajor[k][i] = v[k][i];
	vec.crdvec[ivec++] = vec.crdmajor[k][i];
      }
}

static  
void VecProduct3df(double v0[], double v1[], float res[])
{
  res[0] = (float) (v0[1] * v1[2] - v1[1] * v0[2]);
  res[1] = (float) (v0[2] * v1[0] - v1[2] * v0[0]);
  res[2] = (float) (v0[0] * v1[1] - v1[0] * v0[1]);
}

void
xcRotatefv(float fi, float u[]) {
  register int i, j, ivec;
  float cosfi, sinfi;
  float S[3][3], uut[3][3], uutId[3][3], M[4][4];

  if (distfv(u) > MINTOL) { 
    cosfi = cos(fi);
    sinfi = sin(fi);

    for (i=0; i<3; i++)
      for (j=0; j<3; j++) {
	uut[i][j] = u[i] * u[j];
	uutId[i][j] = -uut[i][j];
	if (i==j) uutId[i][j] += 1.0;
      }
    
    S[0][0] =  0.0;
    S[0][1] = -u[2];
    S[0][2] =  u[1];
    S[1][0] =  u[2];
    S[1][1] =  0.0;
    S[1][2] = -u[0];
    S[2][0] = -u[1];
    S[2][1] =  u[0];
    S[2][2] =  0.0;

    /* in OpenGL, matrix are written in column major mode, not in 
       row-major mode as is usually in C; take care about that */
    for(i=0; i<3; i++)
      for(j=0; j<3; j++) 
	M[j][i] = uut[i][j] + cosfi * uutId[i][j] + sinfi * S[i][j];  
    M[3][0] = M[3][1] = M[3][2] = M[0][3] = M[1][3] = M[2][3] = 0.0;
    M[3][3] = 1.0;    
  } else {
    for (i=0; i<4; i++)
      for (j=0; j<4; j++) {
	M[i][j] = 0.0;
	M[i][i] = 1.0;
      }
  }
  ivec = 0;
  for(i=0; i<4; i++)
    for(j=0; j<4; j++) 
      {
	vec.crdmajor[i][j] = (GLdouble) M[i][j];
	vec.crdvec[ivec++] = (GLdouble) vec.crdmajor[i][j];
      }
}
  
void
VecRotTo_AB(void) {
  /* normal must be oriented to Z direction, AxB=C */
  float fi, n[3], u[3], z[3] = {0.0, 0.0, 1.0};

  /* XxY = Z */ 
  if (xcr.celltype == XCR_NOCELL || xcr.celltype == XCR_PRIMCELL) 
    VecProduct3df(vec.prim[0], vec.prim[1], n);
  else
    VecProduct3df(vec.conv[0], vec.conv[1], n);

  VecProduct3f(n, z, u);
  normalizepvfv( u );
  
  fi = acos( n[2] / distfv(n) );
  xcRotatefv( fi, u );  
}

void
VecRotTo_AC(void) {
  /* normal must be oriented to Z direction, AxC=-B; CxA=B */
  float fi, n[3], u[3], z[3] = {0.0, 0.0, 1.0};

  /* ZxX = Y */ 
  if (xcr.celltype == XCR_NOCELL || xcr.celltype == XCR_PRIMCELL) 
    VecProduct3df(vec.prim[2], vec.prim[0], n);
  else
    VecProduct3df(vec.conv[2], vec.conv[0], n);
    
  VecProduct3f(n, z, u);
  normalizepvfv( u );
  
  fi = acos( n[2] / distfv(n) );
  xcRotatefv( fi, u );  
}

void
VecRotTo_BC(void) {
  /* normal must be oriented to Z direction, BxC=-A */
  float fi, n[3], u[3], z[3] = {0.0, 0.0, 1.0};

  /* YxZ = X */ 
  if (xcr.celltype == XCR_NOCELL || xcr.celltype == XCR_PRIMCELL) 
    VecProduct3df(vec.prim[1], vec.prim[2], n);
  else
    VecProduct3df(vec.conv[1], vec.conv[2], n);

  VecProduct3f(n, z, u);
  normalizepvfv( u );
  
  fi = acos( n[2] / distfv(n) );
  xcRotatefv( fi, u );  
}






