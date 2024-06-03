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
 * Source: $XCRYSDEN_TOPDIR/C/3D.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <math.h>
#include <stdio.h>
#include "struct.h"
#include "3D.h"
#include "xcfunc.h"

/* --- function prototypes --- */
void LoadNull(GLdouble matrix[][4]);
void LoadNull44One(GLdouble matrix[][4]); /* not needed yet */
void LoadIdentity(GLdouble matrix[][4]);
void GetMajorMat(void);
void MajorMatToVec(void);
void GetOldMat(void);
void GetRotXYMat( GLdouble fiX, GLdouble fiY );
static void MatMultXY( GLdouble rotx[][4], GLdouble roty[][4] );
void GetRotX( GLdouble fiX );
void GetRotY( GLdouble fiY );
void GetRotZ( GLdouble fiZ );
/* void From2Dto3D(void); ---moved to --> trash.c */
void MakeCylinderCoor(void);
static void GetCylinderPar(int numb, double x21, double y21, double z21);
/* void From3Dto2D(void); ---moved to --> trash.c */
/*static void UpdateCoor( double *x, double *y, double *z);*/

/* --- auxilproc.c --- */
extern double dist3(double x, double y, double z);


/*****************************************************************************/
/* razlaga funkcij, ki se uporabljajo pri prehodu 2D -> 3D.
   Prepisati je potrebno "AtomBond coor" -> "struct 3D". 
   Pri 3D koordinat ne speminjam ampak zasledujem matriko, s katero vsakic
   pomnozim coordinate */

/* From2Dto3D()  ... koordinate iz 2D oblike v 3D obliko
 * From3Dto2D()  ... obratno kot zgoraj
 * GetNewMatrix()... dobis novo rotacijsko matriko, ki ustreza 
 *                    zeljeni rotaciji
 * GetMajorMat() ... Major Matrix je matrika s katero pomozim koordinate,
 *                   da dobim zeljeno sliko (izracuna Major matriko)
 * GetOldMat()   ... Major Matriko prepise v Old Matriko; MajorMat pa
 *                   v nulto matriko;
 * LoadNull(GLdouble matrix[][])     ... load null matrix
 * LoadNull44One                     ... to all elements assign null,
 *                                       except to last elemnt assign 1.0
 * LoadIdentity(GLdouble matrix[][]) ... load identity matrix
 */
	
void
LoadNull(GLdouble matrix[][4])
{
  int i,j;

  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      matrix[i][j] = 0.0;
}


/* assign null to all elements, except to [3][3] assign 1.0 */
void
LoadNull44One(GLdouble matrix[][4])
{
  int i,j;

  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      matrix[i][j] = 0.0;
  matrix[3][3] = 1.0;  
}


void
LoadIdentity(GLdouble matrix[][4])
{
  int i,j;

  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      if ( i != j) matrix[i][j] = 0.0;
      else matrix[i][j] = 1.0;
}


/* ============================================== */
/* multiplication of matrix in COLOUMN MAJOR mode */
/* MajorMat[Kcol,Irow] = SUM(j){ NewMat[Jcol,Irow] * OldMat[Kcol,Jrow]
 *
 * also assign a vector xcGLvec[] from MajorMat[4][4] in COLOUMN MAJOR mode
 */
void
GetMajorMat(void)
{
  int ivec = 0;
  int i,j,k;

  for(k=0; k<4; k++)
    for(i=0; i<4; i++) 
      {
	for(j=0; j<4;j++) 
	  MajorMat[k][i] = MajorMat[k][i] + NewMat[j][i] * OldMat[k][j];
	xcGLvec[ivec++] = MajorMat[k][i];
      }
}
 
 
void 
MajorMatToVec(void)
{
  int ivec = 0;
  int i,j;

  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      xcGLvec[ivec++] = MajorMat[i][j];
}


void
GetOldMat(void)
{
  int i,j;
  
  for(i=0; i<4; i++) 
    for(j=0; j<4;j++) {
      OldMat[i][j] = MajorMat[i][j];
      MajorMat[i][j] = 0.0;
    }
}


/* In C 2D arrays are in row-major mode, but I'm dealing with coloumn major
 * mode; so whan initializing a matrix I must take care of above 
 * transposition
 */
void
GetRotXYMat( GLdouble fiX, GLdouble fiY )
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

  rotx[1][1] = cos(fiX);
  rotx[2][1] = -sin(fiX);
  rotx[1][2] = sin(fiX);
  rotx[2][2] = cos(fiX);

  roty[0][0] = cos(fiY);
  roty[2][0] = sin(fiY);
  roty[0][2] = -sin(fiY);
  roty[2][2] = cos(fiY);

  MatMultXY( rotx, roty );
}


/* assign NewMat = rotx * roty */
static void
MatMultXY( GLdouble rotx[][4], GLdouble roty[][4] )
{
  int i,j,k;

  LoadNull( NewMat );
  for(k=0; k<4; k++)
    for(i=0; i<4; i++)       
	for(j=0; j<4;j++) 
	  NewMat[k][i] = NewMat[k][i] + rotx[j][i] * roty[k][j];
}


void
GetRotX( GLdouble fiX )
{
  NewMat[0][0] = 1.0;
  NewMat[1][0] = 0.0;
  NewMat[2][0] = 0.0;
  NewMat[3][0] = 0.0;

  NewMat[0][1] = 0.0;
  NewMat[1][1] = cos(fiX);
  NewMat[2][1] = -sin(fiX);
  NewMat[3][1] = 0.0;

  NewMat[0][2] = 0.0;
  NewMat[1][2] = sin(fiX);
  NewMat[2][2] = cos(fiX);
  NewMat[3][2] = 0.0;

  NewMat[0][3] = 0.0;
  NewMat[1][3] = 0.0;
  NewMat[2][3] = 0.0;
  NewMat[3][3] = 1.0;
}

void
GetRotY( GLdouble fiY )
{
  NewMat[0][0] = cos(fiY);
  NewMat[1][0] = 0.0;
  NewMat[2][0] = sin(fiY);
  NewMat[3][0] = 0.0;

  NewMat[0][1] = 0.0;
  NewMat[1][1] = 1.0;
  NewMat[2][1] = 0.0;
  NewMat[3][1] = 0.0;

  NewMat[0][2] = -sin(fiY);
  NewMat[1][2] = 0.0;
  NewMat[2][2] = cos(fiY);
  NewMat[3][2] = 0.0;

  NewMat[0][3] = 0.0;
  NewMat[1][3] = 0.0;
  NewMat[2][3] = 0.0;
  NewMat[3][3] = 1.0;
}


void
GetRotZ( GLdouble fiZ )
{
  NewMat[0][0] = cos(fiZ);
  NewMat[1][0] = -sin(fiZ);
  NewMat[2][0] = 0.0;
  NewMat[3][0] = 0.0;

  NewMat[0][1] = sin(fiZ);
  NewMat[1][1] = cos(fiZ);
  NewMat[2][1] = 0.0;
  NewMat[3][1] = 0.0;

  NewMat[0][2] = 0.0;
  NewMat[1][2] = 0.0;
  NewMat[2][2] = 1.0;
  NewMat[3][2] = 0.0;

  NewMat[0][3] = 0.0;
  NewMat[1][3] = 0.0;
  NewMat[2][3] = 0.0;
  NewMat[3][3] = 1.0;
}


void
MakeCylinderCoor(void)
{
  int i, ii;

  /* first cylinder's coordinates for BONDS */
  for(i=1; i <= nbonds; i++)
    {      
      (coor3D + i)->flag = BOND;
      (coor3D + i)->bondFlag = *(bondend + i);
      (coor3D + i)->nat = *(natbond + i);
      (coor3D + i)->sqn = *(sqnbond + i);
      /* (coor3D + i)->list is not copied !!!! */
      (coor3D + i)->x1 = *(xbond + i);
      (coor3D + i)->y1 = *(ybond + i);
      (coor3D + i)->z1 = *(zbond + i);
      GetCylinderPar(i, *(xbond2 + i) - *(xbond + i),
		     *(ybond2 + i) - *(ybond + i),
		     *(zbond2 + i) - *(zbond + i) );
    }

  /* now, cylinder's coordinates for FRAMES, if there is any frame */
  for(i=nbonds + 1; i<= nbonds + nframes; i++)
    {
      ii = i - nbonds;
      (coor3D + i)->flag = FRAME;
      (coor3D + i)->nat = *(frametype + ii);
      (coor3D + i)->sqn = ii;      
      /* (coor3D + i)->list is not copied !!!! */
      (coor3D + i)->x1 = *(xframe + ii);
      (coor3D + i)->y1 = *(yframe + ii);
      (coor3D + i)->z1 = *(zframe + ii);
      GetCylinderPar(i, *(xframe2 + ii) - *(xframe + ii),
		     *(yframe2 + ii) - *(yframe + ii),
		     *(zframe2 + ii) - *(zframe + ii) );
    }
}


/* calculate fibond, xrvb, yrvb, zrvb, bondl */
static void 
GetCylinderPar(int numb, double x21, double y21, double z21)
{  
  double mintol = -1.0 + 1.0e-8;
  
  /* rotated vector is cross product between (0,0,1) & vector submited 
     to GetCylinderPar */
  (coor3D + numb)->bondl = dist3(x21, y21, z21);  

  /* normalize vectror (x21,y21,z21) */
  normalizepv(&x21, &y21, &z21);

  /* maybe (x21, y21, z21) vector is (0.0, 0.0, 1.0) or
                                     (0.0, 0.0,-1.0) */
  if (z21 > mintol) {
    (coor3D + numb)->xrvb = -y21;
    (coor3D + numb)->yrvb = x21;
    (coor3D + numb)->zrvb = 0.0;
  } else {
    (coor3D + numb)->xrvb = 1.0;
    (coor3D + numb)->yrvb = 0.0;
    (coor3D + numb)->zrvb = 0.0;
  }
  (coor3D+ numb)->fibond = RAD2DEG * acos(z21);

  if (ABS(z21) > 1.0000000001) {
    fprintf(stderr,"Z-component of normalized vector greater than 1.0\n");
    fprintf(stderr,"Vector: %f %f %f\n", x21, y21, z21);
    return;
  }
  /* Testing ... */  
  /*
  if ( ABS((coor3D + numb)->xrvb) < 1e-8 && ABS((coor3D + numb)->yrvb) < 1e-8 ) {
    (coor3D + numb)->xrvb = 1.0;
    (coor3D + numb)->yrvb = 0.0;
    (coor3D + numb)->zrvb = 0.0;
    (coor3D+ numb)->fibond = 0.0;
  } 
  */  
}

/*
static void 
UpdateCoor( double *x, double *y, double *z)
{
  double x1, y1, z1;

  x1 = MajorMat[0][0] * *x + MajorMat[1][0] * *y + MajorMat[2][0] * *z;
  y1 = MajorMat[0][1] * *x + MajorMat[1][1] * *y + MajorMat[2][1] * *z;
  z1 = MajorMat[0][2] * *x + MajorMat[1][2] * *y + MajorMat[2][2] * *z;
  *x = x1;
  *y = y1;
  *z = z1;
}
*/
