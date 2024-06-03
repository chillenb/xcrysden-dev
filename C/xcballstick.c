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
 * Source: $XCRYSDEN_TOPDIR/C/xcballstick.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "struct.h"

/* FUNCTION PROTOTYPES */
int MakeSticks1(int i,
		GLdouble *x1, GLdouble *y1, GLdouble *z1, 
		GLdouble *x2, GLdouble *y2, GLdouble *z2, 
		GLdouble *x3, GLdouble *y3, GLdouble *z3, 
		GLdouble *x4, GLdouble *y4, GLdouble *z4);
void MakeArcPoints(void);
int MakeSticks2(int i, int col, int flag);

/* extern function prototypes */
extern double dist3(double x, double y, double z);
extern void Rotate(double *x, double *y, double *z, double cosfi, double sinfi);

#define NROD 10
#define AB 0
#define BA 1

GLdouble rodx[NROD + 3];
GLdouble rody[NROD + 3];
GLdouble arcy[NROD + 1];
GLdouble arcz[NROD + 1];


/* Make a STICKS from (coor + *(iwksp + i))->x1 , .... 
 * routine is called from xcBallStick12D 
 *
 * routine return 1 on success
 *    &    return 0 on failure
 */
int
MakeSticks1(int i,
	    GLdouble *x1, GLdouble *y1, GLdouble *z1, 
	    GLdouble *x2, GLdouble *y2, GLdouble *z2, 
	    GLdouble *x3, GLdouble *y3, GLdouble *z3, 
	    GLdouble *x4, GLdouble *y4, GLdouble *z4)
{
  GLdouble corr, xcorr, ycorr, dis;
  GLdouble rb, fi, cosfi, sinfi;
  GLdouble xa, ya, za, xb1, yb1, zb1, xb2, yb2, zb2;
  short type;

  if ( (coor + i)->bondend == BOND_ATOM_TO_MIDBOND ) {
    xb1 = (coor + i)->x1;
    yb1 = (coor + i)->y1;
    zb1 = (coor + i)->z1;

    xb2 = (coor + i)->x2;
    yb2 = (coor + i)->y2;
    zb2 = (coor + i)->z2;
  } else {
    /*if ( (coor + i)->bondend == BOND_MIDBOND_TO_ATOM ) */
    xb1 = (coor + i)->x2;
    yb1 = (coor + i)->y2;
    zb1 = (coor + i)->z2;

    xb2 = (coor + i)->x1;
    yb2 = (coor + i)->y1;
    zb2 = (coor + i)->z1;
  }

  dis = dist3(xb2 - xb1, yb2 - yb1, zb2 - zb1);

  if ( dis < MINTOL ) {
    fprintf(stderr,"WARNING: zero bond distance\n");
    return 0;
  }
  rb = rball[(coor + i)->nat]; /*   !!!!!!!!!!!!!!!!!!   */

  /* first find an intersection of rod and a ball */
  xa = xb1 + rb * ( xb2 - xb1 ) / dis;
  ya = yb1 + rb * ( yb2 - yb1 ) / dis;
  za = zb1 + rb * ( zb2 - zb1 ) / dis;
  
  /*     ---------------------
   *     BALL-STICK CORRECTION
   *     ---------------------
   *     sphere on the right of the rod  
   */
  if (rb > rrod) corr = rb - sqrt(rb * rb - rrod * rrod);
  else corr = rb;
  if ( fabs(xb2 - xb1) > MINTOL ) {
    /* if ball is on the left of the rod */
    if ( xb1 < xa ) corr = -corr;
    fi = atan( (yb2 - yb1) / (xb2 - xb1) );
    
  } else { 
    /* if bond is nearly vertical than we must take care of this */
    if ( yb2 > yb1 ) fi = -PI12;
    else fi = PI12;
   }
  cosfi = cos(fi);
  sinfi = sin(fi);

  xa += corr * cosfi;
  ya += corr * sinfi;
  
  xcorr = -rrod * sinfi;
  ycorr =  rrod * cosfi;
    
  /* order of polygon drawing is 1, 2, 3, 4 
   * 
   *          1----------------4
   *   (ATOM) |                |
   *          2----------------3
   *                                  
   */   
  if ( (ABS(xb1-xa) > MINTOL) && xb1 < xa ) {
    type = AB;
  } else {
    type = BA;
  }

  if ( type == AB ) {
    /* atom on the left */
    *x1 = xa + xcorr;
    *y1 = ya + ycorr;
    *z1 = za;

    *x2 = xa - xcorr;
    *y2 = ya - ycorr;
    *z2 = za;

    *x3 = xb2 - xcorr;
    *y3 = yb2 - ycorr;  
    *z3 = zb2;

    *x4 = xb2 + xcorr;
    *y4 = yb2 + ycorr;
    *z4 = zb2;
  } else {
    *x4 = xa + xcorr;
    *y4 = ya + ycorr;
    *z4 = za;

    *x3 = xa - xcorr;
    *y3 = ya - ycorr;
    *z3 = za;

    *x2 = xb2 - xcorr;
    *y2 = yb2 - ycorr;  
    *z2 = zb2;

    *x1 = xb2 + xcorr;
    *y1 = yb2 + ycorr;
    *z1 = zb2;
  }

  return XC_OK;
}


void 
MakeArcPoints(void)
{
  int i;
 /* on one side of a bond there will be arc, calculate points for this arc */
  for(i = 0; i < NROD; i++)
    {
      arcy[i] = rrod * cos(i * PI / NROD);
      arcz[i] = rrod * sin(i * PI / NROD);
    }
  arcy[NROD] = -rrod;
  arcz[NROD] = 0.0;    
}


/* Make a STICKS2 from (coor + *(iwksp + i))->x1 , .... 
 * routine is called from xcBallStick12D 
 *
 * routine return 1 on success
 *    &    return 0 on failure
 */
int
MakeSticks2(int i, int col, int flag)
{
  /* i .... index of atom 
   * col .. index for color; 
   */
  int jj;
  double fac = 1.0;
  GLdouble corr, xcorr, ycorr, dis;
  GLdouble rb, fi, psi, xy, cosfi, sinfi, sinpsi;
  GLdouble xa, ya, za, xb1, yb1, zb1, xb2, yb2, zb2;
  double rodz = 0.0;
  short type;

  if ( (coor + i)->bondend == BOND_ATOM_TO_MIDBOND ) {
    xb1 = (coor + i)->x1;
    yb1 = (coor + i)->y1;
    zb1 = (coor + i)->z1;
    
    xb2 = (coor + i)->x2;
    yb2 = (coor + i)->y2;
    zb2 = (coor + i)->z2;
  } else {
    /*if ( (coor + i)->bondend == BOND_MIDBOND_TO_ATOM ) */
    xb1 = (coor + i)->x2;
    yb1 = (coor + i)->y2;
    zb1 = (coor + i)->z2;

    xb2 = (coor + i)->x1;
    yb2 = (coor + i)->y1;
    zb2 = (coor + i)->z1;
  }

  dis = dist3(xb2 - xb1, yb2 - yb1, zb2 - zb1);

  if ( dis < MINTOL ) {
    fprintf(stderr,"zero bond distance");
    return 0;
  }
  if ( flag != PIPEBALL ) {
    rb = rball[(coor + i)->nat]; /*!!!!!!!!!!!!!!!!!!*/
  } else {
    rb = rball[1];
  }
  
  /* first find an intersection of rod and a ball */
  xa = xb1 + rb * ( xb2 - xb1 ) / dis;
  ya = yb1 + rb * ( yb2 - yb1 ) / dis;
  za = zb1 + rb * ( zb2 - zb1 ) / dis;
  
  /*     ---------------------
   *     BALL-STICK CORRECTION
   *     ---------------------
   *     sphere on the right of the rod  
   */
  if (rb > rrod) corr = rb - sqrt(rb * rb - rrod * rrod);
  else corr = rb;
  if ( fabs(xb2 - xb1) > MINTOL ) {
    /* if ball is on the left of the rod */
    if ( xb1 < xa ) corr = -corr;
    fi = atan( (yb2 - yb1) /
	       (xb2 - xb1) );
    
  } else { 
    /* if bond is nearly vertical than we must take care of this */
    if ( yb2 > yb1 ) fi = -PI12;
    else fi = PI12;
   }  
  cosfi = cos(fi);
  sinfi = sin(fi);

  /* origin of the one side of a bond */
  xa += corr * cosfi;
  ya += corr * sinfi;

  xcorr = -rrod * sinfi;
  ycorr =  rrod * cosfi;
		 
  /* calculate angle between line [(xb2,yb2,zb2) - (xb1,yb1,zb1)] and it's
   * XY projection [ y = - (xcorr/ycorr) * x ]
   */
  xy = dist3(xb2 - xb1, yb2 - yb1, 0.0);
  if ( xy > MINTOL ) psi = atan( (zb2 - zb1) / xy );
  /* if psi > 0 --> (zb2 > zb1) arc on zb1 side */
  else psi = 0.0; /* this means we have a bond in pure Z-direction */
  /* if we have bond in pure XY direction psi = 0.0 also */
  sinpsi = sin(psi);

  /*==========================================================*/
  /* ATTENTION:  ARC IS ALWAYS ON THE LOWEST Z-SIDE OF A BOND */
  /*==========================================================*/
  if ( (ABS(xb1-xa) > MINTOL) && xb1 < xa ) {
    type = AB;
  } else {
    type = BA;
  }

  /* if psi > 0 -> bond in front of the ball, not needed to calculate an arc */
  if ( psi > 0.0 ) {
    GLdouble rodz = 0.0; /* this is just to satisfy Rotate() function */
    /* if ball is on the left of the rod we must subtract an arc,
     *  else we must add an arc */	
    if ( xb1 < xb2 ) fac = -1.0;
          
    if ( VPf.unibond) glColor3fv(unibondCol);
    else glColor3fv( atm.col[col] ); 
    glBegin(GL_POLYGON);	
      if ( type == AB ) {
	for(jj = 0; jj <= NROD; jj++) {
	  rodx[jj] = arcz[jj] * sinpsi * fac;
	  rody[jj] = arcy[jj];
	  Rotate((rodx + jj), (rody + jj), &rodz, cosfi, sinfi);
	  rodx[jj] += xa;
	  rody[jj] += ya;	
	  /* now draw a bond */
	  glVertex3d( rodx[jj], rody[jj], za );
	}
	glVertex3d( xb2 - xcorr, yb2 - ycorr, zb2 );  
	glVertex3d( xb2 + xcorr, yb2 + ycorr, zb2 );
      } else {
	glVertex3d( xb2 + xcorr, yb2 + ycorr, zb2 );
	glVertex3d( xb2 - xcorr, yb2 - ycorr, zb2 );  
	for(jj = NROD; jj >= 0; jj--) {
	  rodx[jj] = arcz[jj] * sinpsi * fac;
	  rody[jj] = arcy[jj];
	  Rotate((rodx + jj), (rody + jj), &rodz, cosfi, sinfi);
	  rodx[jj] += xa;
	  rody[jj] += ya;	
	  /* now draw a bond */
	  glVertex3d( rodx[jj], rody[jj], za );
	}
      }
    glEnd();

    glColor3f( 0.0, 0.0, 0.0 );
    glBegin(GL_LINE_LOOP);	
      for(jj = 0; jj <= NROD; jj++) {
	rodx[jj] = arcz[jj] * sinpsi * fac;
	rody[jj] = arcy[jj];
	Rotate((rodx + jj), (rody + jj), &rodz, cosfi, sinfi);
	rodx[jj] += xa;
	rody[jj] += ya;	
	/* now draw a bond */
	glVertex3d( rodx[jj], rody[jj], za );
      }
      glVertex3d( xb2 - xcorr, yb2 - ycorr, zb2 );  
      glVertex3d( xb2 + xcorr, yb2 + ycorr, zb2 );
    glEnd();	
  } else if ( psi < 0.0 ) {
    /* BOND BESIDE THE BALL */
    /* if ball is on the right of the rod, we must subtract an arc on
       the second end of a bond */
    if ( xb2 > xb1 ) fac = -1.0;

    if ( VPf.unibond) glColor3fv(unibondCol);
    else glColor3fv( atm.col[col] );	  	    
    glBegin(GL_POLYGON);
      if ( type == BA ) {
	for(jj = 0; jj <= NROD; jj++) {
	  rodx[jj] = arcz[jj] * sinpsi * fac;
	  rody[jj] = arcy[jj];
	  Rotate((rodx + jj), (rody + jj), &rodz, cosfi, sinfi);
	  rodx[jj] += xb2;
	  rody[jj] += yb2;	
	  /* now draw a bond */
	  glVertex3d( rodx[jj], rody[jj], zb2 );
	}
	glVertex3d( xa - xcorr, ya - ycorr, za );  
	glVertex3d( xa + xcorr, ya + ycorr, za );
      } else {
	glVertex3d( xa + xcorr, ya + ycorr, za );
	glVertex3d( xa - xcorr, ya - ycorr, za );  
	for(jj = NROD; jj >= 0; jj--) {
	  rodx[jj] = arcz[jj] * sinpsi * fac;
	  rody[jj] = arcy[jj];
	  Rotate((rodx + jj), (rody + jj), &rodz, cosfi, sinfi);
	  rodx[jj] += xb2;
	  rody[jj] += yb2;	
	  /* now draw a bond */
	  glVertex3d( rodx[jj], rody[jj], zb2 );
	}
      }
    glEnd();
    
    glColor3f( 0.0, 0.0, 0.0 );
    glBegin(GL_LINE_LOOP);
      for(jj = 0; jj <= NROD; jj++) {
	rodx[jj] = arcz[jj] * sinpsi * fac;
	rody[jj] = arcy[jj];
	Rotate((rodx + jj), (rody + jj), &rodz, cosfi, sinfi);
	rodx[jj] += xb2;
	rody[jj] += yb2;	
	/* now draw a bond */
	glVertex3d( rodx[jj], rody[jj], zb2 );
      }
      glVertex3d( xa - xcorr, ya - ycorr, za );  
      glVertex3d( xa + xcorr, ya + ycorr, za );
    glEnd();
  } else { 
    /* if ( psi == 0.0 ) */      
    /* order of polygon drawing is 1, 2, 3, 4 
     * 
     *         1----------------4        4---------------1
     *          |              |     or   |             | 
     *         2----------------3        3---------------2
     */   
        
    if ( VPf.unibond) glColor3fv(unibondCol);
    else glColor3fv( atm.col[col] );
    glBegin(GL_POLYGON);
      if ( type == AB ) {
	glVertex3d( xa + xcorr, ya + ycorr, za);
	glVertex3d( xa - xcorr, ya - ycorr, za);
	glVertex3d( xb2 - xcorr, yb2 - ycorr, zb2 );
	glVertex3d( xb2 + xcorr, yb2 + ycorr, zb2 );
      } else {
	glVertex3d( xb2 + xcorr, yb2 + ycorr, zb2 );
	glVertex3d( xb2 - xcorr, yb2 - ycorr, zb2 );
	glVertex3d( xa - xcorr, ya - ycorr, za);
	glVertex3d( xa + xcorr, ya + ycorr, za);
      }
    glEnd();

    glColor3f( 0.0, 0.0, 0.0 );
    glBegin(GL_LINE_LOOP);
      glVertex3d( xa + xcorr, ya + ycorr, za);
      glVertex3d( xa - xcorr, ya - ycorr, za);
      glVertex3d( xb2 - xcorr, yb2 - ycorr, zb2 );
      glVertex3d( xb2 + xcorr, yb2 + ycorr, zb2 );
    glEnd();
  }

  return XC_OK;
}

