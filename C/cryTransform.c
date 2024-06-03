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
 * Source: $XCRYSDEN_TOPDIR/C/cryTransform.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <togl.h>
#include <math.h>
#include "struct.h"

extern void vecMultMat(double major[4][4], double new[4][4], 
 		       double old[4][4], double vec[16]); 
extern void vecOldMat(double major[4][4], double old[4][4]);
extern void vecRotXYMat(double new[4][4], double cosfiX, double sinfiX,
			double cosfiY, double sinfiY);

int cryRotateXY(NEW_WIN_CONTEXT *wc, double fiX, double fiY);

int
cryRotateXY(NEW_WIN_CONTEXT *wc, double fiX, double fiY)
{
  double cosfiX, sinfiX;
  double cosfiY, sinfiY;
  
  cosfiX = cos(fiX);
  sinfiX = sin(fiX);
  cosfiY = cos(fiY);
  sinfiY = sin(fiY);
  
  vecRotXYMat( wc->vec.crdnew, cosfiX, sinfiX, cosfiY, sinfiY);
  vecOldMat( wc->vec.crdmajor, wc->vec.crdold );
  vecMultMat( wc->vec.crdmajor, wc->vec.crdnew, 
	      wc->vec.crdold, wc->vec.crdvec ); 
  return XC_OK;
}
