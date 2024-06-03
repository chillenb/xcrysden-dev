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
 * Source: $XCRYSDEN_TOPDIR/C/readPlvertex.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#define XC_CPP_NO_STDIO
#include "struct.h"
#include "isosurf.h"

void ReadPlvertex123( int type, int islide );

extern PLANEVERTEX ***plvertex;

void
ReadPlvertex123( int type, int islide ) {
  register int i, j;
  /* gridvertex[x][y][z];
     PLANE1 == vec0 x vec1;
     PLANE2 == vec2 x vec0;
     PLANE3 == vec1 x vec2;
  */

  /* 
     make islide within [0..nslide-1]
  */
  islide--;

  if ( type == ISOOBJ_PLANE1 ) {
    for(i=0; i<newgrd.nx; i++)
      for(j=0; j<newgrd.ny; j++) {
	plvertex[type][i][j].p[0] = gridvertex[i][j][islide].p.x;
	plvertex[type][i][j].p[1] = gridvertex[i][j][islide].p.y;
	plvertex[type][i][j].p[2] = gridvertex[i][j][islide].p.z;
	plvertex[type][i][j].val  = gridvertex[i][j][islide].val;
      }
  }
  else if ( type == ISOOBJ_PLANE2 ) {
    for(i=0; i<newgrd.nz; i++)
      for(j=0; j<newgrd.nx; j++) {
	plvertex[type][i][j].p[0] = gridvertex[j][islide][i].p.x;
	plvertex[type][i][j].p[1] = gridvertex[j][islide][i].p.y;
	plvertex[type][i][j].p[2] = gridvertex[j][islide][i].p.z;
	plvertex[type][i][j].val  = gridvertex[j][islide][i].val;
      }
  }
  else if ( type == ISOOBJ_PLANE3 ) {
    for(i=0; i<newgrd.ny; i++)
      for(j=0; j<newgrd.nz; j++) {
	plvertex[type][i][j].p[0] = gridvertex[islide][i][j].p.x;
	plvertex[type][i][j].p[1] = gridvertex[islide][i][j].p.y;
	plvertex[type][i][j].p[2] = gridvertex[islide][i][j].p.z;
	plvertex[type][i][j].val  = gridvertex[islide][i][j].val;
      }
  }
}
