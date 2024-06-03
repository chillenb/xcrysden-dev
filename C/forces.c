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
 * Source: $XCRYSDEN_TOPDIR/C/forces.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <GL/gl.h>
#include <stdio.h>
#include "xcGLparam.h"
#include "vector.h"
#include "struct.h"
#include "memory.h"
#include "xcfunc.h"

RenderVectors *forceVectors = NULL;
extern ForceVector FV;

void
MallocForceVectors(void) {  
  if (natoms>0){
    forceVectors = (RenderVectors *) 
      xcRealloc(forceVectors, sizeof(RenderVectors) * (natoms + 1));
  } else {
    forceVectors = (RenderVectors *) 
      xcRealloc(forceVectors, sizeof(RenderVectors));
  }
}

void
SetForceVectorsCoor( ForceVector *fvPtr, double (*fvec)[3], RenderVectors *rvec ) 
{
  register int i, j;
  
  if (fvPtr->rod_thickf < MINTOL)  fvPtr->rod_thickf = 1.0;
  if (fvPtr->arr_thickf < MINTOL ) fvPtr->arr_thickf = 1.0;
  if (fvPtr->arr_lenf   < MINTOL ) fvPtr->arr_lenf   = 0.1;

  for (i=1; i<=natoms; i++) 
    {

      rvec[i].vecthick = fvPtr->rod_thickf * rrod;
      rvec[i].arrthick = fvPtr->arr_thickf * rvec[i].vecthick;

      if ( distdv(fvec[i]) < MINTOL ) 
	{
	  for (j=0; j<3; j++) 
	    LOAD_NULL_V(3,rvec[i].coor[j]);

	  rvec[i].vecx  = 0.0; 
	  rvec[i].vecy  = 0.0; 
	  rvec[i].vecz  = 0.0;
	  rvec[i].vecfi = 0.0; 
	  rvec[i].vecl  = 0.0;

	  rvec[i].arrx  = 0.0; 
	  rvec[i].arry  = 0.0; 
	  rvec[i].arrz  = 0.0;
	  rvec[i].arrfi = 0.0; 
	  rvec[i].arrl  = 0.0;
	} 
      else 
	{   
	  for (j=0; j<3; j++) 
	    {
	      rvec[i].coor[j][0] = 0.0;
	      rvec[i].coor[j][1] = (1.0 - fvPtr->arr_lenf) * fvec[i][j];
	      rvec[i].coor[j][2] = fvec[i][j];
	    }
	  
	  GetSelCylinderPar(rvec[i].coor[0][1] - rvec[i].coor[0][0],
			    rvec[i].coor[1][1] - rvec[i].coor[1][0],
			    rvec[i].coor[2][1] - rvec[i].coor[2][0],
			    &rvec[i].vecx, &rvec[i].vecy, &rvec[i].vecz,
			    &rvec[i].vecfi, &rvec[i].vecl);
	  GetSelCylinderPar(rvec[i].coor[0][2] - rvec[i].coor[0][1], 
			    rvec[i].coor[1][2] - rvec[i].coor[1][1],
			    rvec[i].coor[2][2] - rvec[i].coor[2][1],
			    &rvec[i].arrx, &rvec[i].arry, &rvec[i].arrz,
			    &rvec[i].arrfi, &rvec[i].arrl);	  
	}
    }
}

void
xcRenderVectorForces(void)
{
  int i;

  glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, FV.color );

  for(i=1; i<=natoms; i++) {
    /* here the factor 2.0 is just for round-off error security */
    if ( forceVectors[i].vecl < 2.0*MINTOL ) continue;
    glPushMatrix();
      glTranslated( xat[i], yat[i], zat[i] );
      xcSolidVector( forceVectors[i] );
    glPopMatrix();
  }
  /* LOAD THE STRUCT MATERIALS BACK; 
   * this has to be changed in future; 
   * each render function will have to take care of it-self
   */
  LoadStructMaterial();
}
