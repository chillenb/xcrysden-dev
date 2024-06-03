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
 * Source: $XCRYSDEN_TOPDIR/C/voronoi.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <GL/gl.h>
#include <stdio.h>
#include "xcGLparam.h"
#include "struct.h"
#include "wigner.h"
#include "xcfunc.h"


void
xcRenderVoronoi(void)
{
  int groupn = XCR_CELL_PC;
  register int ip, iv, i, j, k, l, m;
  float tr[3];
  int mode = WIGNER_PRIM;
  double (*vector)[4] = vec.prim;
  GLint     two_side_voronoi = GL_FALSE;
  GLint     shade_model;
  GLboolean two_side, cull_face;
  
  /* dummy is dummy so far */
  LoadVoronoiMaterial(VORONOI_DEFAULT);
  glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, ws_attrib.color );

  /* enable transparency */
  if ( ws_attrib.transparency ) {
    glEnable( GL_BLEND );
    glDepthMask( GL_FALSE );
    glBlendFunc( blend_voronoi.sfunc, blend_voronoi.dfunc );
  }

  /* SHADEMODEL */
  glGetIntegerv( GL_SHADE_MODEL, &shade_model );
  glShadeModel( GL_FLAT );

  /* unable GL_LIGHT_MODEL_TWO_SIDE */
  glGetBooleanv( GL_LIGHT_MODEL_TWO_SIDE, &two_side );
  glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, two_side_voronoi );

  /* enable CULL FACE */
  glGetBooleanv( GL_CULL_FACE, &cull_face );
  if (!cull_face) 
    glEnable( GL_CULL_FACE );

  /* RENDER - so far just one wigner-seitz*/
  if ( xcr.celltype == XCR_CONVCELL )  {
    mode   = WIGNER_CONV;
    vector = vec.conv;
    groupn = xcr.groupn;
  } 
  
  /* WIGNER_NODEMODE_SELECT */
  if ( ws_attrib.nodemode[mode] == WIGNER_NODEMODE_SELECT )
    for (i=0; i<ws_attrib.nnodes[mode]; i++) {
      for (j=0; j<3; j++)
	tr[j] = ws_attrib.nodes[mode][i][0]*vector[0][j] +
	  ws_attrib.nodes[mode][i][1]*vector[1][j] +
	  ws_attrib.nodes[mode][i][2]*vector[2][j];
      glPushMatrix();
        glTranslated( tr[0]-mx, tr[1]-my, tr[2]-mz);
	for (ip=0; ip < wsp.npoly; ip++) {
	  glBegin( GL_POLYGON );
	  for (iv=0; iv < wsp.nvert[ip]; iv++) {
	    glNormal3fv( wsp.norm[ip][iv] );
	    glVertex3fv( wsp.poly[ip][iv] );
	  }
	  glEnd();
	}
      glPopMatrix();
    }
  else if ( ws_attrib.nodemode[mode] == WIGNER_NODEMODE_EVERY ) {
    register int ni, nj, nk; 
    register float ii, jj, kk, trx, try, trz;
    /*
      don't render on every node, but just on outher regions,
      because inner regions will be covered
    */
    ni = nj = nk = 1;
    if (xcr.nunit[0] > 0 ) ni = 2;
    if (xcr.nunit[1] > 0 ) nj = 2;
    if (xcr.nunit[2] > 0 ) nk = 2;
    for (i=0; i<ni; i++)
      for (j=0; j<nj; j++)
	for (k=0; k<nk; k++)	  
	  for (l=0; l<xcr.npos[groupn]; l++) {
	    ii = jj = kk = 0.0;
	    if ( i > 0 ) ii = xcr.nunit[0];
	    if ( j > 0 ) jj = xcr.nunit[1];
	    if ( k > 0 ) kk = xcr.nunit[2];
	      trx = ii + xcr.cellpos[groupn][l][0];
	      try = jj + xcr.cellpos[groupn][l][1];
	      trz = kk + xcr.cellpos[groupn][l][2];
	      if ( trx < (double) xcr.nunit[0] + MINTOL &&
		   try < (double) xcr.nunit[1] + MINTOL &&
		   trz < (double) xcr.nunit[2] + MINTOL ) {
		/* from fractional coords -> Cartesian */
		for (m=0; m<3; m++)
		  tr[m] = trx*vector[0][m] + try*vector[1][m] + 
		    trz*vector[2][m];
		glPushMatrix();
		  glTranslated( tr[0]-mx, tr[1]-my, tr[2]-mz);
		  for (ip=0; ip < wsp.npoly; ip++) {
		    glBegin( GL_POLYGON );
		    for (iv=0; iv < wsp.nvert[ip]; iv++) {
		      glNormal3fv( wsp.norm[ip][iv] );
		      glVertex3fv( wsp.poly[ip][iv] );
		    }
		    glEnd();
		  }
		  glPopMatrix();
	      }
	  }
  }
  
  /* DISABLE TRANSPARENCY */  
  glDisable(GL_BLEND); 
  glDepthMask(GL_TRUE);  

  /* go back to old shademodel */
  glShadeModel( shade_model );

  /* go to old GL_LIGHT_MODEL_TWO_SIDE */
  /* obstaja funkcija LoadLightModelTwoSide(int type) */
  glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, two_side );

  /* disable cull face */
  if ( !cull_face ) 
    glDisable( GL_CULL_FACE );

  /* LOAD THE STRUCT MATERUALS BACK; 
   * this has to be changed in future; 
   * each render function will have to take care of it-self
   */
  LoadStructMaterial();
} 
