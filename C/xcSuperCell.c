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
* Source: $XCRYSDEN_TOPDIR/C/xcSuperCell.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <togl.h>
#include <tk.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "bz.h"
#include "vector.h"
#include "xcGLparam.h"
#include "xcfunc.h"

#define SC_NULL 0
#define SC_INIT 1

void (*xcSuperCell)(void) = 0x0;
void xcSuperCellFunc( void (*Func)(void) );
int XC_SuperCellCmd(ClientData clientData,Tcl_Interp *interp,
		    int argc, const char *argv[]);
void SuperCell_RenderPVec(void);
void SuperCell_RenderPVec_and_SCVec(void);
void SetUnitCellCage( double vec[4][4], CellCage *cage );
static void SetVectorsCoor( double vec[4][4], RenderVectors rvec[3] );
static void SuperCell_RenderVectors( RenderVectors vec[3] );
static void SuperCell_VectorProduct(GLfloat vp[3], double v1[], double v2[]);


/*****************************************************************************/
void
xcSuperCellFunc( void (*Func)(void) )
{
  xcSuperCell = Func;
}
/*****************************************************************************/


/*****************************************************************************
 * xc_supercell <toglName> init                                              *
 *                         testit                                            *
 *                         clear                                             *
 *****************************************************************************/
int 
XC_SuperCellCmd(ClientData clientData,Tcl_Interp *interp,
		int argc, const char *argv[])
{
  struct Togl *togl;
  static int sc_state = SC_NULL;

  if ( argc != 3 ) {
    Tcl_SetResult(interp, "Usage: xc_supercell <toglName> (init|update|clear)", TCL_STATIC);
    return TCL_ERROR;
  }

  /* find togl associated with toglName */
  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "couldn't find %s togl widget", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* xc_supercell available just for periodic systems */
  if ( xcr.dim < 1 ) {
    Tcl_SetResult(interp, "\"xc_supercell\" available just for periodic systems", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( strcmp( argv[2], "init" ) == 0 ) {
    sc_state      = SC_INIT;
    VPf.supercell = GL_TRUE;
    SetVectorsCoor( vec.prim, rnd_pvec );
    /* if periodicity of system is greater then 1 then display cellcage */
    if (xcr.dim > 1) SetUnitCellCage( vec.prim, &prim_cage );
    xcSuperCellFunc( SuperCell_RenderPVec );
  }
  else if ( strcmp( argv[2], "testit" ) == 0 ) {
    if ( sc_state == SC_NULL ) {
      Tcl_SetResult(interp, "\"xc_supercell init\" should be called before \"xc_supercell testit\"", TCL_STATIC);
      return TCL_ERROR;
    }
    /* by that time, the primitive vectors should already be new supercell
       vectors */
    SetVectorsCoor( vec.prim, rnd_scvec );
    /* if periodicity of system is greater then one then display cage */
    if (xcr.dim > 1) SetUnitCellCage( vec.prim, &sc_cage );
    xcSuperCellFunc( SuperCell_RenderPVec_and_SCVec );
  } else if ( strcmp( argv[2], "clear" ) == 0 ) {
    sc_state      = SC_NULL;
    VPf.supercell = GL_FALSE;
  } else {
    char rss[1024];
    snprintf(rss, sizeof(rss),"unknown option \"%s\", must be one of init, testit or clear", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* --- RENDER --- */
  Togl_PostRedisplay(togl);
  return TCL_OK;
}


void
SuperCell_RenderPVec(void)
{
  GLint     shade_model;
  GLboolean two_side, cull_face;

  glPushMatrix();
  glTranslated( -mx, -my, -mz );
  LoadCageOrVecMaterial( GLPAR_PRIMVEC );
  SuperCell_RenderVectors( rnd_pvec );

  /* if periodicity of system is greater then one then display cage */
  if (xcr.dim > 1) { 
    LoadCageOrVecMaterial( GLPAR_PRIMCAGE );
    SetCageOGLState( xcr.dim, blend_cellcage, 
		     &shade_model, &two_side, &cull_face );
    xcSolidCage( prim_cage );
    DisableCageOGLState( xcr.dim, shade_model, two_side, cull_face );
  }
  /* LOAD THE STRUCT MATERUALS BACK; 
   * this has to be changed in future; 
   * each render function will have to take care of it-self
   */
  LoadStructMaterial();

  glPopMatrix();
}


void
SuperCell_RenderPVec_and_SCVec(void)
{
  GLint     shade_model;
  GLboolean two_side, cull_face;

  glPushMatrix();
  glTranslated( -mx, -my, -mz );
  LoadCageOrVecMaterial( GLPAR_PRIMVEC );
  SuperCell_RenderVectors( rnd_pvec );
  LoadCageOrVecMaterial( GLPAR_SCVEC );
  SuperCell_RenderVectors( rnd_scvec );

  /* if periodicity of system is greater then one then display cage */
  if (xcr.dim > 1) {
    LoadCageOrVecMaterial( GLPAR_PRIMCAGE );
    SetCageOGLState( xcr.dim, blend_cellcage, 
		     &shade_model, &two_side, &cull_face );
    xcSolidCage( prim_cage );
    LoadCageOrVecMaterial( GLPAR_SCCAGE );
    xcSolidCage( sc_cage );
    DisableCageOGLState( xcr.dim, shade_model, two_side, cull_face );
  }

  /* LOAD THE STRUCT MATERUALS BACK; 
   * this has to be changed in future; 
   * each render function will have to take care of it-self
   */
  LoadStructMaterial();

  glPopMatrix();
}


void
SetUnitCellCage( double vec[4][4], CellCage *cage )
{
  register int i, j, k, npoly;
  int vindex[6][4][3] = {
    { /* bottom face */ 
      {0, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0} },
    { /* front face */
      {1, 0, 0}, {1, 1, 0}, {1, 1, 1}, {1, 0, 1} },
    { /* top face */
      {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, {0, 0, 1} },
    { /* back face */
      {0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0} },
    { /* left face */
      {0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {0, 0, 1} },
    { /* right face */
      {1, 1, 1}, {1, 1, 0}, {0, 1, 0}, {0, 1, 1} }
  };
  int nindex[6][2] = {
    /* bottom face */
    {1, 0},
    /* front face */
    {1, 2},
    /* top face */
    {0, 1},
    /* back face */
    {2, 1},
    /* left face */
    {0, 2},
    /* right face */
    {2, 0}
  };

  npoly = 6;
  if ( xcr.dim == 2) npoly = 1;
  if ( VPf.isospacesel3D ) npoly = 6;

  for (i=0; i<npoly; i++)
    for (j=0; j<4; j++) {
      for (k=0; k<3; k++)
	cage->vertex[i][j][k] = 
	  (GLfloat) vindex[i][j][0] * (GLfloat) vec[0][k] +
	  (GLfloat) vindex[i][j][1] * (GLfloat) vec[1][k] +
	  (GLfloat) vindex[i][j][2] * (GLfloat) vec[2][k];

      SuperCell_VectorProduct( &cage->normal[i][j][0], 
			       vec[nindex[i][0]], vec[nindex[i][1]] );
      normalizepvfv( &cage->normal[i][j][0] );
    }
}


static void
SetVectorsCoor( double vec[4][4], RenderVectors rvec[3] ) 
{
  register int i, j;

  for (i=0; i<xcr.dim; i++) {
    for (j=0; j<3; j++) {
      rvec[i].coor[j][0] = 0.0;
      rvec[i].coor[j][1] = (1.0 - VECTOR_ARROWSIZE) * vec[i][j];
      rvec[i].coor[j][2] = vec[i][j];
    }

    rvec[i].vecthick = VECTOR_THICKF * rrod;
    rvec[i].arrthick = VECTOR_ARRTHICKF * rvec[i].vecthick;

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


static void
SuperCell_RenderVectors( RenderVectors vec[3] )
{
  register int i;

  for (i=0; i<xcr.dim; i++)
    xcSolidVector( vec[i] );
}  


static void
SuperCell_VectorProduct(GLfloat vp[3], double v1[], double v2[])
{

  vp[0] = (float) v1[1] * v2[2] - v2[1] * v1[2];
  vp[1] = (float) v1[2] * v2[0] - v2[2] * v1[0];
  vp[2] = (float) v1[0] * v2[1] - v2[0] * v1[1];
}
