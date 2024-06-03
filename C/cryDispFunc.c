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
* Source: $XCRYSDEN_TOPDIR/C/cryDispFunc.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <togl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "xcGLparam.h"
#include "struct.h"
#include "isosurf.h"
#include "memory.h"
#include "molsurf.h"
#include "vector.h"
#include "bz.h"
#include "wigner.h"
#include "lighting.h"
#include "xcfunc.h"

extern XCantialias    antialias;


extern NEW_WIN_CONTEXT *FindWinContextByTogl(struct Togl *togl);

void cryDisplayFS(struct Togl *togl);
static GLfloat _getSize (NEW_WIN_CONTEXT *wc);
static void _cryDisplayFS(struct Togl *togl, NEW_WIN_CONTEXT *wc);
void cryRenderSurface( NEW_WIN_CONTEXT *wc , int fs_mode);
void cryRenderRecCell( NEW_WIN_CONTEXT *wc , int fs_mode);
static void RenderBZ(GLenum what);

void crySet_Trans_Shade_Lightmodel( int transp, GLenum shade, 
				    GLint two_side, GLboolean cull );

/*GLfloat line_color[4] = { 1.0, 1.0, 1.0, 1.0 };*/

/*
  assign
  cry_dispfunc <togl> fermisurf|incomming ...
*/
int 
CRY_DispFuncCmd(ClientData clientData,Tcl_Interp *interp,
		int argc, char *argv[])
{
  struct Togl *togl;
  NEW_WIN_CONTEXT *wc;

  if ( argc != 3 && argc != 5 ) {
    Tcl_SetResult(interp, "Usage: cry_dispfunc <togl> fermisurf ?-antialias 0|1?", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  wc = FindWinContextByTogl( togl );
  LoadIdentity( wc->vec.crdmajor );
  vecMatToVec(  wc->vec.crdmajor, wc->vec.crdvec );

  /* FERMISURFACE */
  if ( strncmp(argv[2], "fermis", 6) == 0 ) {
    wc->xcDisplay = cryDisplayFS;
    /* next lines are just for now */
    glEnable( GL_DITHER );
    glEnable( GL_DEPTH_TEST );      
    glEnable( GL_LIGHTING );
    glEnable( GL_LIGHT0 );
    glDepthFunc(GL_LEQUAL);
    /*glClearDepth( 1.0 );*/
    glClearColor( wc->bg[0], wc->bg[1], wc->bg[2], wc->bg[3] ); /* T.K. */
    LoadBlendfunc_And_Frontface();
    LoadIsoMaterial( MAT_ONELEVEL | MAT_ISOSURF );

    /*LoadLights();*/
  }

  if (argc == 5) {
    int anti_alias;
    if ( strncmp(argv[3], "-antialias", 6) == 0 ) {
      if ( Tcl_GetInt(interp, argv[4], &anti_alias) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\"", argv[4]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      wc->VPf.antialias = anti_alias;
    } else {
      char rss[1024];
      snprintf(rss, sizeof(rss), 
	       "unknown option %s, must be -antialias", argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  return TCL_OK;
}


static GLfloat _getSize (NEW_WIN_CONTEXT *wc)
{
  wc->MVf.structsize = sqrt( wc->ss.maxX*wc->ss.maxX +
			     wc->ss.maxY*wc->ss.maxY +
			     wc->ss.maxZ*wc->ss.maxZ );  
  return wc->MVf.structsize / wc->tr.zoom;
}

/* 
 * when we enter this procedure we just Made a Display, that means that 
 * everything must be made before that
 */
void
cryDisplayFS(struct Togl *togl)
{
  NEW_WIN_CONTEXT *wc;  
  wc = FindWinContextByTogl( togl );

  /* debugging */
  if (wc == (NEW_WIN_CONTEXT *)NULL) return;

  /*
    line_color[0] = 1.0 - wc->bg[0];
    line_color[1] = 1.0 - wc->bg[1];
    line_color[2] = 1.0 - wc->bg[2];
  */

  if ( ! wc->VPf.antialias ) 
    {
      /* 
	 NO ANTIALIASING 
      */
      _cryDisplayFS (togl, wc);
    } 
  else 
    {
      /* 
	 ANTIALIASING 
      */
      GLint viewport[4];
      GLfloat sx, sy, size, sizeX, sizeY, aspect;
      GLfloat scale, dx, dy;
      int min, max, count, i, j;
      int width  = Togl_Width (togl);
      int height = Togl_Height(togl);
      enum {
	XORG, YORG, WID, HT
      }; 

      size   = _getSize (wc);
      aspect = (GLfloat) width / (GLfloat) height;
      if ( aspect>1 ) {
	sizeX = aspect * size;
	sizeY = size;
      } else {
	sizeX = size;
	sizeY = size / aspect;
      }
      glGetIntegerv(GL_VIEWPORT, viewport);

      sx = 2.0 * sizeX / viewport[WID];
      sy = 2.0 * sizeX / viewport[WID];

      min    = -antialias.degree;
      max    = -min + 1;
      count  = -min + max;
      count *= count;

      /* uniform scaling, less than one pixel wide */
      scale = -antialias.offset / min;

      glClear(GL_ACCUM_BUFFER_BIT);

      for (j = min; j < max; j++) {
	for (i = min; i < max; i++) 
	  {
	    dx = sx * scale * i;
	    dy = sy * scale * j;
	    glMatrixMode(GL_PROJECTION);
	    glLoadIdentity();	    
	    glOrtho(-sizeX + dx, sizeX + dx, 
		    -sizeY + dy, sizeY + dy, 
		    -wc->MVf.structsize, 2*wc->MVf.structsize); /* see the crySetProj.c file !!! */

	    glMatrixMode(GL_MODELVIEW);
	    _cryDisplayFS (togl, wc);
	    glAccum(GL_ACCUM, 1.0 / (GLfloat)count);
	  }
      }
      glAccum(GL_RETURN, 1.0);
    }

  Togl_SwapBuffers(togl);
}

static void _cryDisplayFS(struct Togl *togl, NEW_WIN_CONTEXT *wc) 
{
  /*GLdouble size = (GLdouble) _getSize (wc);*/

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  
  /* INSERT HERE ROTATION MATRIX */
  glLoadIdentity();
  glTranslated(0.0, 0.0, -0.5*wc->MVf.structsize); /* see crySetProj.c */
  glMultMatrixd( wc->vec.crdvec );

  /* fog */
  xcFog (togl, wc->VPf.fog, GL_FALSE);

  glDisable(GL_BLEND); 
  glDepthMask(GL_TRUE);  

  if ( wc->VPf.nsurface )    cryRenderSurface( wc , FS_SINGLE );
  if ( wc->VPf.dispLattice ) cryRenderRecCell( wc , FS_SINGLE );
  /*--------------------*/
  glFlush();      
}


void cryRenderRecCell(  NEW_WIN_CONTEXT *wc , int fs_mode ) 
{
  MOL_SURF *m = (MOL_SURF *) wc->VPf.surfacePtr[0];

  /* so far *******************************************************/
  LoadVoronoiMaterial(VORONOI_BZ);
  /*glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_voronoi[0].ambient );*/
  glLineWidth( 5.0 );
  /****************************************************************/

  glPushMatrix();
  if ( fs_mode == FS_SINGLE )
    glTranslated( wc->MVf.o_shift[0], wc->MVf.o_shift[1], wc->MVf.o_shift[2] );

  if ( m->fs.celltype == XCR_PARAPIPEDAL ) {
    /****************************/
    /* PARAPIPEDAL LATTICE CAGE */
    /****************************/
    CellCage *cage = (CellCage*) wc->recprim_cage;

    /* BUG */
    if ( cage == (CellCage *)NULL ) return;

    if ( wc->VPf.dispLatType == CELL_WIRE ) 
      {
	crySet_Trans_Shade_Lightmodel( 0, GL_FLAT, 1, GL_FALSE );
	glColor4fv (m->fs.wirecellcolor);
	xcWireCage( *cage );
      }
    else if ( wc->VPf.dispLatType == CELL_SOLID ) 
      {
	crySet_Trans_Shade_Lightmodel( 1, GL_FLAT, 0, GL_TRUE );
	glBlendFunc( blend_cellcage.sfunc, blend_cellcage.dfunc );
	if ( m->fs.solidcellcolor[0] > -1e-8 ) {
	  glMaterialfv ( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, m->fs.solidcellcolor );
	}
	xcSolidCage( *cage );
	glDisable(GL_BLEND); 
	glDepthMask(GL_TRUE);  
      }
    else if ( wc->VPf.dispLatType == CELL_SOLID_AND_WIRE ) 
      {
	crySet_Trans_Shade_Lightmodel( 0, GL_FLAT, 1, GL_FALSE );
	glColor4fv (m->fs.wirecellcolor);
	xcWireCage( *cage );

	crySet_Trans_Shade_Lightmodel( 1, GL_FLAT, 0, GL_TRUE );
	glBlendFunc( blend_cellcage.sfunc, blend_cellcage.dfunc );
	if ( m->fs.solidcellcolor[0] > -1e-8 ) {
	  glMaterialfv ( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, m->fs.solidcellcolor );
	}
	glEnable(GL_POLYGON_OFFSET_FILL); glPolygonOffset(1.0, 1.0);
	xcSolidCage( *cage );
	glPolygonOffset(0.0, 0.0); glDisable(GL_POLYGON_OFFSET_FILL);
	glDisable(GL_BLEND); 
	glDepthMask(GL_TRUE);  
      }
  }

  else if ( m->fs.celltype == XCR_BZ ) {
    GLenum what;
    /************************/
    /* first BRILLOUIN ZONE */
    /************************/
    if ( wc->VPf.dispLatType == CELL_WIRE ) 
      {
	what = GL_LINE_LOOP;
	crySet_Trans_Shade_Lightmodel( 0, GL_FLAT, 1, GL_FALSE );
	glDisable(GL_LIGHTING);
	glColor4fv (m->fs.wirecellcolor);
	RenderBZ( what );
	glEnable(GL_LIGHTING);
      } 
    else if ( wc->VPf.dispLatType == CELL_SOLID ) 
      {
	what = GL_POLYGON;
	crySet_Trans_Shade_Lightmodel( 1, GL_FLAT, 0, GL_TRUE );
	glBlendFunc( blend_cellcage.sfunc, blend_cellcage.dfunc );
	if ( m->fs.solidcellcolor[0] > -1e-8 ) {
	  glMaterialfv ( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, m->fs.solidcellcolor );
	}
	RenderBZ( what );
	glDisable(GL_BLEND); 
	glDepthMask(GL_TRUE);  
      }
    else if ( wc->VPf.dispLatType == CELL_SOLID_AND_WIRE ) 
      {
	what = GL_LINE_LOOP;
	crySet_Trans_Shade_Lightmodel( 0, GL_FLAT, 1, GL_FALSE );
	glDisable(GL_LIGHTING);
	glColor4fv (m->fs.wirecellcolor);
	RenderBZ( what );
	glEnable(GL_LIGHTING);		  

	what = GL_POLYGON;
	crySet_Trans_Shade_Lightmodel( 1, GL_FLAT, 0, GL_TRUE );
	glBlendFunc( blend_cellcage.sfunc, blend_cellcage.dfunc );
	if ( m->fs.solidcellcolor[0] > -1e-8 ) {
	  glMaterialfv ( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, m->fs.solidcellcolor );
	}
	glEnable(GL_POLYGON_OFFSET_FILL); glPolygonOffset(1.0, 1.0);
	RenderBZ( what );
	glPolygonOffset(0.0, 0.0); glDisable(GL_POLYGON_OFFSET_FILL);
	glDisable(GL_BLEND); 
	glDepthMask(GL_TRUE);  
      } 
  }
  glPopMatrix();
  /* DEBUGGING */
  /*
    glPushMatrix();
    glTranslated( m->lowcoor[0], m->lowcoor[1], m->lowcoor[2] );

    if (1) {
    int      i, j;
    double   vec[4][4];
    CellCage cage;
    for(i=0; i<3; i++)
    for(j=0; j<3; j++)
    vec[i][j] = m->vec[i][j];
    SetUnitCellCage ( vec, &cage );
    xcWireCage(  cage );
    xcSolidCage( cage );
    }
    glPopMatrix();
  */
}

static void RenderBZ(GLenum what) {
  register int ip, iv;

  for (ip=0; ip < bz[BZ_PRIMCELL].npoly; ip++) {
    glBegin( what );
    for (iv=0; iv < bz[BZ_PRIMCELL].nvert[ip]; iv++) {
      glNormal3fv( bz[BZ_PRIMCELL].norm[ip][iv] );
      glVertex3fv( bz[BZ_PRIMCELL].poly[ip][iv] );
    }
    glEnd();
  }      
}


/*===========================================================================*/

void cryRenderSurface( NEW_WIN_CONTEXT *wc , int fs_mode ) 
{
  register int index, is, i, ix, iy, iz, j3;
  int          revertnormals = 0;
  GLint        two_side;
  GLfloat      tr_vec[3];
  MOL_SURF     *m;
  ISOSURFACE   *iso;
  GLenum       front = GL_FRONT, back  = GL_BACK;

  LoadIsoMaterial( MAT_ONELEVEL | MAT_ISOSURF );
  glLineWidth( 2.0 );
  for (is=0; is<wc->VPf.nsurface; is++) {
    if (!wc->VPf.surface[is]) continue;

    m   = (MOL_SURF *) wc->VPf.surfacePtr[is];
    iso = FindIsoSurf( m->isosurf_index );
    if (!iso) fprintf(stderr, 
		      "Isosurface # %d not found !!!\n", m->isosurf_index);

    two_side = 1.0;
    crySet_Trans_Shade_Lightmodel( m->dispt.transparent,
				   m->dispt.shademodel,
				   two_side, GL_FALSE );
    if ( m->dispt.transparent )
      glBlendFunc( blend_isosurf.sfunc, blend_isosurf.dfunc );

    /* RENDER */
    /*
      if ( m->colorscheme == MOLS_MONOCHROME ) {}
    */

    /*
      FRONTFACE
    */    
    if ( m->frontface == GL_CCW ) {
      front = GL_BACK;
      back  = GL_FRONT;
    }
    glMaterialfv(front, GL_AMBIENT_AND_DIFFUSE, m->monocolor);
    glMaterialfv(back,  GL_AMBIENT_AND_DIFFUSE, m->back_monocolor);

    /* REVERTNORMALS */
    if ( m->revertnormals && !iso->revertnormals ) {
      revertnormals      = 1;
      iso->revertnormals = 1;    
    } else if (!m->revertnormals && iso->revertnormals) {
      revertnormals      = 1;
      iso->revertnormals = 0;
    }
    if ( revertnormals ) 
      {
	for (i=0; i<iso->ntriangl; i++) {
	  for (j3=0; j3<3; j3++) {
	    index = iso->tri2verIN[i].ver_ind[j3];
	    iso->normal[index].x *= -1.0;
	    iso->normal[index].y *= -1.0;
	    iso->normal[index].z *= -1.0;
	  }
	}
      }


    glPushMatrix();
    if ( fs_mode == FS_SINGLE )
      glTranslated( wc->MVf.o_shift[0], wc->MVf.o_shift[1], wc->MVf.o_shift[2] );
    for(ix=0; ix<m->isoexpand.irepvec[0]; ix++)
      for(iy=0; iy<m->isoexpand.irepvec[1]; iy++)
	for(iz=0; iz<m->isoexpand.irepvec[2]; iz++) {
	  tr_vec[0] = 
	    (float) ix * (float) m->isoexpand.rep_vec[0][0] +
	    (float) iy * (float) m->isoexpand.rep_vec[1][0] +
	    (float) iz * (float) m->isoexpand.rep_vec[2][0];
	  tr_vec[1] = 
	    (float) ix * (float) m->isoexpand.rep_vec[0][1] +
	    (float) iy * (float) m->isoexpand.rep_vec[1][1] +
	    (float) iz * (float) m->isoexpand.rep_vec[2][1];
	  tr_vec[2] = 			  
	    (float) ix * (float) m->isoexpand.rep_vec[0][2] +
	    (float) iy * (float) m->isoexpand.rep_vec[1][2] +
	    (float) iz * (float) m->isoexpand.rep_vec[2][2];
	  glPushMatrix();
	  glTranslated( tr_vec[0], tr_vec[1], tr_vec[2] );
	  if ( m->dispt.drawstyle == ISOSURF_SOLID ) {
	    glBegin( GL_TRIANGLES );  
	    for (i=0; i<iso->ntriangl; i++) {
	      if ( !iso->triangl_status[i] ) continue;
	      for (j3=0; j3<3; j3++) {
		index = iso->tri2verIN[i].ver_ind[j3];
		/*
		  if ( m->colorscheme != MOLS_MONOCHROME )
		  glMaterialfv(GL_FRONT_AND_BACK, 
		  GL_AMBIENT_AND_DIFFUSE, 
		  iso->color[index]);
		*/
		glNormal3fv( &(iso->normal[index].x) );
		glVertex3fv( &(iso->vertex[index].x) );
	      }
	    } 
	    glEnd();  
	  } else {
	    GLenum what;
	    /* ISOSURF_WIRE && ISOSURF_DOT */
	    if ( m->dispt.drawstyle == ISOSURF_WIRE ) what=GL_LINE_LOOP;
	    else {
	      glPointSize(3.0);
	      what = GL_POINTS;
	    }
	    for (i=0; i<iso->ntriangl; i++) {
	      if ( !iso->triangl_status[i] ) continue;
	      glBegin( what );
	      for (j3=0; j3<3; j3++) {
		index = iso->tri2verIN[i].ver_ind[j3];		    
		/*
		  if ( m->colorscheme != MOLS_MONOCHROME )
		  glMaterialfv(GL_FRONT_AND_BACK, 
		  GL_AMBIENT_AND_DIFFUSE, 
		  iso->color[index]);
		*/
		glNormal3fv( &(iso->normal[index].x) );
		glVertex3fv( &(iso->vertex[index].x) );
	      }
	      glEnd(); 
	    }
	  }
	  glPopMatrix();
	}
    glPopMatrix();
  }
}


void
crySet_Trans_Shade_Lightmodel( int transp, GLenum shade, 
			       GLint two_side, GLboolean cull )
{
  GLint shade_model, two_side_model;
  GLboolean cull_face;

  /* TRANSPARENCY */
  if ( transp ) {
    glEnable( GL_BLEND );
    glDepthMask( GL_FALSE );
  } else {
    glDisable(GL_BLEND); 
    glDepthMask(GL_TRUE);  
  }

  /* SHADEMODEL */
  glGetIntegerv( GL_SHADE_MODEL, &shade_model );
  if ( shade_model != shade ) 
    glShadeModel( shade );    

  /* LIGHTMODEL */
  glGetIntegerv( GL_LIGHT_MODEL_TWO_SIDE, &two_side_model );
  if ( two_side_model != two_side ) 
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, two_side );

  /* CULL FACE */
  glGetBooleanv( GL_CULL_FACE, &cull_face );
  if ( cull_face && !cull ) 
    glDisable( GL_CULL_FACE );
  else if ( !cull_face && cull ) {
    glEnable(GL_CULL_FACE );
    glCullFace( GL_BACK );
  }
}

