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
 * Source: $XCRYSDEN_TOPDIR/C/isorender.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#define XC_CPP_NO_STDIO
#include <GL/gl.h>
#include "xcGLparam.h"
#include "struct.h"
#include "isosurf.h"
#include "molsurf.h"
#include "xcfunc.h"

extern ISO_ATTRIB isoDisp, isoplaneDisp[MAX_ISOOBJECTS];
extern PLANEVERTEX ***plvertex;

void xcRenderIsosurf(int obj);
void xcRenderColorplane(int obj);
void xcRenderIsoLine2D(int obj);
static void DrawColorPlane(int obj);
static void DrawIsoLine2D(int obj);

/* --- lighting.c --- */
extern void LoadStructMaterial(void);
extern void LoadIsoMaterial(int type);
extern void LoadLightModel(int type);

/* isosurface.c */
extern ISOSURFACE *FindIsoSurf(int index);
extern void AddToIsoSurfList(ISOSURFACE *g);
extern int NewIsoSurf(void);

void
xcRenderIsosurf(int obj)
{
  register int j, ix, iy, iz, j3;
  int is_solid, index, ind[3] = {0, 1, 2};
  GLuint i;
  GLint shade_model;
  GLboolean switch_shademodel = GL_FALSE;
  GLboolean switch_lightmodel = GL_FALSE;
  GLfloat tr_vec[3];

  /* TRANSPARENCY */
  if ( isoDisp.transparent == ISOSURF_TRANSP_ON ) {
    glEnable( GL_BLEND );
    glDepthMask( GL_FALSE );
    glBlendFunc( blend_isosurf.sfunc, blend_isosurf.dfunc );
  }

  /* SHADEMODEL */
  glGetIntegerv( GL_SHADE_MODEL, &shade_model );
  if ( shade_model != isoDisp.shademodel ) { 
    switch_shademodel = GL_TRUE;
    glShadeModel( isoDisp.shademodel );
  } 

  /* LIGHTMODEL */
  if ( lightmodel.two_side[0] != lightmodel.two_side_iso[0] ) {
    switch_lightmodel = GL_TRUE;
    glLightModelfv( GL_LIGHT_MODEL_TWO_SIDE, lightmodel.two_side_iso );
  }

  /* isosurface render-style */
  if ( isoDisp.drawstyle == ISOSURF_SOLID ) {
    is_solid = MAT_SOLID;
  } else {
    is_solid = 0;
  }
  
  /****************************************************/
  /* load isosurf's material properties               */
  for (j=0; j < isodata.nlevel; j++) {
    /* FRONTFACE */    
    /* glGetIntegerv ( GL_FRONT_FACE, &front_face );
       if ( frontface_isosurf[j] != front_face ) {
       switch_frontface = GL_TRUE;
       glFrontFace( frontface_isosurf );
       }
    */
    
    if ( isodata.nlevel == 1 ) {
      LoadIsoMaterial( MAT_ONELEVEL | MAT_ISOSURF );
    } else if ( isodata.nlevel == 2 && j == 0 ) {
      LoadIsoMaterial( MAT_POSLEVEL | MAT_ISOSURF );
    } else if ( isodata.nlevel == 2 && j == 1 ) {
      LoadIsoMaterial( MAT_NEGLEVEL | MAT_ISOSURF | is_solid );
    }

    /* t.k. */
    if ( isodata.cell_orientation == XC_RIGHT ) {
      ind[0] = 2;
      ind[1] = 1;
      ind[2] = 0;
    }
       
    /* RENDER */
    glPushMatrix();
      glTranslated( -mx + isodata.points[obj][0][0], 
		    -my + isodata.points[obj][0][1], 
		    -mz + isodata.points[obj][0][2] );
      for(ix=0; ix<isoexpand[obj].irepvec[0]; ix++)
	for(iy=0; iy<isoexpand[obj].irepvec[1]; iy++)
	  for(iz=0; iz<isoexpand[obj].irepvec[2]; iz++) {
	    tr_vec[0] = 
	      (float) ix * (float) isoexpand[obj].rep_vec[0][0] +
	      (float) iy * (float) isoexpand[obj].rep_vec[1][0] +
	      (float) iz * (float) isoexpand[obj].rep_vec[2][0];
	    tr_vec[1] = 
	      (float) ix * (float) isoexpand[obj].rep_vec[0][1] +
	      (float) iy * (float) isoexpand[obj].rep_vec[1][1] +
	      (float) iz * (float) isoexpand[obj].rep_vec[2][1];
	    tr_vec[2] = 
	      (float) ix * (float) isoexpand[obj].rep_vec[0][2] +
	      (float) iy * (float) isoexpand[obj].rep_vec[1][2] +
	      (float) iz * (float) isoexpand[obj].rep_vec[2][2];
	    glPushMatrix();
	      glTranslated( tr_vec[0], tr_vec[1], tr_vec[2] );
	      if ( isoDisp.drawstyle == ISOSURF_SOLID ) {
		glBegin( GL_TRIANGLES );  
		  for (i=0; i<surf.numT[j]; i++) {
		    for (j3=0; j3<3; j3++) {
		      index = tri2verIN[j][i].ver_ind[ind[j3]];		    
		      glNormal3fv( &(normal[j][index].x) );
		      glVertex3fv( &(vertex[j][index].x) );
		    }
		  } 
		glEnd();  
	      } else {
		GLenum what;
		/* ISOSURF_WIRE && ISOSURF_DOT */
		if ( isoDisp.drawstyle == ISOSURF_WIRE ) {
		  glLineWidth (VPf.WF3Dlinewidth);
		  what=GL_LINE_LOOP;
		}
		else {
		  glPointSize (VPf.PLradius);
		  what = GL_POINTS;
		}
		for (i=0; i<surf.numT[j]; i++) {
		  glBegin( what );
		    for (j3=0; j3<3; j3++) {
		      index = tri2verIN[j][i].ver_ind[j3];		    
		      glNormal3fv( &(normal[j][index].x) );
		      glVertex3fv( &(vertex[j][index].x) );
		    }
		  glEnd(); 
		}
	      }
	    glPopMatrix();
	  }
    glPopMatrix();

    /* GO to PREVIOUS FRONTFACE */
    /* if ( switch_frontface ) glFrontFace( front_face ); */
  }
  
  /* DISABLE TRANSPARENCY */
   if ( isoDisp.transparent == ISOSURF_TRANSP_ON ) {      
     glDisable(GL_BLEND); 
     glDepthMask(GL_TRUE);  
   }

   /* GO BACK TO PREVIOUS SHADE_MODEL */
   if ( switch_shademodel ) glShadeModel( shade_model ); 
   
   /* GO BACK TO PREVIOUS LIGHT_MODEL */
   if ( switch_lightmodel ) 
     glLightModelfv( GL_LIGHT_MODEL_TWO_SIDE, lightmodel.two_side );

   /* LOAD THE STRUCT MATERIALS BACK; 
    * this has to be changed in future; 
    * each render function will have to take care of it-self
     */
   LoadStructMaterial();  
}



void
xcRenderColorplane(int obj)
{
  int ix, iy, iz;  
  GLfloat tr_vec[3];
  GLint shade_model;
  GLboolean light_model;
  GLboolean switch_shademodel = GL_FALSE, switch_lightmodel = GL_FALSE;
  float specular[4] = { 0.0, 0.0, 0.0, 1.0 };
  
   /*
     #ifdef GL_EXT_vertex_array
     if (use_vertex_arrays) {
     glDrawArraysEXT( GL_TRIANGLES, 0, numverts );
     }
     else {
     #endif
     */
  /* Load Materials */
  LoadIsoMaterial( MAT_COLORPLANE );

  if ( isoplaneDisp[obj].lighting ) {
    glEnable( GL_LIGHTING );
  } else {
    glDisable( GL_LIGHTING );
  }

  /* TRANSPARENCY */
  if ( isoplaneDisp[obj].transparent == ISOSURF_TRANSP_ON ) {
    glEnable( GL_BLEND );
    glDepthMask( GL_FALSE );
    glBlendFunc( blend_colorplane.sfunc, blend_colorplane.dfunc );
  }

  /* SHADEMODEL */
  glGetIntegerv( GL_SHADE_MODEL, &shade_model );
  if ( shade_model != GL_SMOOTH ) { 
    switch_shademodel = GL_TRUE;
    glShadeModel( GL_SMOOTH );
  }

  /* LIGHTMODEL */
  glGetBooleanv( GL_LIGHT_MODEL_TWO_SIDE, &light_model );
  if ( light_model != GL_TRUE ) {
    switch_lightmodel = GL_TRUE;
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );
  }

  /* Polygon-Offset: 
   *
   *  Put POLYGONs behind the ISOLINEs so that isolines will be nicely 
   * visible !!!
   */
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.0, 1.);

  glPushMatrix();
    glTranslated( -mx + isodata.points[obj][0][0], 
		  -my + isodata.points[obj][0][1],
		  -mz + isodata.points[obj][0][2] );
    glEnable( GL_COLOR_MATERIAL );
    /* --- this is temporate --- */
    glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specular );
    /* ---   end temporate   --- */    
    glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
    
    for (ix=0; ix<isoexpand[obj].irepvec[0]; ix++)
      for (iy=0; iy<isoexpand[obj].irepvec[1]; iy++)
	for (iz=0; iz<isoexpand[obj].irepvec[2]; iz++) {
	  tr_vec[0] = 
	    (float) ix * (float) isoexpand[obj].rep_vec[0][0] +
	    (float) iy * (float) isoexpand[obj].rep_vec[1][0] +
	    (float) iz * (float) isoexpand[obj].rep_vec[2][0];
	  tr_vec[1] = 
	    (float) ix * (float) isoexpand[obj].rep_vec[0][1] +
	    (float) iy * (float) isoexpand[obj].rep_vec[1][1] +
	    (float) iz * (float) isoexpand[obj].rep_vec[2][1];
	  tr_vec[2] = 
	    (float) ix * (float) isoexpand[obj].rep_vec[0][2] +
	    (float) iy * (float) isoexpand[obj].rep_vec[1][2] +
	    (float) iz * (float) isoexpand[obj].rep_vec[2][2];
	  glPushMatrix();
	    glTranslated( tr_vec[0], tr_vec[1], tr_vec[2] );
	    DrawColorPlane( obj );
	  glPopMatrix();
	}
    glDisable( GL_COLOR_MATERIAL );
   glPopMatrix();

   /* disable polygon-offset */
   glPolygonOffset(0.0, 0.0);
   glDisable(GL_POLYGON_OFFSET_FILL);

   /* DISABLE TRANSPARENCY */
   if ( isoplaneDisp[obj].transparent == ISOSURF_TRANSP_ON ) {      
     glDisable(GL_BLEND); 
     glDepthMask(GL_TRUE);  
   }

   /* GO BACK TO PREVIOUS SHADE MODEL */
   if ( switch_shademodel ) glShadeModel( shade_model ); 

   /* GO BACK TO PREVIOUS LIGHT_MODEL */
   if ( switch_lightmodel ) 
     glLightModelf( GL_LIGHT_MODEL_TWO_SIDE, light_model );

   glEnable( GL_LIGHTING );

   /* LOAD THE STRUCT MATERUALS BACK; 
    * this has to be changed in future; 
    * each render function will have to take care of it-self
    */
   LoadStructMaterial();

   /*
     #ifdef GL_EXT_vertex_array
     }
     #endif
     */
}


static void
DrawColorPlane(int obj)
{
  register int n, m, j, k;
  float midC[4], midP[3];

  
  n = newgrd.nx - 1;
  m = newgrd.ny - 1;
  if ( obj == ISOOBJ_PLANE2 ) {
    n = newgrd.nz - 1;
    m = newgrd.nx - 1;
  }
  if ( obj == ISOOBJ_PLANE3 ) {
    n = newgrd.ny - 1;
    m = newgrd.nz - 1;
  }
  
  for (j=0; j<n; j++) 
    for (k=0; k<m; k++) 
      {	
	/* midpoint */
	midC[0] = 0.25*(plvertex[obj][j][k].col[0] + plvertex[obj][j+1][k].col[0] + plvertex[obj][j+1][k+1].col[0] + plvertex[obj][j][k+1].col[0]);
	midC[1] = 0.25*(plvertex[obj][j][k].col[1] + plvertex[obj][j+1][k].col[1] + plvertex[obj][j+1][k+1].col[1] + plvertex[obj][j][k+1].col[1]);
	midC[2] = 0.25*(plvertex[obj][j][k].col[2] + plvertex[obj][j+1][k].col[2] + plvertex[obj][j+1][k+1].col[2] + plvertex[obj][j][k+1].col[2]);
	midC[3] = 0.25*(plvertex[obj][j][k].col[3] + plvertex[obj][j+1][k].col[3] + plvertex[obj][j+1][k+1].col[3] + plvertex[obj][j][k+1].col[3]);
      
	midP[0] = 0.5*(plvertex[obj][j][k].p[0] + plvertex[obj][j+1][k+1].p[0]);
	midP[1] = 0.5*(plvertex[obj][j][k].p[1] + plvertex[obj][j+1][k+1].p[1]);
	midP[2] = 0.5*(plvertex[obj][j][k].p[2] + plvertex[obj][j+1][k+1].p[2]);
		  
	glBegin( GL_TRIANGLE_FAN );
	glNormal3fv( &(isodata.colnml[obj].x) );

	/* order is: midpoint -> (j,k) -> (j,k+1) -> (j+1,k+1) -> (j+1,k) -> (j,k) */
	glColor4fv( midC );
	glVertex3fv( midP );
	     
	glColor4fv ( plvertex[obj][j  ][k  ].col );
	glVertex3fv( plvertex[obj][j  ][k  ].p );
	
	glColor4fv ( plvertex[obj][j  ][k+1].col );
	glVertex3fv( plvertex[obj][j  ][k+1].p );
	
	glColor4fv ( plvertex[obj][j+1][k+1].col );
	glVertex3fv( plvertex[obj][j+1][k+1].p );

	glColor4fv ( plvertex[obj][j+1][k  ].col );
	glVertex3fv( plvertex[obj][j+1][k  ].p );

	glColor4fv ( plvertex[obj][j  ][k  ].col );
	glVertex3fv( plvertex[obj][j  ][k  ].p );
	glEnd();
      }
}



void
xcRenderIsoLine2D(int obj)
{
  int ix, iy, iz;  
  GLenum depthf;
  GLfloat tr_vec[3];
  GLboolean switch_depthf = GL_FALSE;
  float specular[4] = { 0.0, 0.0, 0.0, 1.0 };

  LoadIsoMaterial( MAT_COLORPLANE );

  /* DISABLE LIGHTING */
  glDisable( GL_LIGHTING );

  /* DEPTH FUNC */
  glGetIntegerv(GL_DEPTH_FUNC, (GLint *) &depthf);
  if ( depthf != GL_LEQUAL ) {
    switch_depthf = GL_TRUE;
    glDepthFunc(GL_LEQUAL);
  }

  /* SHADEMODEL */
  /*
    glGetIntegerv( GL_SHADE_MODEL, &shade_model );
    if ( shade_model != GL_FLAT ) { 
    switch_shademodel = GL_TRUE;
    glShadeModel( GL_FLAT );
  }
  */

  glPushMatrix();
    glTranslated( -mx + isodata.points[obj][0][0], 
		  -my + isodata.points[obj][0][1],
		  -mz + isodata.points[obj][0][2] );
    glEnable( GL_COLOR_MATERIAL );
    /* --- this is temporate --- */
    glMaterialfv( GL_FRONT, GL_SPECULAR, specular );
    /* ---   end temporate   --- */    
    glColorMaterial( GL_FRONT, GL_AMBIENT_AND_DIFFUSE );
    
    for (ix=0; ix<isoexpand[obj].irepvec[0]; ix++)
      for (iy=0; iy<isoexpand[obj].irepvec[1]; iy++)
	for (iz=0; iz<isoexpand[obj].irepvec[2]; iz++) {
	  tr_vec[0] = 
	    (float) ix * (float) isoexpand[obj].rep_vec[0][0] +
	    (float) iy * (float) isoexpand[obj].rep_vec[1][0] +
	    (float) iz * (float) isoexpand[obj].rep_vec[2][0];
	  tr_vec[1] = 
	    (float) ix * (float) isoexpand[obj].rep_vec[0][1] +
	    (float) iy * (float) isoexpand[obj].rep_vec[1][1] +
	    (float) iz * (float) isoexpand[obj].rep_vec[2][1];
	  tr_vec[2] = 
	    (float) ix * (float) isoexpand[obj].rep_vec[0][2] +
	    (float) iy * (float) isoexpand[obj].rep_vec[1][2] +
	    (float) iz * (float) isoexpand[obj].rep_vec[2][2];
	  glPushMatrix();
	    glTranslated( tr_vec[0], tr_vec[1], tr_vec[2] );
	    DrawIsoLine2D( obj );
	  glPopMatrix();
	}
    glDisable( GL_COLOR_MATERIAL );
   glPopMatrix();

   /* GO BACK TO PREVIOUS DEPTHFUNC */
   if ( switch_depthf ) glDepthFunc( depthf ); 

   /* GO BACK TO PREVIOUS SHADE MODEL */
   /*
     if ( switch_shademodel ) glShadeModel( shade_model ); 
   */
   glEnable( GL_LIGHTING );

   /* LOAD THE STRUCT MATERUALS BACK; 
    * this has to be changed in future; 
    * each render function will have to take care of it-self
    */
   LoadStructMaterial();
}


static void
DrawIsoLine2D(int obj)
{
  register int nlevel, nseg, j, k;
  float *color; 
  unsigned short dash;
  int dashf;
  int stipple_enabled = 0;

  nlevel = isoline2D[obj].nlevel;

  color = &isoline2D[obj].color[0][0];
  dash  = isoline2D[obj].dash[0];
  dashf = isoline2D[obj].dashfactor[0];
  if ( dash != LINESTIPPLE_SOLID ) {
    glLineStipple( dashf, dash );
    glEnable( GL_LINE_STIPPLE );
    stipple_enabled = 1;
  }

  glLineWidth( isoline2D[obj].linewidth );
  glColor4fv( color );
  for (j=0; j<nlevel; j++) {
    nseg = isoline2D[obj].iseg[j];
    /* color */
    if ( !IsEqualf( isoline2D[obj].color[j][0], color[0], 1.e-6 ) || \
	 !IsEqualf( isoline2D[obj].color[j][1], color[1], 1.e-6 ) || \
	 !IsEqualf( isoline2D[obj].color[j][2], color[2], 1.e-6 ) ) {
      color = &isoline2D[obj].color[j][0];
      glColor4fv( isoline2D[obj].color[j] );
    }
    /* stipple */
    if ( isoline2D[obj].dash[j] != dash || 
	 isoline2D[obj].dashfactor[j] != dashf ) {
      dash  = isoline2D[obj].dash[j];
      dashf = isoline2D[obj].dashfactor[j];
      glLineStipple( dashf, dash );
    }
    if ( stipple_enabled && dash == LINESTIPPLE_SOLID ) {
      glDisable( GL_LINE_STIPPLE );
      stipple_enabled = 0;
    }
    if ( !stipple_enabled && dash != LINESTIPPLE_SOLID ) {
      glEnable( GL_LINE_STIPPLE );
      stipple_enabled = 1;
    }      
    glBegin( GL_LINES );
      glNormal3fv( &(isodata.colnml[obj].x) ); /* front face */
      for (k=0; k<nseg; k++) {
	glVertex3fv( &isoline2D[obj].segment[j][k].p1.x );
	glVertex3fv( &isoline2D[obj].segment[j][k].p2.x );
      }
    glEnd();
  }
  glDisable( GL_LINE_STIPPLE );
}


/*===========================================================================*/

void xcRenderSurface(void) 
{
  register int is, i, ix, iy, iz, nx, ny, nz, j3;
  int          index;
  GLint        shade_model;
  GLboolean    switch_shademodel = GL_FALSE;
  GLboolean    switch_lightmodel = GL_FALSE;
  GLfloat      tr_vec[3];
  MOL_SURF     *m;
  ISOSURFACE   *iso;

  LoadBlendfunc_And_Frontface();

  for (is=0; is<VPf.nsurface; is++) {
    if (!VPf.surface[is]) continue;
    
    m   = (MOL_SURF *) VPf.surfacePtr[is];
    iso = FindIsoSurf( m->isosurf_index );
    if (!iso) {
      fprintf(stderr, "Isosurface # %d not found !!!\n", m->isosurf_index);
      continue;
    }
    
    /* TRANSPARENCY */
    if ( m->dispt.transparent == ISOSURF_TRANSP_ON ) {
      glEnable( GL_BLEND );
      glDepthMask( GL_FALSE );
      glBlendFunc( blend_isosurf.sfunc, blend_isosurf.dfunc );
    }

    /* SHADEMODEL */
    glGetIntegerv( GL_SHADE_MODEL, &shade_model );
    if ( shade_model != m->dispt.shademodel ) { 
      switch_shademodel = GL_TRUE;
      glShadeModel( m->dispt.shademodel );
    } 

    /* LIGHTMODEL */
    /* for testing(temporarily): */
    m->lightm.two_side_iso[0] = 1.0;
    if ( lightmodel.two_side[0] != m->lightm.two_side_iso[0] ) {
      switch_lightmodel = GL_TRUE;
      glLightModelfv( GL_LIGHT_MODEL_TWO_SIDE, m->lightm.two_side_iso );
    }    

    LoadIsoMaterial( MAT_ONELEVEL | MAT_ISOSURF );

    /* RENDER */
    glPushMatrix(); /* PUSH-1 */
    glTranslated( -mx, -my, -mz );
    /*
      for(ix=0; ix<m->isoexpand.irepvec[0]; ix++)
	for(iy=0; iy<m->isoexpand.irepvec[1]; iy++)
	  for(iz=0; iz<m->isoexpand.irepvec[2]; iz++) {
    */
    /* for now: */
    /*----------*/
      if ( m->colorscheme == MOLS_MONOCHROME ) {
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, m->monocolor);
	glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE,  m->monocolor);
      }

      nx = xcr.nunit[0];
      ny = xcr.nunit[1];
      nz = xcr.nunit[2];
      if ( xcr.dim < 3 ) nz = 1;
      if ( xcr.dim < 2 ) ny = 1;
      if ( xcr.dim < 1 ) nx = 1;

      for(ix=0; ix<nx; ix++)
	for(iy=0; iy<ny; iy++)
	  for(iz=0; iz<nz; iz++) {
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
	    glPushMatrix(); /* PUSH-2 */
	      glTranslated( tr_vec[0], tr_vec[1], tr_vec[2] );
	      if ( m->dispt.drawstyle == ISOSURF_SOLID ) {

		glBegin( GL_TRIANGLES );  
		  for (i=0; i<iso->ntriangl; i++) {
		    for (j3=0; j3<3; j3++) 
		      {
			index = iso->tri2verIN[i].ver_ind[j3];
			if ( m->colorscheme != MOLS_MONOCHROME )
			  {
			    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,  iso->color[index]);
			  }
			glNormal3fv( &(iso->normal[index].x) );
			glVertex3fv( &(iso->vertex[index].x) );
		      }
		  } 
		glEnd();  
	      } else {
		GLenum what;
		/* ISOSURF_WIRE && ISOSURF_DOT */
		if ( m->dispt.drawstyle == ISOSURF_WIRE ) {
		  glLineWidth (VPf.WF3Dlinewidth);
		  what=GL_LINE_LOOP;
		} else {
		  glPointSize (VPf.PLradius);
		  what = GL_POINTS;
		}
		for (i=0; i<iso->ntriangl; i++) {
		  glBegin( what );
		    for (j3=0; j3<3; j3++) {
		      index = iso->tri2verIN[i].ver_ind[j3];		    
		      if ( m->colorscheme != MOLS_MONOCHROME )
			glMaterialfv(GL_FRONT_AND_BACK, 
				     GL_AMBIENT_AND_DIFFUSE, 
				     iso->color[index]);
		      glNormal3fv( &(iso->normal[index].x) );
		      glVertex3fv( &(iso->vertex[index].x) );
		    }
		  glEnd(); 
		}
	      }
	    glPopMatrix(); /* POP-2 */
	  }
    glPopMatrix(); /* POP-1 */
  
    /* GO to PREVIOUS FRONTFACE */
    /* if ( switch_frontface ) glFrontFace( front_face ); */
    
    /* DISABLE TRANSPARENCY */
    if ( m->dispt.transparent == ISOSURF_TRANSP_ON ) {      
      glDisable(GL_BLEND); 
      glDepthMask(GL_TRUE);  
    }
    
    /* GO BACK TO PREVIOUS SHADE_MODEL */
    if ( switch_shademodel ) glShadeModel( shade_model ); 
    
    /* GO BACK TO PREVIOUS LIGHT_MODEL */
    if ( switch_lightmodel ) 
      glLightModelfv( GL_LIGHT_MODEL_TWO_SIDE, lightmodel.two_side );
  }
   /* LOAD THE STRUCT MATERIALS BACK; 
    * this has to be changed in future; 
    * each render function will have to take care of it-self
     */
   LoadStructMaterial();  
}
