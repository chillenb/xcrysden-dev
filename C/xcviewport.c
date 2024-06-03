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
 * Source: $XCRYSDEN_TOPDIR/C/xcviewport.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <stdio.h>
#include <math.h>
#include <GL/glu.h>
#include "struct.h"

extern OrthoProj ort;
extern RasterFontSize rf;
PERSPECTIVE persp;

/* --- functions prototypes --- */
void xcViewPort(void);
void MaybeClipAndMakeProjection(void);
void Screen2Model_Coords(int xs, int ys, float *xm, float *ym);
/*
void xcPerspective(PERSPECTIVE *per, GLdouble fovy, 
		   GLdouble near_factor, GLdouble size);
*/

/* this proc make ViewPort Transformation */
void
xcViewPort(void)
{
  double size, scrF, xorig, yorig, x, y;

  /*
   * srcF determine how big looks structure in the screen; TR.ZOOM
   * taken into account
   */
  
  /* scrF = tr.zoom / VPf.scrf; Vpf.scrf is no more used */
  scrF = tr.zoom;
  
  xorig = (double) VPf.width / 2.0  + tr.xtransl;
  yorig = (double) VPf.height / 2.0 + tr.ytransl;

  if ( VPf.width > VPf.height) {
    size = (double) VPf.height * scrF;
  } else {
    size = (double) VPf.width * scrF;
  }
  VPf.x = xorig - size / 2.0;    
  VPf.y = yorig - size / 2.0;  
  VPf.size = size;


  /* some reshapes of window requires new cliping of projection */

  MaybeClipAndMakeProjection();


  x = VPf.x;
  y = VPf.y;
  VPf.sizex = VPf.sizey = VPf.size;
  if (VPf.x < 0) {
    VPf.x     = 0;
    VPf.sizex += x;
  }
  if (VPf.y < 0) {
    VPf.y      = 0;
    VPf.sizey += y;
  }
  if (VPf.x + VPf.sizex > VPf.width)  VPf.sizex += VPf.width - (VPf.sizex + VPf.x);
  if (VPf.y + VPf.sizey > VPf.height) VPf.sizey += VPf.height - (VPf.sizey + VPf.y);

  /*
  fprintf(stderr,"VIEWPORT::");
  fprintf(stderr,"Vpf: x=%f, y=%f; xorig=%f, yorig=%f, xtransl=%f, ytransl=%f\n",
	 VPf.x, VPf.y, xorig, yorig, tr.xtransl, tr.ytransl);
  fprintf(stderr,"viewPort= %d x %d + %d x %d\n", 
	 (int) VPf.x, (int) VPf.y, (int) VPf.sizex, (int) VPf.sizey);
  fprintf(stderr,"ortho=    %f x %f + %f x %f\n\n", ort.minx, ort.miny, ort.maxx, ort.maxy);
  */

  glViewport( (GLint) VPf.x, (GLint) VPf.y, 
	      (GLsizei) VPf.sizex, (GLsizei) VPf.sizey);

  /* 
   * VPf.VPfactor is ration between structuresize in screen units and 
   * structuresize in A; TR.ZOOM is taken into account !!!
   */

  VPf.VPfactor = VPf.size / (2.0 * ort.size);

  /* how large is font in Angstrom units (we want negative values):: */

  rf.w2 = (float) -rf.wid / (2.0 * VPf.VPfactor);
  rf.h2 = (float) -rf.height / (2.0 * VPf.VPfactor);

  /* 
   * SCREEN UNITS - ANGSTROMS FACTOR !!!!
   *
   * for each direction (X and Y) one factor -- just in case, if X & Y
   * glOrtho differ
   */

  if ( VPf.projmade ) {
    VPf.scrAnX = VPf.sizex / (ort.maxx - ort.minx); /* ort.minx is negativ */
    VPf.scrAnY = VPf.sizey / (ort.maxy - ort.miny); 
  }

  /* 
     printf("viewPort= %d x %d +%d + %d\n", VPf.x, VPf.y, VPf.size, VPf.size);
  */
}

void
MaybeClipAndMakeProjection(void)
{
  double cr = 0.0;
  double size = (double) VPf.size;
/*    GLdouble fovy, aspect, near, far; */

  /**************************************
                   Y-up
               +------------+
               |    ort.size|
       X-left  |      +-----| X-right
               |            |
               +------------+
                   Y-down

  ****************************************/

  ort.minx = ort.miny = -ort.size;
  ort.maxx = ort.maxy =  ort.size;

  ort.minz = -1.5 * ort.size;
  ort.maxz =  1.5 * ort.size;

  /* clip X-left */
  if (VPf.x < 0.0 && size > MINTOL) {
    cr = ABS(VPf.x) / size;
    ort.minx = -(1.0 - 2.0*cr) * ort.size; 
  }

  /* clip X-right */
  if (VPf.size + VPf.x > VPf.width && size > MINTOL) {
    cr = ABS(VPf.size+VPf.x - VPf.width) / size;
    ort.maxx = (1.0 - 2.0*cr) * ort.size;
  }
  
  /* clip Y-up */
  if (VPf.size + VPf.y > VPf.height && size > MINTOL) {
    cr = ABS(VPf.size+VPf.y - VPf.height) / size;
    ort.maxy = (1.0 - 2.0*cr) * ort.size;
  }

  /* clip Y-down */
  if (VPf.y < 0.0 && size > MINTOL) {
    cr = ABS(VPf.y) / size;
    ort.miny = -(1.0 - 2.0*cr) * ort.size; 
  }

  /*fprintf(stderr, "DEBUG> ort: minZ=%f, maxZ=%f\n", ort.minz, ort.maxz);*/

  /* calculate perspective from ort.XXXX */
/*    fovy   = 10.0; */
/*    aspect = (ort.maxx - ort.minx) / (ort.maxy - ort.miny); */
/*    near   = (ort.maxx - ort.minx) / tan( fovy ); */
/*    far    = near + ort.maxz - ort.minz; */
  /*      gluPerspective( fovy, aspect, near, far); */

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if ( ! VPf.perspective ) {
    glOrtho(ort.minx, ort.maxx, 
	    ort.miny, ort.maxy, 
	    ort.minz+ort.size, ort.maxz+ort.size);
  } else {
    /*
      xcPerspective (&persp, VPf.perspective_fovy, 
      VPf.perspective_fovy, ort.size);
      gluPerspective (persp.fovy, 1.0 / persp.aspect, persp.near, persp.far);
    */
        
    /*
      persp.near = VPf.perspective_fovy * ort.size;
      persp.far  = persp.near + ort.size;
    */
    double _size = VPf.perspective_size * ort.size;
    persp.near   = VPf.perspective_fovy * _size;
    persp.far    = persp.near + VPf.perspective_far * _size;
    
    /* persp.shiftZ is translated in xcDisplay3D */

    persp.shiftZ = -(persp.near + _size);
      
    glFrustum (ort.minx, ort.maxx, ort.miny, ort.maxy, persp.near, persp.far);
  }

  glMatrixMode(GL_MODELVIEW);
  glClearColor( bg[0], bg[1], bg[2], bg[3] ); /* T.K. */
  glLoadIdentity();
}


void
Screen2Model_Coords(int xs, int ys, float *xm, float *ym)
{
  float xorig, yorig;
  
  xorig = (float) VPf.width / 2.0  + tr.xtransl;
  yorig = (float) VPf.height / 2.0 - tr.ytransl;
  
  /*
   *xm = ((float) xs - xorig) / VPf.VPfactor;
   *ym = (yorig - (float) ys) / VPf.VPfactor;
   */
  
  *xm = ((float) xs - xorig) / VPf.scrAnX;
  *ym = (yorig - (float) ys) / VPf.scrAnY;
}

/*
  GLdouble *zm; 
  GLdouble modelMatrix[16];
  GLdouble projMatrix[16];
  GLdouble viewport[4];

  glGetDoublev( GL_MODELVIEW_MATRIX, modelMatrix );
  glGetDoublev( GL_PROJECTION_MATRIX, projMatrix );
  glGetIntegerv( GL_VIEWPORT, viewport);
  gluUnProject( (GLdouble) xs, (GLdouble) ys, 0.0,
		modelMatrix, projMatrix, viewport,
		xm, ym, zm );

*/
/*
void xcPerspective(PERSPECTIVE *per, GLdouble fovy, 
		   GLdouble near_factor, GLdouble size) {
  GLdouble near, fi = 0.5 * fovy / RAD2DEG;
  
  per->fovy   = fovy;
  per->aspect = (GLdouble) (VPf.sizex - VPf.x)/(GLdouble)(VPf.sizey - VPf.y);
  near        = 0.75*size * cos(fi) / sin(fi) / tr.zoom;
  per->near   = near_factor*near;
  per->far    = near + 1.5*size;
  per->shiftZ = -(near + 0.75*size) / tr.zoom;
}
*/
