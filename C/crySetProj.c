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
 * Source: $XCRYSDEN_TOPDIR/C/crySetProj.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <togl.h>
#include <math.h>
#include "struct.h"
#include "xcfunc.h"

/* crySetProj.c */
void crySetProjection(NEW_WIN_CONTEXT *wc, struct Togl *togl);

void
crySetProjection( NEW_WIN_CONTEXT *wc, struct Togl *togl )
{
  double aspect, size, size_a;
  OrthoProj o;
  int width  = Togl_Width (togl);
  int height = Togl_Height(togl);

  wc->VPf.size = (double) width;
  aspect = (double) width / (double) height;

  if (wc->VPf.stropened) {
    o.minx = wc->ss.minX / wc->tr.zoom;
    o.maxx = wc->ss.maxX / wc->tr.zoom;

    o.miny = wc->ss.minY / wc->tr.zoom;
    o.maxy = wc->ss.maxY / wc->tr.zoom;
    
    o.minz = 2.0 * wc->ss.minZ / wc->tr.zoom;
    o.maxz = 2.0 * wc->ss.maxZ / wc->tr.zoom;
    
    if ( aspect > 1 ) {
      o.minx *= aspect;
      o.maxx *= aspect;
    } else {
      wc->VPf.size = wc->VPf.height;
      o.miny /= aspect;
      o.maxy /= aspect;
    }
  } else {
    o.minx = o.miny = o.minz = -1.0;
    o.maxx = o.maxy = o.maxz = -1.0;
  }

  o.size = o.maxx;
  if ( o.maxy > o.size ) o.size = o.maxy;
  if ( o.maxz > o.size ) o.size = o.maxz;

  if (o.size < 1e-20) o.size = 1.0;
  wc->VPf.VPfactor = wc->VPf.size / (2.0 * o.size);

  wc->MVf.structsize = sqrt( wc->ss.maxX*wc->ss.maxX +
			     wc->ss.maxY*wc->ss.maxY +
			     wc->ss.maxZ*wc->ss.maxZ );

  size = wc->MVf.structsize / wc->tr.zoom;

  Togl_MakeCurrent( togl );
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if ( aspect>1 ) {
    size_a = aspect * size;
    /*glOrtho(-size_a, size_a, -size, size, -size, 2*size);*/
    glOrtho(-size_a, size_a, -size, size, -wc->MVf.structsize, 2*wc->MVf.structsize);
  } else {
    size_a = size / aspect;
    /*glOrtho(-size, size, -size_a, size_a, -size, 2*size);*/
    glOrtho(-size, size, -size_a, size_a, -wc->MVf.structsize, 2*wc->MVf.structsize);
  }
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glClearColor( wc->bg[0], wc->bg[1], wc->bg[2], wc->bg[3] ); /* T.K. */
  glViewport(0, 0, width, height);

  /*glLoadIdentity();*/
  LoadLights();
  {
    int il; 
    for ( il=0; il<GLPAR_MAXLIGHT; il++) {
      if ( light[il].isenabled ) {
	/*CalcLightPosition (il, (GLdouble)_getSize(wc));*/
	CalcLightPosition (il, size);
      }
    }
  }
}
