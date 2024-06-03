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
 * Source: $XCRYSDEN_TOPDIR/C/fog.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#define XC_CPP_NO_STDIO
#include <togl.h>
#include "struct.h"
#include "xcfunc.h"

extern OrthoProj   ort; 
extern PERSPECTIVE persp;
extern XCfog       fog;
extern XCantialias antialias;
extern struct Togl *mesa_togl;

/* fog.c */
void xcFog(struct Togl *togl, GLboolean make_fog, GLboolean perspective);
void xcAntiAlias(GLboolean make_antialias);

void 
xcFog(struct Togl *togl, GLboolean make_fog, GLboolean perspective) {
  NEW_WIN_CONTEXT *wc;
  GLclampf        *background;
  
  if ( togl == mesa_togl ) {
    background = bg;
  } else {
    wc = FindWinContextByTogl( togl );
    background = wc->bg;
  }

  if (make_fog)
    {
      GLint fog_mode = GL_LINEAR;

      switch (fog.mode) 
	{
	case XC_FOG_EXP:
	  fog_mode = GL_EXP;
	  break;
	case XC_FOG_EXP2:
	  fog_mode = GL_EXP2;
	  break;
	}
      
      glEnable (GL_FOG);
      glFogi   (GL_FOG_MODE,    fog_mode);
      glFogf   (GL_FOG_DENSITY, fog.density);

      if (fog.colormode == XC_FOG_BGCOLOR) 
	glFogfv (GL_FOG_COLOR, background);
      else
	glFogfv (GL_FOG_COLOR, fog.color);
      
      if (!perspective) 
	{
	  glFogf   (GL_FOG_START,  fog.ort_start_f * ort.size);
	  glFogf   (GL_FOG_END,    fog.ort_end_f   * ort.size);
	}
      else
	{
	  glFogf   (GL_FOG_START, -fog.persp_f1 * persp.shiftZ);
	  glFogf   (GL_FOG_END,   -fog.persp_f1 * persp.shiftZ + fog.persp_f2 * ort.size);	    
	}    
    }
  else
    {
      glDisable (GL_FOG);
    }
}
      
  

void xcAntiAlias(GLboolean make_antialias) {
  /* this function is used for Lighting-Off anti-aliasing */

  if (make_antialias)
    {
      glEnable    (GL_LINE_SMOOTH);
      glEnable    (GL_BLEND);
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glHint      (GL_LINE_SMOOTH_HINT, GL_NICEST);
      glEnable    (GL_POINT_SMOOTH);
      /*  
	  glEnable    (GL_POLYGON_SMOOTH);
	  glHint      (GL_POLYGON_SMOOTH_HINT, GL_FASTEST);
      */      
    }
  else
    {
      glDisable (GL_BLEND);
      glDisable (GL_LINE_SMOOTH);
      glDisable (GL_POINT_SMOOTH);
      /* glDisable (GL_POLYGON_SMOOTH); */
    }
}

