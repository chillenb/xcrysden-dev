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
* Source: $XCRYSDEN_TOPDIR/C/cryDispFuncMultiFS.c
* ------                                                                    *
* Copyright (c) 1996-2004 by Anton Kokalj                                   *
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

extern void cryRenderSurface( NEW_WIN_CONTEXT *wc , int fs_mode);
extern void cryRenderRecCell( NEW_WIN_CONTEXT *wc , int fs_mode);

void cryDisplayMultiFS(struct Togl *togl);
static GLfloat _getSize (NEW_WIN_CONTEXT *wc);
static void _cryDisplayMultiFS(struct Togl *togl, NEW_WIN_CONTEXT *wc);

/*
  assign
  cry_dispfuncmultiFS <togl> ...
*/
int 
CRY_DispFuncMultiFSCmd(ClientData clientData,Tcl_Interp *interp,
		       int argc, char *argv[])
{
  Togl *togl;
  NEW_WIN_CONTEXT *wc;

  if ( argc < 2 || argc%2 ) {
    Tcl_SetResult(interp, "Usage: cry_dispfuncmultiFS <togl> ?-togllist? ?-antialias 0|1?", TCL_STATIC);
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

  wc->xcDisplay = cryDisplayMultiFS;

  wc->ss.minX = 0.0;
  wc->ss.maxX = 0.0;

  wc->ss.minY = 0.0;
  wc->ss.maxY = 0.0;

  wc->ss.minZ = 0.0;
  wc->ss.maxZ = 0.0;

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


  if (wc->fermiContext.toglVector)
    {
      free(wc->fermiContext.toglVector);
      wc->fermiContext.toglVector = NULL;
    }

  if (argc > 2) {
    int i, it, anti_alias;
    for (i=2; i<argc; i+=2)
      {
	/*
	  -TOGLLIST
	*/
	if ( strncmp(argv[i], "-togll", 6) == 0 ) {
	  int argcList;
	  const char **argvList;
	  struct Togl **toglVector;
	  Tcl_SplitList(interp, argv[i+1], &argcList, &argvList);
	  toglVector = (struct Togl**) malloc ( sizeof(struct Togl*) * argcList );
	  for (it=0; it<argcList; it++)
	    {
	      NEW_WIN_CONTEXT *ith_wc;
	      if ( Togl_GetToglFromName(interp, argvList[it], &(toglVector[it])) == TCL_ERROR ) {
		char rss[1024];
		snprintf(rss, sizeof(rss), 
			 "couldn't find %s togl widget", argvList[it]);
		Tcl_SetResult(interp, rss, TCL_VOLATILE);
		return TCL_ERROR;
	      }

	      ith_wc = FindWinContextByTogl( toglVector[it] );

	      if (ith_wc->VPf.stropened) 
		{
		  wc->VPf.stropened = 1;
		  wc->ss.minX = MIN(wc->ss.minX, ith_wc->ss.minX);
		  wc->ss.maxX = MAX(wc->ss.maxX, ith_wc->ss.maxX);

		  wc->ss.minY = MIN(wc->ss.minY, ith_wc->ss.minY);
		  wc->ss.maxY = MAX(wc->ss.maxY, ith_wc->ss.maxY);

		  wc->ss.minZ = MIN(wc->ss.minZ, ith_wc->ss.minZ);
		  wc->ss.maxZ = MAX(wc->ss.maxZ, ith_wc->ss.maxZ);
		}		
	    }
	  wc->fermiContext.ntogl      = argcList;
	  wc->fermiContext.toglVector = toglVector;
	}

	/*
	  -ANTIALIAS
	*/
	else if ( strncmp(argv[i], "-antialias", 6) == 0 ) {
	  if ( Tcl_GetInt(interp, argv[i+1], &anti_alias) == TCL_ERROR ) {
	    char rss[1024];
	    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\"", argv[i+1]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	  wc->VPf.antialias = anti_alias;
	} else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), 
		   "unknown option %s, must be -antialias or -togllist", argv[i]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
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
cryDisplayMultiFS(struct Togl *togl)
{
  NEW_WIN_CONTEXT *wc;  
  wc = FindWinContextByTogl( togl );
  if (!wc->fermiContext.toglVector) return; /* no togls associated with this wc */

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
      _cryDisplayMultiFS (togl, wc);
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
	    _cryDisplayMultiFS (togl, wc);
	    glAccum(GL_ACCUM, 1.0 / (GLfloat)count);
	  }
      }
      glAccum(GL_RETURN, 1.0);
    }

  Togl_SwapBuffers(togl);
}

static void _cryDisplayMultiFS(struct Togl *togl, NEW_WIN_CONTEXT *wc) 
{
  /*GLdouble size = (GLdouble) _getSize (wc);*/
  int it;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  

  /* INSERT HERE ROTATION MATRIX */
  glLoadIdentity();
  glTranslated(0.0, 0.0, -0.5*wc->MVf.structsize); /* see crySetProj.c */
  glMultMatrixd( wc->vec.crdvec );

  /* fog */
  xcFog (togl, wc->VPf.fog, GL_FALSE);

  glDisable(GL_BLEND); 
  glDepthMask(GL_TRUE);  

  for (it=0; it<wc->fermiContext.ntogl; it++)
    {
      NEW_WIN_CONTEXT *ith_wc = FindWinContextByTogl( wc->fermiContext.toglVector[it] );
      if ( ith_wc->VPf.nsurface )    cryRenderSurface( ith_wc , FS_MULTI );
      if ( ith_wc->VPf.dispLattice ) cryRenderRecCell( ith_wc , FS_MULTI );
    }
  /*--------------------*/
  glFlush();      
}
