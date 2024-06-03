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
* Source: $XCRYSDEN_TOPDIR/C/cryTogl.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/


#include <togl.h>
#include "struct.h"
#include "xcfunc.h"

extern struct Togl *mesa_togl;

extern void crySetProjection( NEW_WIN_CONTEXT *wc, struct Togl *togl );
extern NEW_WIN_CONTEXT *FindWinContextByTogl(struct Togl *togl);

/* ------------------------------------------------------------------------
 *
 * cry_toglzoom $togl zoom-factor
 *
 * ------------------------------------------------------------------------ */
int
CRY_ToglZoomCb(ClientData clientData, Tcl_Interp *interp,
	       int argc, char *argv[])
{
  Togl *togl;
  double zoom;

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( argc != 3 ) {
    Tcl_SetResult(interp, "Usage: cry_toglzoom toglName zoom-factor", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( Tcl_GetDouble(interp, argv[2], &zoom) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\"", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if (zoom < -(1-1e-8) ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"zoom factor %f out of range", zoom);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( togl != mesa_togl ) {
    NEW_WIN_CONTEXT *wc = FindWinContextByTogl (togl);
    wc->tr.zoom *= (1.0 + zoom);
    crySetProjection (wc, togl);      
  } else {
    tr.zoom *= (1.0 + zoom);
    xcViewPort();
    LoadLights();
  }

  /* now update the display */
  Togl_PostRedisplay(togl);

  return TCL_OK; 
}
