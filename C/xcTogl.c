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
* Source: $XCRYSDEN_TOPDIR/C/xcTogl.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <togl.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "struct.h"
#include "xcfunc.h"


struct Togl *mesa_togl = NULL;
extern void (*xcDisplay)(struct Togl *togl);
extern NEW_WIN_CONTEXT *FindWinContextByTogl(struct Togl *togl);
extern OrthoProj ort; 
extern realTimeMove makeMovie;

extern int togl_exists;

/*
 * Togl widget create callback.  This is called by Tcl/Tk when the widget has
 * been realized.  Here's where one may do some one-time context setup or
 * initializations.
 */
int
xcToglCreateFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
{
  Togl *togl;
  static int first_time = 1;

  if (objc != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, "pathName");
    return TCL_ERROR;
  }
  if ( Togl_GetToglFromObj(interp, objv[1], &togl) != TCL_OK ) {
    return TCL_ERROR;
  }
  //Togl_GetToglFromName(interp, ".mesa", &mesa);

  /***************************************************************/
  /*   made atomic labels (only when first structure is opened   */
  /***************************************************************/
  /* first make font */
  if (first_time) {
    makeRasterFont();
    makeAtomLabels();
    makeXYZLabels();
    makeTemp3D2DList();
    /* make coor-sist */
    makeCrdList();  
    xcGenDispList();    
    first_time = 0;
  } 

  if ( strcmp(Togl_Ident(togl), ".mesa") == 0 ) {
    xcDisplayFunc(xcDummyDisplay);
    mesa_togl = togl;
  }
  else {
    NEW_WIN_CONTEXT *wc;    
    cryNewToglInit( togl );
    wc = FindWinContextByTogl( togl );
    wc->xcDisplay = xcDummyDisplay;
  }

  togl_exists = XC_TRUE;
  LoadLights();
  LoadStructMaterial();

  return TCL_OK;
}


/*
 * Togl widget reshape callback.  This is called by Tcl/Tk when the widget
 * has been resized.  Typically, we call glViewport and perhaps setup the
 * projection matrix.
 */
int
xcToglReshapeFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
{
  Togl *togl;
  int w, h;

  if (objc != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, "pathName");
    return TCL_ERROR;
  }
  if (Togl_GetToglFromObj(interp, objv[1], &togl) != TCL_OK) {
    return TCL_ERROR;
  }

  w = Togl_Width (togl);
  h = Togl_Height(togl);

  if ( togl == mesa_togl ) {
    /* old fashion style */
    VPf.width  = w;
    VPf.height = h;
    VPf.canvassize = VPf.height;
    if ( VPf.width < VPf.height) VPf.canvassize = VPf.width;

    if (VPf.stropened) {
      xcViewPort();
    } else {
      float aspect = (float) VPf.width / (float) VPf.height;
      glViewport( 0, 0, VPf.width, VPf.height );
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glFrustum(-aspect, aspect, -1.0, 1.0, 1.0, 10.0);
      glMatrixMode(GL_MODELVIEW);
    }
  } else {
    /* new nicer fashion */
    NEW_WIN_CONTEXT *wc;

    /*glViewport( 0, 0, VPf.width, VPf.height );*/
    wc = FindWinContextByTogl( togl );
    if (!wc) return TCL_ERROR;
    wc->VPf.height = w;
    wc->VPf.width  = h;

    crySetProjection( wc, togl );
  }
  return TCL_OK;
}

/*
 * this function is called when togle widget is destroyed
 */
int
xcToglDestroyFunc(ClientData clientData, Tcl_Interp *interp,
		  int objc, Tcl_Obj *const *objv)
{
  Togl *togl;

  if (objc != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, "pathName");
    return TCL_ERROR;
  }
  if (Togl_GetToglFromObj(interp, objv[1], &togl) != TCL_OK) {
    return TCL_ERROR;
  }

  if ( togl == mesa_togl ) {
    /* if this happends, then XCRYSDEN is about to exit. Do nothing. */
    /*
      xcMaybeDestroyLists(); 
      FreeAllVariables();
      mesa_togl = NULL;
    */
    xcTkFontFreeAll();
  } else {
    xcDisplay = NULL;
    DestroyWinContext( togl );
  }

  return TCL_OK;
}

/*
 * after user has stooped rotating/translating redisplay in nicer fashion 
 */
int
xcToglTimerFunc(ClientData clientData, Tcl_Interp *interp,
		int objc, Tcl_Obj *const *objv)
{
  Togl *togl;
  static int toglIsInteractive = 0;

  if (objc != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, "pathName");
    return TCL_ERROR;
  }
  if (Togl_GetToglFromObj(interp, objv[1], &togl) != TCL_OK) {
    return TCL_ERROR;
  }

  if ( togl == mesa_togl ) {

    /* interative rotation - zooming */

    if ( !toglIsInteractive && (tr.b1motion || tr.b2motion || tr.shiftB1motion) ) {
      toglIsInteractive = 1;
    }
    if ( toglIsInteractive && !(tr.b1motion || tr.b2motion || tr.shiftB1motion ) ) {
      toglIsInteractive = 0;
      Togl_PostRedisplay(togl);
    }

    /* real-time movie making */
    if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_REALTIME_INTERVAL ) {
      createMoviePPMFrame(togl);
    }
  }
  return TCL_OK;
}

/*
 * this function takes care <B1-Motion> binding -> that's ROTATION
 */
int
XC_B1MotionCb(ClientData clientData, Tcl_Interp *interp,
	      int argc, const char *argv[])
{
  int xrotnew, yrotnew; /* current <x> <y> mouse pointer position */
  int dx, dy;
  double fiX = 0, fiY = 0, ssize;
  Togl *togl;

  if (argc != 4) {
    Tcl_SetResult(interp,
		  "Usage: xc_B1motion <toglname> <x> <y>", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( Tcl_GetInt(interp, argv[2], &xrotnew) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer for <x>, but got \"%s\" in \"<toglname> xc_B1motion <x> <y>\" command", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( Tcl_GetInt(interp, argv[3], &yrotnew) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer for <y>, but got \"%s\" in \"<toglname> xc_B1motion <x> <y>\" command", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /********************************************/
  /* .MESA ---------------------------------- */
  /********************************************/
  if ( togl == mesa_togl ) {
    /* if structure is not opened, do noting, just return silently */
    if (!VPf.stropened) return TCL_OK;

    /* if B1-Motion just started, then new:=old and return silently */
    if (tr.b1motion == 0) {
      tr.b1motion = 1;
      tr.xrotold = xrotnew;
      tr.yrotold = yrotnew;
      return TCL_OK;
    } 

    /* how much did Xrot&Yrot differ from last time this routine was called */
    dx = xrotnew - tr.xrotold; 
    dy = yrotnew - tr.yrotold;
    tr.xrotold = xrotnew;
    tr.yrotold = yrotnew;

    if ( xcr.lforce ) {
      if ( MVf.structsize > 1e-5 ) ssize = MVf.structsize * VPf.VPfactor;
      else ssize = ort.size * VPf.VPfactor;
    } else {
      ssize = ort.size * VPf.VPfactor;
    }

    /* moving Mouse in X direction means rotation in Y dir !!!!! */
    /*
     * 4.0 is just some emperical factor to maximize performance of
     * mouse rotation so, that it will follow mouse pointer 
     */
    if ( dx != 0 ) fiY = atan( 4.0 * (double) dx / ssize);    
    if ( dy != 0 ) fiX = atan( 4.0 * (double) dy / ssize);
    /* now rotate */
    xcRotateXY( fiX, fiY );
    /* if dimType == XC_2D, also orientate */
    if ( dimType == XC_2D )
      hpsort_index1(tmp_nobjects, zorient, iwksp);
  } else {
    /*********************/
    /* NEW NICER FASHION */
    /*********************/
    NEW_WIN_CONTEXT *wc;
    wc = FindWinContextByTogl( togl );

    /* if structure is not opened, do noting, just return silently */
    if (!wc->VPf.stropened) return TCL_OK;

    /* if B1-Motion just started, then new:=old and return silently */
    if (wc->tr.b1motion == 0) {
      wc->tr.b1motion = 1;
      wc->tr.xrotold = xrotnew;
      wc->tr.yrotold = yrotnew;
      return TCL_OK;
    } 

    /* how much did Xrot&Yrot differ from last time this routine was called */
    dx = xrotnew - wc->tr.xrotold; 
    dy = yrotnew - wc->tr.yrotold;
    wc->tr.xrotold = xrotnew;
    wc->tr.yrotold = yrotnew;
    ssize = wc->MVf.structsize * wc->VPf.VPfactor;

    /* moving Mouse in X direction means rotation in Y dir !!!!! */
    /*
     * 4.0 is just some emperical factor to maximize performance of
     * mouse rotation so, that it will follow mouse pointer 
     */
    if ( dx != 0 ) fiY = atan( 4.0 * (double) dx / ssize);        
    if ( dy != 0 ) fiX = atan( 4.0 * (double) dy / ssize);
    cryRotateXY( wc, fiX, fiY );
  }
  Togl_PostRedisplay(togl);

  return TCL_OK;
}


/*
 * this function takes care <Shift-B1-Motion> binding -> that's ZOOMING
 *
 * It only register Shift-B1-Motion, so that the display-mode goes to
 * crude mode.
 *
 */
int
XC_ShiftB1MotionCb(ClientData clientData, Tcl_Interp *interp,
		   int argc, const char *argv[])
{
  Togl *togl;

  if (argc != 2) {
    Tcl_SetResult(interp,
		  "Usage: xc_ShiftB1motion <toglname> ", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /********************************************/
  /* .MESA ---------------------------------- */
  /********************************************/
  if ( togl == mesa_togl ) {
    /* if structure is not opened, do noting, just return silently */
    if (!VPf.stropened) return TCL_OK;

    /* if Shift-B1-Motion just started register that */
    if (tr.shiftB1motion == 0) {
      tr.shiftB1motion = 1;
    } 
  } else {
    /*********************/
    /* NEW NICER FASHION */
    /*********************/
    NEW_WIN_CONTEXT *wc;
    wc = FindWinContextByTogl( togl );

    /* if structure is not opened, do noting, just return silently */
    if (!wc->VPf.stropened) return TCL_OK;

    /* INSERT CODE FOR REGISTERING ShiftB1Motion HERE */
  }

  Togl_PostRedisplay(togl);

  return TCL_OK;
}


/*
 * this function takes care <B2-Motion> binding -> that's TRANSLATION
 */
int
XC_B2MotionCb(ClientData clientData, Tcl_Interp *interp,
	      int argc, const char *argv[])
{
  int trX, trY; /* cuurent <x> <y> mouse pointer position */
  int dx, dy;
  Togl *togl;

  if (argc != 4) {
    Tcl_SetResult(interp,
		  "Usage: xc_B2motion <toglname> <x> <y>", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( Tcl_GetInt(interp, argv[2], &trX) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer for <x>, but got \"%s\" in \"<toglname> xc_B2motion <x> <y>\" command", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( Tcl_GetInt(interp, argv[3], &trY) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer for <y>, but got \"%s\" in \"<toglname> xc_B2motion <x> <y>\" command", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }


  /********************************************/
  /* .MESA ---------------------------------- */
  /********************************************/
  if ( togl == mesa_togl ) {
    /* if structure is not opened, do noting, just return silently */
    if (!VPf.stropened) return TCL_OK;

    /* if B2-Motion just started, then new -> old and return silently */
    if (tr.b2motion == 0) {
      tr.b2motion = 1;
      tr.trXold = trX;
      tr.trYold = trY;
      return TCL_OK;
    } 
    /* how much did trX & trY differ from last time this routine was called */
    dx = trX - tr.trXold; 
    dy = trY - tr.trYold;
    tr.trXold = trX;
    tr.trYold = trY;
    /*printf("dx=%d, dy=%d", dx, dy);
      fflush(stdout);*/

    tr.xtransl += dx; /* this is used in ViewPort to determine where the
		       * structure has been translated */
    tr.ytransl -= dy; /* in X11 Y-coor is upside down */

    /* ViewPort */
    xcViewPort();
    /* Update a Display */
  } else {
    /*********************/
    /* NEW NICER FASHION */
    /*********************/
    NEW_WIN_CONTEXT *wc;
    wc = FindWinContextByTogl( togl );

    /* if structure is not opened, do noting, just return silently */
    if (!wc->VPf.stropened) return TCL_OK;

    /* INSERT CODE FOR TRANSLATION HERE */
  }

  Togl_PostRedisplay(togl);

  return TCL_OK;
}


/*
 * this function takes care of BUTTON_RELEASE events
 */
int
XC_ButtonReleaseCb(ClientData clientData, Tcl_Interp *interp,
		   int argc, const char *argv[])
{
  Togl *togl;

  if (argc != 3) {
    Tcl_SetResult(interp,
		  "Usage: xc_Brelease <toglname> <B1|B2|Shift-B1>", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /********************************************/
  /* .MESA ---------------------------------- */
  /********************************************/
  if ( togl == mesa_togl ) {  
    if ( strcmp(argv[2],"B1") == 0 || strcmp(argv[2],"Shift-B1") == 0 ) {
      tr.b1motion = 0;
      tr.shiftB1motion = 0;
    }
    else if ( strcmp(argv[2],"B2") == 0 ) tr.b2motion = 0;
    /*else if ( strcmp(argv[2],"Shift-B1") == 0 ) tr.shiftB1motion = 0;*/
    else {
      Tcl_SetResult(interp,
		    "Usage: <toglname> xc_Brelease <B1|B2|Shift-B1>", TCL_STATIC);
      return TCL_ERROR;
    }
  } else {
    /*********************/
    /* NEW NICER FASHION */
    /*********************/
    NEW_WIN_CONTEXT *wc;
    wc = FindWinContextByTogl( togl );

    if ( strcmp(argv[2],"B1") == 0 || strcmp(argv[2],"Shift-B1") == 0 )  {
      wc->tr.b1motion = 0;
      wc->tr.shiftB1motion = 0;
    }
    else if ( strcmp(argv[2],"B2") == 0 ) 
      wc->tr.b2motion = 0;
    /* else if ( strcmp(argv[2],"Shift-B1") == 0 ) wc->tr.shiftB1motion = 0; */
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss), "wrong mode %s, must be B1, B2, or Shift-B1", argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      /* Tcl_SetResult(interp,
	 "Usage: <toglname> xc_Brelease <B1|B2|Shift-B1>", TCL_STATIC);*/
      return TCL_ERROR;
    }
  }

  return TCL_OK;
}
