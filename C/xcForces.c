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
* Source: $XCRYSDEN_TOPDIR/C/xcForces.c
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
#include "vector.h"
#include "xcfunc.h"

ForceVector FV = { NULL, XC_FORCE_TRESHHOLD, XC_FORCE_LENGTHFACTOR, 0.0};

extern RenderVectors *forceVectors;

typedef struct {
  char  *name;
  float (*ScaleFunc)(float value);
} parseSFunc;

static parseSFunc scaleFuncTable[] = {
  {"linear", xcLinf},
  {"log",    xcLogf},
  {"log10",  xcLog10f},
  {"sqrt",   xcSqrtf},
  {"root3",  xcRoot3f},
  {"exp",    xcExpf},
  {"exp2",   xcExp2f},
  {NULL, NULL}
};

/* --- function prototypes --- */
static void Usage(Tcl_Interp *interp);
int XC_ForcesCmd(ClientData clientData,Tcl_Interp *interp,
		 int argc, const char *argv[]);
void BuildForceVectors(ForceVector *fvPtr);

static void Usage(Tcl_Interp *interp) {
  Tcl_SetResult(interp, "Usage: xc_forces <toglName> OPTIONS\nwhere options are:\n    on|off\nor\n    scalefunction <scalefunction>\n    threshold      <threshold>\n    lengthfactor  <lengthfactor>\n", TCL_STATIC);
}


/*****************************************************************************
 * xc_forces <toglName> on|off
 *                      scalefunction <scalefunction> (see scaleFuncTable
 *                      threshold     <threshold>
 *                      lengthfactor  <lengthfactor>
 *****************************************************************************/
int 
XC_ForcesCmd(ClientData clientData,Tcl_Interp *interp,
	     int argc, const char *argv[])
{
  struct Togl *togl;
  char         c;
  int          i, len;
  double       num;

  if ( argc < 3 ) {
    Usage(interp);
    return TCL_ERROR;
  }

  /* find togl associated with toglName */
  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  c=argv[2][0];
  len=strlen(argv[2]);

  if ((c == 'o') && (strncmp(argv[2], "on", len) == 0) && (len >= 2)) {
    VPf.force = GL_TRUE;
    BuildForceVectors(&FV);
    UpdateProjection();
    return TCL_OK;
  } else if ((c == 'o') && (strncmp(argv[2], "off", len) == 0) 
	     && (len >= 2)) {
    VPf.force = GL_FALSE;
    return TCL_OK;
  } else if ((c == 's') && (strncmp(argv[2], "scalefunction", len) == 0)) {
    i=-1;
    FV.ScaleFunc = NULL;
    while ( scaleFuncTable[++i].name != NULL ) {
      if ( strcmp( scaleFuncTable[i].name, argv[3] ) == 0 ) {
	fprintf(stderr,"Force scale func: %s\n", scaleFuncTable[i].name);
	FV.ScaleFunc = scaleFuncTable[i].ScaleFunc;
      }
    }
    if (FV.ScaleFunc == NULL) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"unknown scale function %s, must be one of ...\n",
	       argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  } else if ((c == 't') && (strncmp(argv[2], "threshold", len) == 0)) {
    if ( Tcl_GetDouble(interp, argv[3], &num) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted double, but got %s\n", argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    FV.threshold = (float) num;
  } else if ((c == 'l') && (strncmp(argv[2], "lengthfactor", len) == 0)) {
    if ( Tcl_GetDouble(interp, argv[3], &num) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted double, but got %s\n", argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    FV.lengthfactor = (float) num;
  } else {
    Usage(interp);
    return TCL_ERROR;
  }

  BuildForceVectors(&FV);
  UpdateProjection();    
  return TCL_OK;
}


static double (*force)[3] = NULL;
void
BuildForceVectors(ForceVector *fvPtr)
{
  float size, sizeOrig, sizeMin = 99999.9;
  int i, j;
  float xyz[3], max[3] ={0.0, 0.0, 0.0};

  if (natoms>0) {
    force = (double (*)[3]) realloc(force,  sizeof(double [3]) * (natoms+1) );
  } else {
    return;
  }

  /*
    for (i=1; i<=natoms; i++) {
    size=(float)distdv(fv[i]);
    if( size > threshold && size > 1.0e-10 ) {
    size = ScaleFunc(size);
    if (size < sizeMin) sizeMin=size;
    }
    }
  */

  /* if force is lower then the given threshold set it to ZERO */
  for (i=1; i<=natoms; i++) 
    {
      /* first clean the force */
      for (j=0; j<3; j++) force[i][j] = 0.0;

      sizeOrig = size = (float)distdv(fv[i]);

      /*
	fprintf(stderr,"Original force: %7.4f %7.4f %7.4f; size = %7.4f\n",
	fv[i][0],fv[i][1],fv[i][2], size);
      */

      if( size > fvPtr->threshold ) 
	{

	  size = fvPtr->ScaleFunc(size) * fvPtr->lengthfactor;
	  /* if the ScaleFunc is logaritmic the size will be negative,
	     then put the size equal to size-sizeMin and augment that
	     somehow */
	  if (size < 0.0) {
	    sizeMin=fvPtr->ScaleFunc(fvPtr->threshold);
	    size=(size-sizeMin);
	  } 
	  for (j=0; j<3; j++)
	    force[i][j] = (fv[i][j]/sizeOrig)*(double)size;

	  /* 
	     fprintf(stderr," scaled force: %7.4f %7.4f %7.4f; size = %7.4f\n",
	     force[i][0],force[i][1],force[i][2], distdv(force[i]) );
	  */

	  xyz[0] = ABS(xat[i] + force[i][0]);
	  xyz[1] = ABS(yat[i] + force[i][1]);
	  xyz[2] = ABS(zat[i] + force[i][2]);

	  max[0] = xyz[0] > max[0] ? xyz[0] : max[0];
	  max[1] = xyz[1] > max[1] ? xyz[1] : max[1];
	  max[2] = xyz[2] > max[2] ? xyz[2] : max[2];
	}
    }

  fvPtr->max_size = distfv(max);

  SetForceVectorsCoor(fvPtr, force, forceVectors);
}
