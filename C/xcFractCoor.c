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
* Source: $XCRYSDEN_TOPDIR/C/xcFractCoor.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <tk.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "xcGLparam.h"
#include "struct.h"
#include "xcfunc.h"

int XC_FractCoorCmd(ClientData clientData, Tcl_Interp *interp,
		    int argc, const char *argv[]);

/* ------------------------------------------------------------ *
 *                                                              *
 * xc_fractcoor -ctype prim|conv -coor {x y z}               *
 *                                                              *
 * ------------------------------------------------------------ */
int 
XC_FractCoorCmd(ClientData clientData, Tcl_Interp *interp,
		int argc, const char *argv[])
{
  register int i;
  double (*rvec)[4];
  float frcoor[4];
  char *result = Tcl_Alloc( sizeof(char) * 128);
  GetComOption coor  = { 0, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0} };

  if ( argc != 5 ) {
    Tcl_SetResult(interp, "Usage: xc_fractcoor -ctype prim|coor -coor {x y z}", TCL_STATIC);
    return TCL_ERROR;
  }

  for (i=1; i<argc; i+=2) {
    if ( strcmp(argv[i], "-ctype") == 0 ) {
      if ( strcmp(argv[i+1], "prim") == 0 )
	rvec = vec.recprim;
      else if ( strcmp(argv[i+1], "conv") == 0 )
	rvec = vec.recconv;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss), "unknown cell type %s, must be prim or conv", 
		 argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else if ( strcmp(argv[i], "-coor") == 0 ) {
      if ( !xcSplitList(XC_GET_XYZ, interp, &argv[i+1], &coor) ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "parse error; when parsing -coor {%s}; should be -coor {x y z}", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss),"unknown option %s, must be -ctype or -coor", 
	       argv[i]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  GetFractionalCoor(rvec, coor.vec, frcoor);

  sprintf(result, "%f %f %f", frcoor[0], frcoor[1], frcoor[2]); 
  Tcl_SetResult(interp, result, TCL_DYNAMIC);
  return TCL_OK;
}


void
GetFractionalCoor(double ivec[][4], float coor[4], float frcoor[])
{
  frcoor[0] = ivec[0][0]*coor[0] + ivec[0][1]*coor[1] + ivec[0][2]*coor[2];  
  frcoor[1] = ivec[1][0]*coor[0] + ivec[1][1]*coor[1] + ivec[1][2]*coor[2];
  frcoor[2] = ivec[2][0]*coor[0] + ivec[2][1]*coor[1] + ivec[2][2]*coor[2];
}

