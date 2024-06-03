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
* Source: $XCRYSDEN_TOPDIR/C/xcF3toI4.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <tk.h>
#include <stdio.h>
#include <math.h>
#include "struct.h"
#include "xcfunc.h"

int XC_F3toI4Cmd(ClientData clientData, Tcl_Interp *interp,
		 int argc, const char *argv[]);

/* ------------------------------------------------------------ *
 *                                                              *
 * xc_f3toi4 f1 f2 f3                                           *
 *                                                              *
 * ------------------------------------------------------------ */
int 
XC_F3toI4Cmd(ClientData clientData, Tcl_Interp *interp,
	     int argc, const char *argv[])
{
  double a1, a[3];
  double ddv, tol = 1.0e-5;
  int i, j, idv;
  int ifc[3] = {1, 1, 1};
  char *result = Tcl_Alloc( sizeof(char) * 256);

  if ( argc != 4 ) {
    Tcl_SetResult(interp, "Usage: xc_f3toi4 f1 f2 f3", TCL_STATIC);
    return TCL_ERROR;
  }

  for (i=0; i<3; i++)
    if ( Tcl_GetDouble( interp, argv[i+1], a + i) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "wanted double, but got %s", argv[i+1]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;		  
    }


  /**************************************/
  /* find an integer multiplier for a[] */  
  for (j=0; j<3; j++) {
    i = 0;
    while( 1 ) {      
      a1 = ABS( (double) ++i * a[j] );
      if ( ABS(a1 - iround(a1)) < tol ) {
	ifc[j] = i;
	break;
      }
      else if (i > iround(1.0/tol)) {
	/* some combination of triplet-numbers can produce realy huge 
	   numbers, so this loop can never end, lets say that the 
	   greater allowed I is iround(1.0/tol) */	
	tol *= 10.0;
	i = 0;
      }
    }
  }

  /* now numbers should be multiplied by ifc[0]*ifc[1]*ifc[2] */
  idv = ifc[0]*ifc[1]*ifc[2];
  ddv = (double) idv;
  /* breakpoint("breakpoint"); */
  /* we must return back the integer coordinates */
  sprintf(result, "%d %d %d %d", 
	  iround(a[0] * ddv), iround(a[1] * ddv), iround(a[2] * ddv), idv);
  Tcl_SetResult(interp, result, TCL_DYNAMIC);
  return TCL_OK;
}



