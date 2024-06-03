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
* Source: $XCRYSDEN_TOPDIR/C/xcBz.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <tk.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "xcfunc.h"
#include "bz.h"

typedef struct {
  int points;   /* 0/1 whether to render special points */
  int vectors;  /* 0/1 whether to render reciprocal vectors */
  int last_points; /* what was last state for points (used by BzRenderBZ - 
		      to recover when points should be rendered, if the 
		      points have been previously deleted); 
		      just for special-points */
  int last_vectors; /* the same for vectors */
} BZ_STATE;
BZ_STATE bz_state[2];

static char canvasname[2][256];
static Tk_Window canvaswin[2];
static float canvassize[2];
static float canvas_offset[2] = { BZ_OFFSET, BZ_OFFSET };

int XC_BzCmd(ClientData clientData,Tcl_Interp *interp, int argc, const char *argv[]);
static int  BzInitBZ( int typ );
static void BzViewPortBZ(Tcl_Interp *interp, int typ);
static int  BzRenderBZ(Tcl_Interp *interp, int typ);
static void BzCreateSelLine( Tcl_Interp *interp, Tcl_CmdInfo info, char *argv[], int typ, int i );
static void BzUpdateSelLine( Tcl_Interp *interp, Tcl_CmdInfo info, char *argv[], int typ, int i, char *tag );
static void BzRotateBZ( int typ, float x, float y, float z );
static void BzDegRotateBZ( int typ, float fiX, float fiY, float fiZ );
static void BzGetCoor( int typ, int pID, float coor[] );
static int  BzGetISS( int typ );
static int  BzGetMinMult( int a, int b );
/*****************************************************************************
 * xc_bz exists (prim|conv)
 *       init (prim|conv) <canvasName> ?options_in_one_piece?
 *       render   <canvasName>
 *       viewport <canvasName>
 *       rotate   <canvasName> dx dy dz
 *       zoom     <canvasName> factor
 *       degrotate <canvasName fiX fiY fiZ
 *       get      <canvasName>  (n_edgepoint | n_linepoint | n_polypoint | 
 *                               npoly | npoint | last_selected | nselected)
 *       select   <canvasName> pointnumber
 *       deselect <canvasName>
 *       deselectall <canvasName>
 *       state    <canvasName> points  (0|1)
 *                             vectors (0|1)
 *       iss      <canvasName>
 *****************************************************************************/
int
XC_BzCmd(ClientData clientData, Tcl_Interp *interp,
	 int argc, const char *argv[])
{ 
  Tk_Window tkwin = (Tk_Window) clientData;
  int i, typ;
  char *result = Tcl_Alloc( sizeof(char) * 64);
  /* at least three arguments must be present */

  /* debuging */
  /*    printf("COMMAND:: ",NULL); */
  /*    for(i=0; i<argc; i++) */
  /*      printf("%s ",argv[i]); */
  /*    printf("\n",NULL); */
  /*    fflush(stdout); */

  if ( argc < 3 ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"not enough arguments, must be one of xc_bz (exists|init|viewport|render|rotate|degrotate|zoom|get|select|deselect|state|iss); while executing %s %s %s ...", argv[0], argv[1], argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  /**********
   * EXISTS *
   **********/
  if ( strncmp(argv[1], "exis", 4) == 0 ) {    
    if ( strcmp(argv[2], "prim") == 0 )
      sprintf( result, "%d", xcr.lprimbz );      
    else
      sprintf( result, "%d", xcr.lconvbz );      
    Tcl_SetResult(interp, result, TCL_DYNAMIC);
  }
  /******** 
   * INIT *
   ********/
  else if ( strcmp(argv[1], "init") == 0 ) {
    int i;
    if ( strncmp(argv[2], "prim", 4) == 0 || 
	 strncmp(argv[2], "conv", 4) == 0 ) {
      typ = BZ_PRIMCELL;
      if ( strncmp(argv[2], "conv", 4) == 0 )
	typ = BZ_CONVCELL;
      bz[typ].isrendered = 0;
      bz[typ].recvec_is_rendered = 0;
      for (i=0; i<BZ_MAXALLPOINTS; i++) {
	bz[typ].linepoint_index[i][0]   = -1; /* because 0 means 1st polynom */
	bz[typ].linepoint_index[i][1]   = -1;
      }
      if ( !BzInitBZ(typ) ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"error when initializing rendering of Brillouin Zone; while executing %s %s %s ...", 
		 argv[0], argv[1], argv[2]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }

      bz_state[typ].points  = 1;
      bz_state[typ].vectors = 0;
      bz_state[typ].last_vectors = 0;
    }
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss),"unknown Brillouin Zone type \"%s\", must be \"prim\" or \"conv\"; while executing %s %s %s ...", 
	       argv[2], argv[0], argv[1], argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    /* now create canvas */
    if ( argc == 4 ) {
      if ( Tcl_VarEval(interp, "canvas ", argv[3], (char *)NULL) 
	   == TCL_ERROR )
	return TCL_ERROR;
    } else {
      if ( Tcl_VarEval(interp, "canvas ", argv[3], " ",
		       argv[4], (char *)NULL) == TCL_ERROR ) 
	return TCL_ERROR;
    }    
    /* if ( Tcl_VarEval(interp, "pack ", argv[3]," -fill both -expand 1", 
       (char *)NULL) == TCL_ERROR )
       return TCL_ERROR;
       if ( Tcl_VarEval(interp, "tkwait visibility ", argv[3], 
       (char *)NULL) == TCL_ERROR )
       return TCL_ERROR;	 */
    canvaswin[typ]  = Tk_NameToWindow( interp, argv[3], tkwin );
    sprintf(canvasname[typ], "%s", argv[3]);
    /* this is for point-selection */
    bz[typ].nselected = 0;
    bz[typ].n_sline   = 0;
    for (i=0; i<BZ_MAXALLPOINTS; i++) {
      bz[typ].selected_is_rendered[i] = 0;
      bz[typ].sline_is_rendered[i]    = 0;
      bz[typ].slinetype[i]            = 0;
    }
  } else {
    /* does command refer to prim or conv BZ ??? */
    if ( strcmp( argv[2], canvasname[BZ_PRIMCELL] ) == 0 )
      typ = BZ_PRIMCELL;
    else if ( strcmp( argv[2], canvasname[BZ_CONVCELL] ) == 0 )
      typ = BZ_CONVCELL;
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wrong name of bz_canvas \"%s\", must be one of \"%s\", \"%s\"; while executing %s %s %s", 
	       argv[2], canvasname[BZ_PRIMCELL], canvasname[BZ_CONVCELL],
	       argv[0], argv[1], argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    /************ 
     * VIEWPORT *
     ************/
    if ( strcmp( argv[1], "viewport" ) == 0 )
      BzViewPortBZ(interp, typ);
    /********** 
     * RENDER *
     **********/
    else if ( strcmp( argv[1], "render" ) == 0 )
      BzRenderBZ(interp, typ);
    /********** 
     * ROTATE *
     **********/
    else if ( strcmp( argv[1], "rotate" ) == 0 ) {
      double dx, dy, dz;
      float x, y, z;
      if ( Tcl_GetDouble(interp, argv[3], &dx) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "expected double, but got \"%s\"; while executing %s %s %s %s %s %s", 
		 argv[3], argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      if ( Tcl_GetDouble(interp, argv[4], &dy) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "expected double, but got \"%s\"; while executing %s %s %s %s %s %s", 
		 argv[4], argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      if ( Tcl_GetDouble(interp, argv[5], &dz) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "expected double, but got \"%s\"; while executing %s %s %s %s %s %s", 
		 argv[5], argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      x = (float) dx;
      y = (float) dy;
      z = (float) dz;      
      BzRotateBZ(typ, x, y, z);
      BzViewPortBZ(interp, typ);
      BzRenderBZ(interp, typ);
    }
    /************* 
     * DEGROTATE *
     *************/
    else if ( strcmp( argv[1], "degrotate" ) == 0 ) {
      double dfiX, dfiY, dfiZ;
      float fiX, fiY, fiZ;
      if ( Tcl_GetDouble(interp, argv[3], &dfiX) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "expected double, but got \"%s\"; while executing %s %s %s %s %s %s", 
		 argv[3], argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      if ( Tcl_GetDouble(interp, argv[4], &dfiY) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "expected double, but got \"%s\"; while executing %s %s %s %s %s %s", 
		 argv[4], argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      if ( Tcl_GetDouble(interp, argv[5], &dfiZ) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "expected double, but got \"%s\"; while executing %s %s %s %s %s %s", 
		 argv[5], argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      fiX = (float) dfiX / RAD2DEG;
      fiY = (float) dfiY / RAD2DEG;
      fiZ = (float) dfiZ / RAD2DEG;      
      BzDegRotateBZ(typ, fiX, fiY, fiZ);
      BzViewPortBZ(interp, typ);
      BzRenderBZ(interp, typ);
    }
    /********
     * ZOOM *
     ********/    
    else if ( strcmp( argv[1], "zoom") == 0 ) {
      double f;
      if ( Tcl_GetDouble(interp, argv[3], &f) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "expected double, but got \"%s\"; while executing %s %s %s %s", 
		 argv[3], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      /* greater offset -> minor BZ */
      printf("zoom step == %f\n", f);
      canvas_offset[typ] *= (float) (1.0 - f);
      BzViewPortBZ(interp, typ);
      BzRenderBZ(interp, typ);
    }
    /*******
     * GET *
     *******/
    else if ( strcmp( argv[1], "get") == 0 ) {
      if ( strncmp( argv[3], "n_edg", 5) == 0 ) {
	sprintf(result, "%d", bz[typ].n_edgepoint);
      }
      else if ( strncmp( argv[3], "n_lin", 5) == 0 ) {
	sprintf(result, "%d", bz[typ].n_linepoint);
      }   
      else if ( strncmp( argv[3], "n_pol", 5) == 0 ) {
	sprintf(result, "%d", bz[typ].n_polypoint);
      }   
      else if ( strncmp( argv[3], "npoly", 5) == 0 ) {
	sprintf(result, "%d", bz[typ].npoly);
      }   
      else if ( strncmp( argv[3], "npoin", 5) == 0 ) {
	sprintf(result, "%d", bz[typ].npoint);
      }
      else if ( strncmp( argv[3], "nsele", 5) == 0 ) {
	sprintf(result, "%d", bz[typ].nselected);
      }
      else if ( strncmp( argv[3], "last_", 5) == 0 ) {
	if ( bz[typ].nselected > 0 )
	  sprintf(result, "%d", bz[typ].selectedID[ bz[typ].nselected - 1]);
	else
	  sprintf(result, "-1");
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss), "unknown option \"%s\", must be one of n_edgepoint, n_linepoint, n_polypoint, npoly, npoint, nselected, last_selected; while executing %s %s %s %s",
		 argv[3], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      Tcl_SetResult(interp, result, TCL_DYNAMIC);
    }
    /**********
     * SELECT *
     **********/
    else if ( strcmp( argv[1], "select" ) == 0 ) {
      int pID;
      float *coor;
      if ( Tcl_GetInt(interp, argv[3], &pID) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "expected int, but got \"%s\"; while executing %s %s %s %s", argv[3], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      if ( pID < 0 || pID >= bz[typ].npoint) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "selected point #\"%s\" is out of range, must be within [0,%d); while executing %s %s %s %s", 
		 argv[3], bz[typ].npoint, argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      bz[typ].selectedID[ bz[typ].nselected ] = pID;
      /* get the coordinates of selected points */
      BzGetCoor( typ, pID, bz[typ].reccoor[bz[typ].nselected] );
      bz[typ].nselected++;

      if ( bz[typ].nselected > 1 ) {
	int i, j;
	int isl = bz[typ].n_sline;
	/* get line coordinates */
	for (j=0; j<2; j++) {
	  for (i=0; i<3; i++) {
	    bz[typ].sline[isl][j*3+i] = 
	      bz[typ].point[ bz[typ].selectedID[bz[typ].nselected-2+j] ][i];
	    bz[typ].rot_sline[isl][j*3+i] =
	      bz[typ].point[ bz[typ].selectedID[bz[typ].nselected-2+j] ][i];
	  }
	  bz[typ].slinetype[isl] |= 
	    bz[typ].pointtype[ bz[typ].selectedID[bz[typ].nselected-2+j] ];
	  if ( bz[typ].pointtype[ bz[typ].selectedID[bz[typ].nselected-2+j] ] == BZ_LINEPOINT ) {
	    bz[typ].sline_index[isl][0] = bz[typ].linepoint_index[ bz[typ].selectedID[bz[typ].nselected-2+j] ][0];
	    bz[typ].sline_index[isl][1] = bz[typ].linepoint_index[ bz[typ].selectedID[bz[typ].nselected-2+j] ][1];
	  }
	}
	bz[typ].nZpoint++;
	for (i=0; i<3; i++)
	  bz[typ].slcenter[isl][i] = 
	    (bz[typ].sline[isl][i] + bz[typ].sline[isl][i+3]) / 2.0;
	bz[typ].n_sline++;
      }
      BzDegRotateBZ(typ, 0.0, 0.0, 0.0);
      BzViewPortBZ(interp, typ);
      BzRenderBZ(interp, typ);
      /* update the Bz($can,nselected) tcl variable, 
	 and the coor of selected points; 
	 thatwhy we will return their values 
      */
      coor = bz[typ].reccoor[bz[typ].nselected-1];
      sprintf(result, "%d %f %f %f", 
	      bz[typ].nselected, coor[0], coor[1], coor[2]);
      Tcl_SetResult(interp, result, TCL_DYNAMIC);
    }
    /************
     * DESELECT *
     ************/
    else if ( strcmp( argv[1], "deselect" ) == 0 ) {
      char tag[10];
      if ( bz[typ].nselected > 1 ) {
	bz[typ].n_sline--;
	bz[typ].nZpoint--;
	bz[typ].sline_is_rendered[ bz[typ].n_sline ] = 0;
	sprintf(tag, "sln%d", bz[typ].n_sline);
	/* delete selline */
	if ( Tcl_VarEval( interp, canvasname[typ], " delete ", 
			  tag, (char *)NULL ) == TCL_ERROR )
	  return TCL_ERROR;
      }
      if ( bz[typ].nselected > 0 ) {
	bz[typ].nselected--;
	bz[typ].selected_is_rendered[ bz[typ].nselected ] = 0;
	sprintf(tag, "spt%d", bz[typ].nselected);
	if ( Tcl_VarEval( interp, canvasname[typ], " delete ", 
			  tag, (char *)NULL ) == TCL_ERROR )
	  return TCL_ERROR;
	BzRenderBZ(interp, typ);
	/* update the Bz($can,nselected) tcl variable, 
	   thatwhy we will return its value 
	*/
      }
      sprintf(result, "%d", bz[typ].nselected);
      Tcl_SetResult(interp, result, TCL_DYNAMIC);
    }
    /***************
     * DESELECTALL *
     ***************/
    else if ( strcmp( argv[1], "deselectall" ) == 0 ) {
      /* first handle sellines */
      for (i=0; i<bz[typ].n_sline; i++)
	bz[typ].sline_is_rendered[i] = 0;
      bz[typ].n_sline = 0;
      bz[typ].nZpoint = bz[typ].npoint;      
      /* delete selline */
      if ( Tcl_VarEval( interp, canvasname[typ], " delete ", 
			"selline", (char *)NULL ) == TCL_ERROR )
	return TCL_ERROR;

      /* now handle selected points */
      for (i=0; i<bz[typ].nselected; i++)
	bz[typ].selected_is_rendered[i] = 0;      
      bz[typ].nselected = 0;
      if ( Tcl_VarEval( interp, canvasname[typ], " delete ", 
			"selpoint", (char *)NULL ) == TCL_ERROR )
	return TCL_ERROR;

      BzRenderBZ(interp, typ);
      /* update the Bz($can,nselected) tcl variable, 
	 thatwhy we will return its value 
      */    
      strcpy(result, "0");
      Tcl_SetResult(interp, result, TCL_DYNAMIC);
    }
    /*********
     * STATE *
     *********/
    else if ( strcmp( argv[1], "state" ) == 0 ) {
      int state;
      if ( Tcl_GetInt(interp, argv[4], &state) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "expected int, but got \"%s\"; while executing %s %s %s %s %s", argv[3], argv[0], argv[1], argv[2], argv[3], argv[4]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      if ( strncmp( argv[3], "point", 5 ) == 0 ) {
	bz_state[typ].points = state;
      }
      else if ( strncmp( argv[3], "vecto", 5 ) == 0 ) {
	bz_state[typ].vectors = state;
	/* vectors are larger than BZ */
	/* INSERT: update bz[typ].max */
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", must be one of points, vectors; while executing %s %s %s %s %s", 
		 argv[3], argv[0], argv[1], argv[2], argv[3], argv[4]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      BzDegRotateBZ( typ, 0.0, 0.0, 0.0 );
      BzViewPortBZ( interp, typ);
      BzRenderBZ( interp, typ );
    }
    /*******
     * ISS *
     *******/
    else if ( strcmp( argv[1], "iss" ) == 0 ) {
      i = BzGetISS(typ);
      sprintf(result, "%d", i);
      Tcl_SetResult(interp, result, TCL_DYNAMIC);
    }
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss),"unknown command \"xc_bz %s\", must be one of xc_bz (init|viewport|render|rotate|degrotate|select|deselect|deselectall|state|get|iss); while executing %s %s %s ...", argv[1],
	       argv[0], argv[1], argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  return TCL_OK;
}


static int 
BzInitBZ( int typ ) 
{
  int   i, j, k, l, jj, isequal, npoi = 0;
  float sumx, sumy, sumz;
  float edgepoint[BZ_MAXPOINTS][3], linepoint[BZ_MAXPOINTS][3];
  int   lindex[BZ_MAXPOINTS];

  /********************************/
  /* check if BZ[typ] does exists */
  if ( !xcr.lprimbz && !xcr.lconvbz ) {
    printf("BZ:: prim=%d, conv=%d\n", xcr.lprimbz, xcr.lconvbz); 
    return XC_ERROR;
  }

  /***************************/
  /* copy reciprocal vectors */
  for (i=0; i<xcr.dim; i++)
    for (j=0; j<3; j++) {
      if ( typ == BZ_PRIMCELL)
	bz[typ].recvec[i][j] = (float) vec.recprim[i][j];
      else
	bz[typ].recvec[i][j] = (float) vec.recconv[i][j];
      bz[typ].rot_recvec[i][j] = bz[typ].recvec[i][j];
    }

  /***************************************************************************/
  /* copy polygons to rot_poly */
  for (i=0; i<bz[typ].npoly; i++)
    for (j=0; j<bz[typ].nvert[i]; j++)
      for (k=0; k<3; k++)
	bz[typ].rot_poly[i][j][k] = bz[typ].poly[i][j][k];

  /* get special points */
  bz[typ].n_edgepoint = 0;
  bz[typ].n_linepoint = 0;
  bz[typ].n_polypoint = 0;
  /* get centerpoint */
  bz[typ].point[npoi][0] = 0.0;
  bz[typ].point[npoi][1] = 0.0;
  bz[typ].point[npoi][2] = 0.0;
  bz[typ].rot_points[npoi][0] = 0.0;
  bz[typ].rot_points[npoi][1] = 0.0;
  bz[typ].rot_points[npoi][2] = 0.0;
  bz[typ].pointtype[npoi] = BZ_CENTERPOINT;
  npoi++;

  for (i=0; i<bz[typ].npoly; i++) {
    sumx = 0.0;
    sumy = 0.0;
    sumz = 0.0;
    for (j=0; j<bz[typ].nvert[i]; j++) {
      isequal = 0;
      sumx += bz[typ].poly[i][j][0];
      sumy += bz[typ].poly[i][j][1];
      sumz += bz[typ].poly[i][j][2];

      /* check if the EDGEPOINT is already registered */
      for (k=0; k<bz[typ].n_edgepoint; k++) {
	if ( IsSamePointvf( bz[typ].poly[i][j],
			    edgepoint[bz[typ].n_edgepoint-1],
			    1.e-7) ) {
	  isequal = 1;
	}
      }
      if ( !isequal ) {
	for (l=0; l<3; l++) {
	  edgepoint[bz[typ].n_edgepoint][l] =
	    bz[typ].poly[i][j][l];
	  bz[typ].point[npoi][l] = bz[typ].poly[i][j][l];
	  bz[typ].rot_points[npoi][l] =
	    bz[typ].poly[i][j][l];
	}
	bz[typ].pointtype[npoi] = BZ_EDGEPOINT;
	npoi++;
	bz[typ].n_edgepoint++;
      }

      /* get new LINEPOINT */
      if ( j < bz[typ].nvert[i] - 1) {	 
	jj = j + 1;
      } else {
	jj = 0;
      }
      for (l=0; l<3; l++) {
	linepoint[bz[typ].n_linepoint][l] = 
	  ( bz[typ].poly[i][j][l] +  bz[typ].poly[i][jj][l] ) / 2.0;
      }
      /* check if the linepoint is already registered */
      isequal = 0;
      for (k=0; k<bz[typ].n_linepoint; k++) {
	if ( IsSamePointvf( linepoint[bz[typ].n_linepoint],
			    linepoint[k], 1.e-7) ) {
	  isequal = 1;
	  bz[typ].linepoint_index[ lindex[k] ][1] = i;
	}
      }
      if ( !isequal ) {
	for (l=0; l<3; l++) {
	  bz[typ].point[npoi][l] =
	    linepoint[bz[typ].n_linepoint][l];
	  bz[typ].rot_points[npoi][l] =
	    linepoint[bz[typ].n_linepoint][l];
	}
	bz[typ].pointtype[npoi] = BZ_LINEPOINT;
	lindex[bz[typ].n_linepoint] = npoi;
	bz[typ].linepoint_index[npoi][0] = i;	  
	npoi++;
	bz[typ].n_linepoint++;	
      }
    }

    /* get new POLYPOINT */
    bz[typ].point[npoi][0] = sumx / (float) bz[typ].nvert[i];
    bz[typ].point[npoi][1] = sumy / (float) bz[typ].nvert[i];
    bz[typ].point[npoi][2] = sumz / (float) bz[typ].nvert[i];

    bz[typ].rot_points[npoi][0] = sumx / (float) bz[typ].nvert[i];
    bz[typ].rot_points[npoi][1] = sumy / (float) bz[typ].nvert[i];
    bz[typ].rot_points[npoi][2] = sumz / (float) bz[typ].nvert[i];

    bz[typ].n_polypoint++;
    bz[typ].polyindex[npoi] = i;
    printf("npoi=%d; polygon: %d\n", npoi, i);
    bz[typ].pointtype[npoi] = BZ_POLYPOINT;
    npoi++;
  }

  /* initialize rotation matrix */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      bz[typ].rotmat[i][j] = 0.0;
  bz[typ].rotmat[0][0] = 1.0;
  bz[typ].rotmat[1][1] = 1.0;
  bz[typ].rotmat[2][2] = 1.0;

  /* get number of Z-points */
  bz[typ].npoint  = npoi;
  bz[typ].nZpoint = npoi;

  /*    for (i=0; i<bz[typ].n_linepoint; i++) */
  /*      printf("LINE:: i=%d, lindex=%d, line_index={%d, %d}\n", */
  /*  	   i, lindex[i], */
  /*  	   bz[typ].linepoint_index[lindex[i]][0], */
  /*  	   bz[typ].linepoint_index[lindex[i]][1]); */
  /*    fflush(stdout); */

  return XC_OK;
}


static void
BzViewPortBZ(Tcl_Interp *interp, int typ)
{
  float w, h, f;
  int i, j;

  w = (float) Tk_Width( canvaswin[typ] );
  h = (float) Tk_Height( canvaswin[typ] );
  canvassize[typ] = w;
  if ( h < canvassize[typ] )
    canvassize[typ] = h;

  f = 0.5 * (1.0 - canvas_offset[typ]) * canvassize[typ] / bz[typ].max;

  /* WARNING: take care of X y-coor, 
     because y-coor in X are turned upside down */
  for (i=0; i<xcr.dim; i++) {    
    bz[typ].vp_recvec[i][0] = w/2.0;
    bz[typ].vp_recvec[i][1] = h/2.0;
    bz[typ].vp_recvec[i+3][0] = w/2.0 + f * bz[typ].rot_recvec[i][0] / 2.0;
    bz[typ].vp_recvec[i+3][1] = h/2.0 - f * bz[typ].rot_recvec[i][1] / 2.0;    
    bz[typ].vp_recvec[i+6][0] = w/2.0 + f * bz[typ].rot_recvec[i][0];
    bz[typ].vp_recvec[i+6][1] = h/2.0 - f * bz[typ].rot_recvec[i][1];    
  }
  for(i=0; i<bz[typ].npoly; i++)
    for(j=0; j<bz[typ].nvert[i]; j++) {     
      bz[typ].vp_poly[i][j][0] = w/2.0 + f * bz[typ].rot_poly[i][j][0];
      bz[typ].vp_poly[i][j][1] = h/2.0 - f * bz[typ].rot_poly[i][j][1];
    }
  for(i=0; i<bz[typ].npoint; i++) {
    bz[typ].vp_points[i][0] = w/2.0 + f * bz[typ].rot_points[i][0];
    bz[typ].vp_points[i][1] = h/2.0 - f * bz[typ].rot_points[i][1];
  }
  for (i=0; i<bz[typ].n_sline; i++) {
    bz[typ].vp_sline[i][0] = w/2.0 + f * bz[typ].rot_sline[i][0];
    bz[typ].vp_sline[i][1] = h/2.0 - f * bz[typ].rot_sline[i][1];
    bz[typ].vp_sline[i][2] = w/2.0 + f * bz[typ].rot_sline[i][3];
    bz[typ].vp_sline[i][3] = h/2.0 - f * bz[typ].rot_sline[i][4];
  }
}


static void
BzRotateBZ( int typ, float x, float y, float z )
{
  float fiX, fiY, fiZ;

  /* 4.0 is just some emperical factor to maximize performance of
   * mouse rotation so, that it will follow mouse pointer 
   */
  if ( x != 0 ) fiY = atan( 4.0 * x / canvassize[typ] );    
  else fiY = 0.0;

  if ( y != 0 ) fiX = atan( 4.0 * y / canvassize[typ] );
  else fiX = 0.0;

  if ( z != 0 ) fiZ = atan( 4.0 * z / canvassize[typ] );
  else fiZ = 0.0;

  BzDegRotateBZ( typ, fiX, fiY, fiZ );
}


static void
BzDegRotateBZ( int typ, float fiX, float fiY, float fiZ )
{
  int i,j;
  float cosfiX, cosfiY, cosfiZ, sinfiX, sinfiY, sinfiZ;
  /* 2.2 is just to mark where a cos() or sin() will go */
  float rotx[3][3] = {
    { 1.0, 0.0, 0.0},
    { 0.0, 2.2, 2.2},
    { 0.0, 2.2, 2.2},
  };
  float roty[3][3] = {
    { 2.2, 0.0, 2.2},
    { 0.0, 1.0, 0.0},
    { 2.2, 0.0, 2.2},
  };
  float rotz[3][3] = {
    { 2.2, 2.2, 0.0},
    { 2.2, 2.2, 0.0},
    { 0.0, 0.0, 1.0},
  };
  float rotxy[3][3], rotxyz[3][3], oldrot[3][3];

  cosfiX = cos(fiX);
  sinfiX = sin(fiX);
  cosfiY = cos(fiY);
  sinfiY = sin(fiY);
  cosfiZ = cos(fiZ);
  sinfiZ = sin(fiZ);

  rotx[1][1] = cosfiX;
  rotx[1][2] = -sinfiX;
  rotx[2][1] = sinfiX;
  rotx[2][2] = cosfiX;

  roty[0][0] = cosfiY;
  roty[0][2] = sinfiY;
  roty[2][0] = -sinfiY;
  roty[2][2] = cosfiY;

  rotz[0][0] = cosfiZ;
  rotz[0][1] = -sinfiZ;
  rotz[1][0] = sinfiZ;
  rotz[1][1] = cosfiZ;

  /* now multiply the matrices */
  MatCopy33f(bz[typ].rotmat, oldrot);
  MatToZero33f(rotxy);
  MatMult33f(rotx, roty, rotxy);
  MatToZero33f(rotxyz);
  MatMult33f(rotxy, rotz, rotxyz);
  MatToZero33f(bz[typ].rotmat);
  MatMult33f(rotxyz, oldrot, bz[typ].rotmat);

  /* now update the reciprocal vectors */
  for (i=0; i<xcr.dim; i++)
    MatVecMult33f(bz[typ].rotmat, bz[typ].recvec[i], bz[typ].rot_recvec[i]);

  /* now update the polygons */
  for (i=0; i<bz[typ].npoly; i++)
    for (j=0; j<bz[typ].nvert[i]; j++)
      MatVecMult33f(bz[typ].rotmat, bz[typ].poly[i][j], 
		    bz[typ].rot_poly[i][j]);

  /* now update the points */
  for (i=0; i<bz[typ].npoint; i++)
    MatVecMult33f(bz[typ].rotmat, bz[typ].point[i], bz[typ].rot_points[i] );

  /* now update the selected lines */
  for (i=0; i<bz[typ].n_sline; i++) {
    MatVecMult33f(bz[typ].rotmat,
		  bz[typ].sline[i], bz[typ].rot_sline[i]);
    MatVecMult33f(bz[typ].rotmat,
		  &(bz[typ].sline[i][3]), &(bz[typ].rot_sline[i][3]));
    MatVecMult33f(bz[typ].rotmat, 
		  bz[typ].slcenter[i], bz[typ].rot_slcenter[i]);
  }
}


static void
BzGetCoor( int typ, int pID, float v[] ) 
{
  double (*vc)[4];

  if ( typ == BZ_PRIMCELL )
    vc = (double (*)[4]) vec.prim;
  else
    vc = (double (*)[4]) vec.conv;

  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;
  if ( xcr.dim == 1 ) 
    v[0] = vc[0][0] * bz[typ].point[pID][0];
  else if ( xcr.dim == 2 ) {
    v[0] = 
      vc[0][0] * bz[typ].point[pID][0] + 
      vc[0][1] * bz[typ].point[pID][1];
    v[1] = 
      vc[1][0] * bz[typ].point[pID][0] + 
      vc[1][1] * bz[typ].point[pID][1];
  }			    
  else if ( xcr.dim == 3 ) {
    v[0] = 
      vc[0][0] * bz[typ].point[pID][0] + 
      vc[0][1] * bz[typ].point[pID][1] +
      vc[0][2] * bz[typ].point[pID][2];
    v[1] = 
      vc[1][0] * bz[typ].point[pID][0] + 
      vc[1][1] * bz[typ].point[pID][1] +
      vc[1][2] * bz[typ].point[pID][2];
    v[2] = 
      vc[2][0] * bz[typ].point[pID][0] + 
      vc[2][1] * bz[typ].point[pID][1] +
      vc[2][2] * bz[typ].point[pID][2];
  }

}

#define BZRENDERBZ_ARGC 100

static int 
BzRenderBZ(Tcl_Interp *interp, int typ)
{
  /* extern void * malloc(); */
  Tcl_CmdInfo info;
  int i, j, ii, ip, argc, ipoly, npoi;
  float ps, ps1;
  float pointZ[1000];
  int   iorder[1000];
  char *argv[BZRENDERBZ_ARGC], *tag, *cmd, *color;

  for(i=0; i<BZRENDERBZ_ARGC; i++)
    argv[i] = (char *) malloc( sizeof(char) * 64 );

  cmd = canvasname[typ];  
  strcpy(argv[0], cmd);

  Tcl_ResetResult(interp);
  if ( !Tcl_GetCommandInfo(interp, cmd, &info) ) {
    Tcl_AppendResult(interp, "Unknown command \"", cmd, "\"", NULL);
    return TCL_ERROR;
  }
  /*Tcl_ResetResult( interp );*/

  tag    = (char *) malloc( sizeof(char) * 20);
  /*
    pointZ = (float *) malloc( sizeof(float) * (bz[typ].nZpoint + 1 + xcr.dim));
    iorder = (int *)   malloc( sizeof(int) * (bz[typ].nZpoint + 1 + xcr.dim));
  */

  printf("ADDRESS: pointZ == %p, iorder == %p\n", 
	 (void *) pointZ, (void *) iorder);
  /* copy points to pointZ[] */  
  for(i=1; i<=bz[typ].npoint; i++)
    pointZ[i] = bz[typ].rot_points[i-1][2];
  for(i=bz[typ].npoint+1; i<=bz[typ].nZpoint; i++)
    pointZ[i] = bz[typ].rot_slcenter[i-(1+bz[typ].npoint)][2];
  for(i=bz[typ].nZpoint+1; i<=bz[typ].nZpoint+xcr.dim; i++) {
    ii = i - bz[typ].nZpoint - 1;
    pointZ[i]         = 0.25 * bz[typ].rot_recvec[ii][2];
    pointZ[i+xcr.dim] = 1.0  * bz[typ].rot_recvec[ii][2];
  }

  /* Z-orientate */
  npoi = bz[typ].nZpoint + 2 * xcr.dim;
  hpsort_index1_f( npoi , pointZ, iorder );

  ipoly = 0;
  if ( !bz[typ].isrendered ) {
    printf("CREATE ITEMS\n");
    /****************
     * CREATE ITEMS *
     ****************/
    for(i=0; i<npoi; i++) {
      if ( iorder[i + 1] <= bz[typ].npoint ) {
	ii = iorder[i + 1] - 1;
	printf("\nRENDERING %dth object;  iorder == %d,  pointtype == %d\n", 
	       i, ii, bz[typ].pointtype[ii]);
	fflush(stdout);
	if ( bz[typ].pointtype[ii] == BZ_POLYPOINT ) {
	  /*
	   * first render POLYGON, then render points 
	   */
	  argc = 1;
	  WriteArgv( argv, &argc, "ss", "create", "line" );
	  ip = bz[typ].polyindex[ii];
	  printf("ii=%d; polygon number: %d\n",ii, ip);		
	  for(j=0; j<bz[typ].nvert[ip]; j++)
	    WriteArgv( argv, &argc, "ff", 
		       bz[typ].vp_poly[ip][j][0], bz[typ].vp_poly[ip][j][1]);
	  WriteArgv( argv, &argc, "ff", 
		     bz[typ].vp_poly[ip][0][0], bz[typ].vp_poly[ip][0][1]);

	  sprintf(tag, "polygon p%d", ip);

	  ps    = BZ_THIN_POINT_SIZE;
	  if ( (float) (ipoly + 1) / bz[typ].npoly >= 0.49 ) {
	    ps    = BZ_BOLD_POINT_SIZE;
	  }
	  ipoly++;

	  WriteArgv( argv, &argc, "sssfss", "-fill", "#000000",
		     "-width", ps, "-tags", tag);
	  Tcl_ResetResult( interp );
	  (*info.proc)(info.clientData, interp, argc, (const char**) argv);	
	  //printf("LINE:: polyID:: %s p%d\n",interp->result,ip);

	  if ( xcr.dim == 3) {
	    /* now render POLYPOINT */
	    argc = 1;
	    WriteArgv( argv, &argc, "ssffffffff", "create", "polygon",
		       bz[typ].vp_points[ii][0] - ps,
		       bz[typ].vp_points[ii][1],
		       bz[typ].vp_points[ii][0],
		       bz[typ].vp_points[ii][1] - ps,
		       bz[typ].vp_points[ii][0] + ps,
		       bz[typ].vp_points[ii][1],
		       bz[typ].vp_points[ii][0],
		       bz[typ].vp_points[ii][1] + ps);

	    printf("OBJECT DETAILS: polypoint, argc=%d\n", argc);

	    sprintf(tag, "point pt%d", ii);

	    WriteArgv( argv, &argc, "ssssss", "-fill", "#000000", "-width", 
		       "1", "-tags", tag);
	    xcdebug("before info.proc");
	    (*info.proc)(info.clientData, interp, argc, (const char **) argv);
	    xcdebug("after info.proc");
	  }
	} else {
	  /*
	   * RENDER non POLYPOINT
	   */
	  if ( !( xcr.dim == 1 &&  bz[typ].pointtype[ii] == BZ_LINEPOINT ) ) {
	    /* for xcr.dim == 1 the only BZ_LINEPOINT is BZ_CENTERPOINT */
	    argc = 1;
	    WriteArgv( argv, &argc, "ss", "create", "oval");

	    ps = BZ_THIN_POINT_SIZE;
	    color = "#000000";
	    if ( (float) i / bz[typ].npoint >= 0.49 )
	      ps = BZ_BOLD_POINT_SIZE;	    
	    if ( xcr.dim < 3 && bz[typ].rot_points[ii][2] > -1.0e-6 )
	      ps = BZ_BOLD_POINT_SIZE;
	    if ( bz[typ].pointtype[ii] == BZ_CENTERPOINT ) {
	      ps = BZ_CENTERPOINT_SIZE;
	      color = "#ff0000";
	    }

	    WriteArgv( argv, &argc, "ffff", 
		       bz[typ].vp_points[ii][0] - ps,
		       bz[typ].vp_points[ii][1] - ps,
		       bz[typ].vp_points[ii][0] + ps,
		       bz[typ].vp_points[ii][1] + ps);

	    printf("OBJECT DETAILS: point\n, argc=%d\n", argc);

	    sprintf(tag, "point pt%d", ii );
	    WriteArgv( argv, &argc, "ssssssss", "-fill", color, 
		       "-outline", color, "-width", "1", "-tags", tag);
	    xcdebug("before info.proc");
	    (*info.proc)(info.clientData, interp, argc, (const char**) argv);	
	    xcdebug("after info.proc");		    
	  }
	}
      }
    }
  } else {
    float vp_points[2];
    /**********************
     * UPDATE COORDINATES *
     **********************/
    for(i=0; i<npoi; i++) {
      if ( iorder[i+1] <= bz[typ].npoint ) {
	ii = iorder[i + 1] - 1;
	if ( bz_state[typ].points ) {
	  /* points are turned on */
	  vp_points[0] = bz[typ].vp_points[ii][0];
	  vp_points[1] = bz[typ].vp_points[ii][1];
	} else {
	  vp_points[0] = -100.0;
	  vp_points[1] = -100.0;
	}

	/*  	printf("LOOP --> i=%d, ii=%d, typ ==%d\n",  */
	/*  	       i, ii, bz[typ].pointtype[ii]); */
	/*  	fflush(stdout); */

	if ( bz[typ].pointtype[ii] == BZ_POLYPOINT ) {
	  /* first update POLYGON coordinates, then update points */
	  ip = bz[typ].polyindex[ii];
	  argc = 1;
	  sprintf(tag, "p%d", ip);
	  WriteArgv( argv, &argc, "ss", "coords", tag);

	  /*  	  printf("ii=%d; polygon number: %d\n",ii, ip); */
	  /*  	  fflush(stdout); */

	  for(j=0; j<bz[typ].nvert[ip]; j++)
	    WriteArgv( argv, &argc, "ff",      
		       bz[typ].vp_poly[ip][j][0], bz[typ].vp_poly[ip][j][1]);
	  WriteArgv( argv, &argc, "ff",      
		     bz[typ].vp_poly[ip][0][0], bz[typ].vp_poly[ip][0][1]);
	  (*info.proc)(info.clientData, interp, argc, (const char**) argv);

	  ps  = BZ_THIN_POINT_SIZE;
	  ps1 = BZ_THIN_POLYGON_WIDTH;	  
	  if ( (float) (ipoly + 1) / bz[typ].npoly > 0.5 ) {
	    ps  = BZ_BOLD_POINT_SIZE;
	    ps1 = BZ_BOLD_POLYGON_WIDTH;
	  }
	  ipoly++;

	  /* update the width & outline of polygons */
	  argc = 1;
	  WriteArgv( argv, &argc, "sssf", "itemconfigure", tag, 
		     "-width",  ps1);
	  (*info.proc)(info.clientData, interp, argc, (const char**) argv);	

	  /* change the stack order */
	  argc = 1;
	  WriteArgv( argv, &argc, "sss", "raise", tag, "all");
	  (*info.proc)(info.clientData, interp, argc, (const char**) argv);

	  if ( xcr.dim == 3 ) {
	    /* update POLYPOINT */
	    argc = 1;
	    sprintf(tag,"pt%d",ii);
	    WriteArgv( argv, &argc, "ssffffffff", "coords", tag,
		       vp_points[0] - ps, vp_points[1],
		       vp_points[0],      vp_points[1] - ps,
		       vp_points[0] + ps, vp_points[1],
		       vp_points[0],      vp_points[1] + ps);
	    (*info.proc)(info.clientData, interp, argc, (const char**) argv);

	    /* change the stack order */
	    argc = 1;
	    WriteArgv( argv, &argc, "sss", "raise", tag, "all");
	    (*info.proc)(info.clientData, interp, argc, (const char**) argv);
	  }
	} else {
	  /* 
	   * non POLYPOINT 
	   */
	  if ( !( xcr.dim == 1 &&  bz[typ].pointtype[ii] == BZ_LINEPOINT ) ) {
	    /* for xcr.dim == 1 the only BZ_LINEPOINT is BZ_CENTERPOINT */
	    ps  = BZ_THIN_POINT_SIZE;
	    ps1 = BZ_THIN_POLYGON_WIDTH;
	    color = "#000000";
	    if ( (float) i / bz[typ].npoint > 0.5 ) {
	      ps  = BZ_BOLD_POINT_SIZE;
	      ps1 = BZ_BOLD_POLYGON_WIDTH;
	    }
	    if ( xcr.dim < 3 && bz[typ].rot_points[ii][2] > -1.0e-6 )
	      ps = BZ_BOLD_POINT_SIZE;
	    if ( bz[typ].pointtype[ii] == BZ_CENTERPOINT ) {
	      color = "#ff0000";
	      ps = BZ_CENTERPOINT_SIZE;
	    }

	    argc = 1;
	    sprintf(tag, "pt%d", ii);
	    WriteArgv( argv, &argc, "ssffff", "coords", tag,
		       vp_points[0] - ps, vp_points[1] - ps,
		       vp_points[0] + ps, vp_points[1] + ps);
	    (*info.proc)(info.clientData, interp, argc, (const char**) argv);

	    /* change the stack order */
	    argc = 1;
	    sprintf(tag, "pt%d", ii);
	    WriteArgv( argv, &argc, "sss", "raise", tag, "all");
	    (*info.proc)(info.clientData, interp, argc, (const char**) argv);
	  }
	}
      } else if ( iorder[i+1] > bz[typ].npoint && 
		  iorder[i+1] <= bz[typ].nZpoint ) {
	/* 
	 * it is a non EDGE-LINE seline 
	 */
	ii = iorder[i + 1] - bz[typ].npoint - 1;
	if ( !((bz[typ].slinetype[ii] & BZ_LINEPOINT) && 
	       (bz[typ].slinetype[ii] & BZ_EDGEPOINT)) ) {
	  argc = 1;
	  printf("ii = %d\n",ii);
	  fflush(stdout);
	  if (!bz[typ].sline_is_rendered[ii]) {
	    BzCreateSelLine( interp, info, argv, typ, ii );
	  } else {
	    sprintf(tag, "all");
	    BzUpdateSelLine( interp, info, argv, typ, ii, tag );
	  }
	}
      } else if ( iorder[i+1] > bz[typ].nZpoint ) {
	/*********************
	 * render RECVECTORS *
	 *********************/
	if ( bz_state[typ].vectors ) {
	  if ( !bz[typ].recvec_is_rendered ) {
	    char *arrow, *text;
	    if ( iorder[i+1] <= bz[typ].nZpoint+xcr.dim) {	
	      j = iorder[i+1] - (bz[typ].nZpoint + 1);
	      argc = 1;
	      /*********************************
	       * create 1st part of RECVECTORS *
	       *********************************/	     
	      WriteArgv( argv, &argc, "ssffff", "create", "line",
			 bz[typ].vp_recvec[j][0], bz[typ].vp_recvec[j][1],
			 bz[typ].vp_recvec[j+3][0], bz[typ].vp_recvec[j+3][1]);

	      ps = 2.0;
	      if ( bz[typ].rot_recvec[j][2] < 0.0 )
		ps = 1.0;	      
	      sprintf(tag, "vector v1_%d", j);
	      WriteArgv( argv, &argc, "sssfss", "-fill", "#0000ff", 
			 "-width", ps, "-tags", tag);
	      (*info.proc)(info.clientData, interp, argc, (const char**) argv);	    
	    } else {
	      j = iorder[i+1] - (bz[typ].nZpoint + xcr.dim + 1);
	      /* 2nd part of vector */
	      ps = 2.0;
	      arrow = "14 15 4";
	      if ( bz[typ].rot_recvec[j][2] < 0.0 ) {
		ps = 1.0;
		arrow = "10 11 3";
	      }
	      argc = 1;
	      WriteArgv( argv, &argc, "ssffff", "create", "line",
			 bz[typ].vp_recvec[j+3][0], bz[typ].vp_recvec[j+3][1],
			 bz[typ].vp_recvec[j+6][0], bz[typ].vp_recvec[j+6][1]);

	      sprintf(tag, "vector v2_%d", j);
	      WriteArgv( argv, &argc, "sssssssfss", "-fill", "#0000ff", 
			 "-arrow", "last", "-arrowshape", arrow,
			 "-width", ps, "-tags", tag);
	      (*info.proc)(info.clientData, interp, argc, (const char**) argv);

	      /*
	       * create label of recvectors
	       */
	      argc = 1;
	      WriteArgv( argv, &argc, "ssff", "create", "text",
			 bz[typ].vp_recvec[j+6][0] - 5,
			 bz[typ].vp_recvec[j+6][1]);
	      sprintf(tag, "vector t%d", j);
	      if (j == 0) text = "a*";
	      else if (j == 1) text = "b*";
	      else text = "c*";	/* if (j == 2) */
	      WriteArgv( argv, &argc, "ssssssss", "-fill", "#0000ff", 
			 "-anchor", "ne", "-text", text, "-tags", tag);
	      (*info.proc)(info.clientData, interp, argc, (const char**) argv);
	    }	    
	  } else {
	    /*********************
	     * update RECVECTORS *
	     *********************/
	    char *arrow;	    
	    if ( iorder[i+1] <= bz[typ].nZpoint+xcr.dim) {	
	      j = iorder[i+1] - (bz[typ].nZpoint + 1);
	      /* update 1st part of RECVECTORS */
	      argc = 1;
	      sprintf(tag, "v1_%d", j);
	      WriteArgv( argv, &argc, "ssffff", "coords", tag,
			 bz[typ].vp_recvec[j][0], bz[typ].vp_recvec[j][1],
			 bz[typ].vp_recvec[j+3][0], bz[typ].vp_recvec[j+3][1]);
	      (*info.proc)(info.clientData, interp, argc, (const char**) argv);
	      /* reconfigure recvectors */
	      argc = 1;
	      sprintf(tag, "v1_%d", j);
	      ps = 2.0;
	      if ( bz[typ].rot_recvec[j][2] < 0.0 )
		ps = 1.0;
	      WriteArgv( argv, &argc, "sssf", "itemconfigure", tag, 
			 "-width", ps); 
	      (*info.proc)(info.clientData, interp, argc, (const char**) argv);
	      /*
	       * get the right stack order for 1st part
	       */
	      argc = 1;
	      sprintf(tag, "v1_%d", j);
	      WriteArgv( argv, &argc, "sss", "raise", tag, "all" );
	      (*info.proc)(info.clientData, interp, argc, (const char**) argv);	       
	    } else {
	      j = iorder[i+1] - (bz[typ].nZpoint + xcr.dim + 1);
	      /* update 2nd part of RECVECTORS */
	      argc = 1;
	      sprintf(tag, "v2_%d", j);
	      WriteArgv( argv, &argc, "ssffff", "coords", tag,
			 bz[typ].vp_recvec[j+3][0], bz[typ].vp_recvec[j+3][1],
			 bz[typ].vp_recvec[j+6][0], bz[typ].vp_recvec[j+6][1]);
	      (*info.proc)(info.clientData, interp, argc, (const char**) argv);
	      /* reconfigure recvectors */
	      argc = 1;
	      ps = 2.0;
	      arrow = "14 15 4";
	      if ( bz[typ].rot_recvec[j][2] < 0.0 ) {
		ps = 1.0;
		arrow = "10 11 3";
	      }	      
	      argc = 1;
	      sprintf(tag, "v2_%d", j);
	      WriteArgv( argv, &argc, "sssssf", "itemconfigure", tag, 
			 "-arrowshape", arrow, "-width", ps);
	      (*info.proc)(info.clientData, interp, argc, (const char**) argv);

	      /*
	       * get the right stack order for 2nd part
	       */
	      argc = 1;
	      WriteArgv( argv, &argc, "sss", "raise", tag, "all" );
	      (*info.proc)(info.clientData, interp, argc, (const char**) argv);

	      /*
	       * update label of recvector
	       */
	      argc = 1;
	      sprintf(tag, "t%d", j);
	      WriteArgv( argv, &argc, "ssff", "coords", tag,
			 bz[typ].vp_recvec[j+6][0] - 5, 
			 bz[typ].vp_recvec[j+6][1]);
	      (*info.proc)(info.clientData, interp, argc, (const char**) argv);

	      /*
		printf("RECVEC:: iorder i=%d, ii=%d, v(z)=%f, p(z)=%f\n",
		i,ii,
		bz[typ].rot_recvec[i][2], bz[typ].rot_points[ii][2]);
	      */

	      /*
		printf("RECVEC:: iorder i=%d, ii=%d, v(z)=%f, p(z)=%f\n",
		i, ii,
		bz[typ].rot_recvec[i][2], bz[typ].rot_points[ii][2]);
		fflush(stdout);
	      */

	      argc = 1;
	      sprintf(tag, "t%d", j);
	      WriteArgv( argv, &argc, "sss", "raise", tag, "all" );
	      (*info.proc)(info.clientData, interp, argc, (const char**) argv);	      
	    }
	  }
	}
      }
    }
  }

  if ( bz_state[typ].vectors )
    bz[typ].recvec_is_rendered = 1;

  if ( !bz_state[typ].vectors && bz[typ].recvec_is_rendered ) {
    /*
     * delete recvectors
     */
    /* recvectors should be deleted */
    Tcl_VarEval(interp, argv[0], " delete", " vector", NULL);
    bz[typ].recvec_is_rendered = 0;    
  }

  /* render selected points */
  for (i=0; i<bz[typ].nselected; i++) {
    ii   = bz[typ].selectedID[i];
    argc = 1;
    if (!bz[typ].selected_is_rendered[i]) {
      WriteArgv( argv, &argc, "ss", "create", "oval");
      ps = BZ_BOLD_POINT_SIZE;
      if ( bz[typ].pointtype[ii] == BZ_CENTERPOINT )
	ps = BZ_CENTERPOINT_SIZE;

      color = "#00ff66";
      sprintf(tag, "selpoint spt%d", i);

      WriteArgv( argv, &argc, "ffffssssssss",
		 bz[typ].vp_points[ii][0] - ps, bz[typ].vp_points[ii][1] - ps,
		 bz[typ].vp_points[ii][0] + ps, bz[typ].vp_points[ii][1] + ps,
		 "-fill", color, "-outline", color, 
		 "-width", "1", "-tags", tag);	
      (*info.proc)(info.clientData, interp, argc, (const char**) argv);
      bz[typ].selected_is_rendered[i] = 1;
    } else {
      ps = BZ_BOLD_POINT_SIZE;
      if ( bz[typ].pointtype[ii] == BZ_CENTERPOINT )
	ps = BZ_CENTERPOINT_SIZE;
      sprintf(tag, "spt%d", i);
      WriteArgv( argv, &argc, "ssffff", "coords", tag, 
		 bz[typ].vp_points[ii][0] - ps,
		 bz[typ].vp_points[ii][1] - ps,
		 bz[typ].vp_points[ii][0] + ps,
		 bz[typ].vp_points[ii][1] + ps);
      (*info.proc)(info.clientData, interp, argc, (const char**) argv);

      /* change the stack order */
      argc = 1;
      WriteArgv( argv, &argc, "sss", "raise", tag, "all");
      (*info.proc)(info.clientData, interp, argc, (const char**) argv);
    }
  }

  /************************* 
   * render SELECTED LINES -> EDGE-LINE slines *
   *************************/
  for (i=0; i<bz[typ].n_sline; i++) {
    argc = 1;
    /* BZ_LINE must be rendered properly */
    if ( !bz[typ].sline_is_rendered[i] &&
	 ( (bz[typ].slinetype[i] & BZ_LINEPOINT) &&
	   (bz[typ].slinetype[i] & BZ_EDGEPOINT) ) ) {
      BzCreateSelLine( interp, info, argv, typ, i );
    } else if (bz[typ].sline_is_rendered[i] && 
	       ( (bz[typ].slinetype[i] & BZ_LINEPOINT) &&
		 (bz[typ].slinetype[i] & BZ_EDGEPOINT) ) ) {
      int pID, pID0, pID1;
      float sum0 = 0.0, sum1 = 0.0;

      pID0 = bz[typ].sline_index[i][0];
      pID1 = bz[typ].sline_index[i][1];
      /* if xcr.dim < 3, the second pID will be negative */
      if ( xcr.dim < 3 ) {
	pID = pID0;
      } else {
	for(j=0; j<bz[typ].nvert[pID0]; j++)
	  sum0 += bz[typ].rot_poly[pID0][j][2];
	for(j=0; j<bz[typ].nvert[pID1]; j++)
	  sum1 += bz[typ].rot_poly[pID1][j][2];

	sum0 /= (float) bz[typ].nvert[pID0];
	sum1 /= (float) bz[typ].nvert[pID1];

	pID = pID0;
	if ( sum1 > sum0 )
	  pID = pID1;    

	printf("LINE:: pID=%d, pID1=%d, pID2=%d; %f %f\n",
	       pID, pID0, pID1, sum0, sum1);
	fflush(stdout);
      }
      sprintf(tag, "p%d", pID);

      BzUpdateSelLine( interp, info, argv, typ, i, tag );
    }
  }

  bz[typ].isrendered = 1;
  bz_state[typ].last_points  = bz_state[typ].points;
  bz_state[typ].last_vectors = bz_state[typ].vectors;

  xcdebug ("end of BzRenderBZ #1");
  free((FREE_ARG) tag );
  for(i=0; i<BZRENDERBZ_ARGC; i++) 
    free((FREE_ARG) argv[i] );

  xcdebug ("end of BzRenderBZ #1");

  return TCL_OK;
}


static void
BzCreateSelLine( Tcl_Interp *interp, Tcl_CmdInfo info, 
		 char *argv[], int typ, int i)
{
  char tag[80], aps[3];
  int  argc = 1;

  sprintf(aps, "%d", BZ_BOLD_POLYGON_WIDTH);
  sprintf(tag, "selline sln%d", i);

  WriteArgv( argv, &argc, "ssffffssssssssss", "create", "line", 
	     bz[typ].vp_sline[i][0], bz[typ].vp_sline[i][1],
	     bz[typ].vp_sline[i][2], bz[typ].vp_sline[i][3],
	     "-fill", "#00ff66", "-arrow", "last",
	     "-arrowshape", "14 15 4", "-width", aps, "-tags", tag);
  (*info.proc)(info.clientData, interp, argc, (const char**) argv);

  bz[typ].sline_is_rendered[i] = 1;
}

static void
BzUpdateSelLine( Tcl_Interp *interp, Tcl_CmdInfo info, 
		 char *argv[], int typ, int i, char *tag )
{
  char tg[80];
  int  argc = 1;

  sprintf(tg, "sln%d", i);
  WriteArgv( argv, &argc, "ssffff", "coords", tg,
	     bz[typ].vp_sline[i][0], bz[typ].vp_sline[i][1],
	     bz[typ].vp_sline[i][2], bz[typ].vp_sline[i][3]);
  (*info.proc)(info.clientData, interp, argc, (const char**) argv);

  /* change the stack order */
  argc = 1;
  WriteArgv( argv, &argc, "sss", "raise", tg, tag);
  printf("LINE:: tag == %s\n", tag);
  fflush(stdout);
  (*info.proc)(info.clientData, interp, argc, (const char**) argv);
}


static int  
BzGetISS( int typ )
{
  int imen[3 * BZ_MAXPOINTS], imin = 0;
  double fabs();
  register int iss;
  register double f, point;
  register double tol = (double)(BZ_ISS - 1)/(double)BZ_ISS - (double)(BZ_ISS - 2)/(double)(BZ_ISS - 1);
  register int k, l, i, n;
  register double min;

  n = 3 * BZ_MAXPOINTS;
  for (i=0; i< n; i++)
    imen[i] = 0;

  /* t.k. begin: new algorithm due to a bug reported by P. Blaha */
  for (k=0; k<bz[typ].nselected; k++) 
    {
      n = 3 * k;

      for (l=0; l<3; l++) 
	{
	  min   = 1.0;
	  point = fabs( bz[typ].reccoor[ k ][ l ] );

	  fprintf(stderr, "k=%d, point=%f:  \n", k, point); 

	  if ( point < tol || point > (1.0 - tol) ) 
	    {
	      imen[n + l] = 1;
	    } 
	  else 
	    {
	      for (i=1; i<=BZ_ISS; i++) 
		{
		  f = point * (double) i;
		  f = fabs(f - roundf(f));
		  if ( f < min ) {
		    min  = f;
		    imin = i;
		  }

		  if ( f < tol ) {
		    imen[n + l] = i;
		    break;
		  }
		}

	      /* if imen[] was not found; assign imen[] := imin */
	      if ( imen[n + l] == 0 ) imen[n + l] = imin;      
	    }
	}
    }

  /* /\* t.k. begin: old corrected algorithm *\/ */
  /* for (k=0; k<bz[typ].nselected; k++) { */
  /*   n = 3 * k; */
  /*   for (l=0; l<3; l++) { */
  /*     min = 1.0; */
  /*     point = fabs( bz[typ].reccoor[ k ][ l ] ); */
  /*     if ( point < tol || point > (1.0 - tol) ) { */
  /* 	imen[n + l] = 1; */
  /* 	/\*break;*\/ */
  /*     } else { */
  /* 	for (i=1; i<=BZ_ISS; i++) */
  /* 	  for (j=1; j<i; j++) { */
  /* 	    f = fabs( point - (double) j / (double) i); */
  /* 	    if ( f < min ) { */
  /* 	      min  = f; */
  /* 	      imin = i; */
  /* 	    } */
  /* 	    if ( f < tol ) { */
  /* 	      imen[n + l] = i; */
  /* 	      i = BZ_ISS + 1; /\* to stop also the second loop *\/ */
  /* 	      break; */
  /* 	    } */
  /* 	  } */
  /* 	/\* if imen[] was not found; assign imen[] := imin *\/ */
  /* 	if ( imen[n + l] == 0 ) imen[n + l] = imin;       */
  /*     } */
  /*   } */
  /* } */
  /* /\* t.k. end *\/ */

  /* now find the minimum common multiplier */
  iss = imen[0];

  n = 3 * bz[typ].nselected;
  for (i=1; i<n; i++)
    {
      iss = BzGetMinMult( iss, imen[i] );
    }

  return iss;
}


static int
BzGetMinMult( int a, int b )
{  
  register int i;
  int c;

  /* maybe both integers are ZERO */
  if ( a == 0 && b == 0 ) {
    return 1;
  }
  /* suppose that: a>b */  
  if ( b > a ) {
    c = a;
    a = b;
    b = c;
  }
  if ( b == 0 ) {
    return a;
  }

  c=1;
  /* now: a is greater than b */
  for (i=1; i<=b; i++)
    if ( (a * i % b) == 0 ) {
      c = i;
      break;
    }

  return a * c;
}
