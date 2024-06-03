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
* Source: $XCRYSDEN_TOPDIR/C/xcWigner.c
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
#include "xcGLparam.h"
#include "wigner.h"
#include "memory.h"
#include "xcfunc.h"

extern Options3D is;

int wigner_nodes_malloc[2] = {0, 0};


static int ReadWSNodes( int mode );
static float DetWignerSize(void);

/*****************************************************************************
xc_wigner  <toglName> prim|conv -nodesfile <file>
-transparency 0|1
-color {r g b}
-render after|now (DEFAULT:: after)
<toglName> clear			     
*****************************************************************************/
int 
XC_WignerCmd(ClientData clientData,Tcl_Interp *interp,
	     int argc, const char *argv[])
{
  struct Togl *togl;
  register int i;
  static int first_time = 1;
  int mode, trans;
  int nodemode = WIGNER_NODEMODE_EVERY;
  int render = 0;
  GetGlParam color = { 0, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			   0, 0, 0, 0, 0, 0, 0, 0, 0, 0} };

  if (first_time) {
    LoadVoronoiMaterial( VORONOI_WIGNERSEITZ );
    first_time = 1;
  }

  if ( argc < 3 ) {
    Tcl_SetResult(interp, "Usage: xc_wigner  <toglName> prim|conv -nodesfile <file> -transparency 0|1 -color {r g b} -render now|after\nor\nxc_wigner <toglName> clear", TCL_STATIC);
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

  if ( strcmp(argv[2], "prim") == 0 ) {
    mode = WIGNER_PRIM;
  }
  else if ( strcmp(argv[2], "conv") == 0 ) {
    mode = WIGNER_CONV;
  }
  else if ( strcmp(argv[2], "clear") == 0 ) {
    VPf.wignerseitz = GL_FALSE;
    mode = WIGNER_CLEAR;
    ws_attrib.is_initialized[0] = 0;
    ws_attrib.is_initialized[1] = 0;
    Togl_PostRedisplay(togl);
    return TCL_OK;
  } 
  else {
    char rss[1024];
    snprintf(rss, sizeof(rss), "unknown option \"%s\", must be one of prim, conv, clear", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( mode == WIGNER_PRIM || mode == WIGNER_CONV ) {
    if ( !xcr.lprimvec ) {
      Tcl_SetResult(interp, "can not render Wigner-Seitz cell, since primitive cell vectors were not specified", TCL_STATIC);
      return TCL_ERROR;
    }
    for (i=3; i<argc; i+=2) {
      if ( strncmp(argv[i], "-nodesf", 7) == 0 ) {
	if ( (ws_attrib.fp = fopen(argv[i+1],"r")) == NULL ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"can't open file \"%s\"",argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  fclose( ws_attrib.fp );
	  return TCL_ERROR;
	}
	if ( !ReadWSNodes( mode ) ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "error reading %s file !!!", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  fclose( ws_attrib.fp );
	  return TCL_ERROR;
	}
	fclose( ws_attrib.fp );
	nodemode = WIGNER_NODEMODE_SELECT;
      }
      else if ( strncmp(argv[i], "-transp", 7) == 0 ) {
	if ( Tcl_GetInt( interp, argv[i+1], &trans ) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "wanted integer, but got %s", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	if ( trans < 0 || trans > 1 ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "invalid value %d for -transparency value, should be 0 or 1", trans);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	ws_attrib.transparency = trans;
      }
      else if ( strcmp(argv[i], "-color") == 0 ) {
	if ( !xcSplitList(XC_GET_XYZ, interp, &argv[i+1], &color) ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "parse error; when parsing -color {%s}, should be -color {r g b}", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	ws_attrib.color[0] = color.vec[0];
	ws_attrib.color[1] = color.vec[1];
	ws_attrib.color[2] = color.vec[2];
	/* so far we can submit just RGB and not RGBA, so I will use default
	   ALFA */
	ws_attrib.color[3] = mat_voronoi[0].ambient[3];
      }
      else if ( strcmp(argv[i], "-render") == 0 ) {
	if ( strcmp(argv[i+1], "now" ) == 0 ) render = 1;
	else if ( strcmp(argv[i+1], "after") == 0 ) render = 0;
	else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "unknown value \"%s\" for -render option, must be now or after", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss), "unknown option \"%s\", must be one of -nodesfile, -transparency, -color or -render", argv[i]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    ws_attrib.is_initialized[mode] = 1;
    ws_attrib.nodemode[mode] = nodemode;
    VPf.wignerseitz = GL_TRUE;
    /* determine the size of wigner cells */
    if ( render ) {
      MVf.wignerseitzsize = DetWignerSize();
      /* because Wigner_Seitz cell may be greater than structure, 
	 take care of projection */
      if (is.stickmode && !is.ballmode) xcMakeProjection3D("sticks");
      if (is.ballmode) xcMakeProjection3D("balls");
      if (is.spacefillmode) xcMakeProjection3D("space");
    }
  }

  if (render) Togl_PostRedisplay(togl);

  return TCL_OK;
}


static float
DetWignerSize(void)
{
  int mode = WIGNER_PRIM;
  register int i, j;
  struct Node {
    float min_cellind[3];
    float max_cellind[3];
  } nd;
  double (*vector)[4]; /* because vec.XXXX are (double [4][4]) */
  float vm[3][3], vM[3][3], orig[3];

  if ( xcr.celltype == XCR_CONVCELL ) {
    mode = WIGNER_CONV;
  }

  if ( ws_attrib.nodemode[mode] == WIGNER_NODEMODE_EVERY ) {
    /* min = cell(0,0,0); 
       max = cell(xcr.nunit[0], xcr.nunit[1], xcr.nunit[2] */
    nd.min_cellind[0] = nd.min_cellind[1] = nd.min_cellind[2] = 0.0;
    nd.max_cellind[0] = (float) xcr.nunit[0];
    nd.max_cellind[1] = (float) xcr.nunit[1];
    nd.max_cellind[2] = (float) xcr.nunit[2];    
  } 
  else if ( ws_attrib.nodemode[mode] == WIGNER_NODEMODE_SELECT ) {
    nd.min_cellind[0] = (float) xcr.nunit[0];
    nd.min_cellind[1] = (float) xcr.nunit[1];
    nd.min_cellind[2] = (float) xcr.nunit[2];    
    nd.max_cellind[0] = nd.max_cellind[1] = nd.max_cellind[2] = 0.0;
    for (i=0; i<ws_attrib.nnodes[mode]; i++) {
      for (j=0; j<3; j++) {
	if ( ws_attrib.nodes[mode][i][j] < nd.min_cellind[j] )
	  nd.min_cellind[j] = ws_attrib.nodes[mode][i][j];
	if ( ws_attrib.nodes[mode][i][j] > nd.max_cellind[j] )
	  nd.max_cellind[j] = ws_attrib.nodes[mode][i][j];
      }
    }
  }

  /* 
     determine size of single WS cell !!!
     Size of WS cell can not be greater then size of primcell, so
     use the size of primcell for the determination of WS cell
     for (i=0; i<3; i++)
     orig[i] = ( vec.prim[0][i] + vec.prim[1][i] + vec.prim[2][i] ) / 3.0;
     ws_size = DetermineParapipedSize( vec.prim[0], vec.prim[1], 
     vec.prim[2], orig );
  */

  vector = vec.prim;
  if ( xcr.celltype == XCR_CONVCELL && 
       ws_attrib.nodemode[mode] == WIGNER_NODEMODE_EVERY ) {
    vector = vec.conv;
  }
  for (i=0; i<3; i++) {
    orig[i] = 0.0;
    for (j=0; j<3; j++) {
      vm[i][j] = (nd.min_cellind[i] - 0.5) * vector[i][j];
      vM[i][j] = (nd.max_cellind[i] + 0.5) * vector[i][j];
      orig[i] -= vm[i][j];
    }
  }
  orig[0] += mx;
  orig[1] += my;
  orig[2] += mz;

  return DetermineParapipedSize(vM[0], vM[1], vM[2], orig);
}


static int
ReadWSNodes( int mode )
{
  int nnodes;
  register int i;

  if ( fscanf(ws_attrib.fp,"%d", &nnodes) != 1 )  return XC_ERROR;

  if ( wigner_nodes_malloc[mode] ) {
    xcFree_Matrixf( (float **) ws_attrib.nodes[mode] );
    wigner_nodes_malloc[mode] = 0;
  }

  if ( !wigner_nodes_malloc[mode] ) {
    ws_attrib.nodes[mode] = xcMallocMatrixf( nnodes, 3 );
    wigner_nodes_malloc[mode] = 1;
  }

  for (i=0; i<nnodes; i++)
    if ( fscanf(ws_attrib.fp, "%f %f %f", 
		&ws_attrib.nodes[mode][i][0],
		&ws_attrib.nodes[mode][i][1],
		&ws_attrib.nodes[mode][i][2] ) != 3 ) {
      ws_attrib.nnodes[mode] = 0;
      return XC_ERROR;
    }

  ws_attrib.nnodes[mode] = nnodes;

  return XC_OK;
}
