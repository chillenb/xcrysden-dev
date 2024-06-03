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
* Source: $XCRYSDEN_TOPDIR/C/gl2psPrintTogl.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <stdio.h>
#include <string.h>
#include <togl.h>
#include "gl2ps.h"
#include "xcfunc.h"

extern struct Togl *mesa_togl;
extern void (*xcDisplay)(struct Togl *togl); 

/************************************************************************

Tcl-usage: cry_gl2psPrintTogl $togl $format $sort $options $pointsize $linewidth $stream $filename

Usage: 

CRY_gl2psPrintToglCb format sort options pointsize linewidth filename


************************************************************************/
int 
CRY_gl2psPrintToglCb(ClientData clientData, Tcl_Interp *interp,
		     int argc, const char *argv[])
{
  Togl *togl;
  FILE  *fp;
  int   argcList, i;
  const char  **argvList;
  GLint format, sort, options = GL2PS_NONE;
  GLint buffersize;
  GLint state = GL2PS_OVERFLOW;
  GLfloat pointsize=2.0, linewidth=2.0;
  GLint viewport[4];
  double num;

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if (argc != 8) {
    Tcl_SetResult(interp,
		  "wrong # args: should be \"cry_gl2psPrintTogl toglname format sort options pointsize linewidth filename\"",
		  TCL_STATIC);
    return TCL_ERROR;
  }

  buffersize = Togl_Width (togl) * Togl_Height (togl) * 8;

  /*
    FORMAT
  */
  if ( strcmp(argv[2], "GL2PS_PS") == 0 ) {
    format = GL2PS_PS;
  } else if ( strcmp(argv[2], "GL2PS_EPS") == 0 ) {
    format = GL2PS_EPS;
  } else if ( strcmp(argv[2], "GL2PS_TEX") == 0 ) {
    format = GL2PS_TEX;
  } else if ( strcmp(argv[2], "GL2PS_PDF") == 0 ) {
    format = GL2PS_PDF;
  } else if ( strcmp(argv[2], "GL2PS_SVG") == 0 ) {
    format = GL2PS_SVG;
  } else {
    char rss[1024];
    snprintf(rss, sizeof(rss), "wrong format %s, must be one ofGL2PS_PS, GL2PS_EPS, or GL2PS_TEX", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /*
    SORT
  */
  if ( strcmp(argv[3], "GL2PS_NO_SORT") == 0 ) {
    sort = GL2PS_NO_SORT;
  } else if ( strcmp(argv[3], "GL2PS_SIMPLE_SORT") == 0 ) {
    sort = GL2PS_SIMPLE_SORT;
  } else if ( strcmp(argv[3], "GL2PS_BSP_SORT") == 0 ) {
    sort = GL2PS_BSP_SORT;
  } else {
    char rss[1024];
    snprintf(rss, sizeof(rss), "wrong sort %s, must be one of GL2PS_NO_SORT, GL2PS_SIMPLE_SORT, or GL2PS_BSP_SORT", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /*
    OPTIONS
  */
  Tcl_SplitList(interp, argv[4], &argcList, &argvList);
  for (i=0; i<argcList; i++) {
    if ( strcmp(argvList[i], "GL2PS_DRAW_BACKGROUND") == 0 ) {
      options |= GL2PS_DRAW_BACKGROUND;
    } else if ( strcmp(argvList[i], "GL2PS_SIMPLE_LINE_OFFSET") == 0 ) {
      options |= GL2PS_SIMPLE_LINE_OFFSET;
    } else if ( strcmp(argvList[i], "GL2PS_SILENT") == 0 ) {
      options |= GL2PS_SILENT;
    } else if ( strcmp(argvList[i], "GL2PS_BEST_ROOT") == 0 ) {
      options |= GL2PS_BEST_ROOT;
    } else if ( strcmp(argvList[i], "GL2PS_OCCLUSION_CULL") == 0 ) {
      options |= GL2PS_OCCLUSION_CULL;
    } else if ( strcmp(argvList[i], "GL2PS_NO_TEXT") == 0 ) {
      options |= GL2PS_NO_TEXT;
    } else if ( strcmp(argvList[i], "GL2PS_LANDSCAPE") == 0 ) {
      options |= GL2PS_LANDSCAPE;
    } else if ( strcmp(argvList[i], "GL2PS_NO_PS3_SHADING") == 0 ) {
      options |= GL2PS_NO_PS3_SHADING;
    } else if ( strcmp(argvList[i], "GL2PS_NO_PIXMAP") == 0 ) {
      options |= GL2PS_NO_PIXMAP;
    } else if ( strcmp(argvList[i], "GL2PS_NO_BLENDING") == 0 ) {
      options |= GL2PS_NO_BLENDING;
    } else if ( strcmp(argvList[i], "GL2PS_NONE") != 0 ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "wrong options %s, must be one of GL2PS_DRAW_BACKGROUND, GL2PS_SIMPLE_LINE_OFFSET, GL2PS_SILENT, GL2PS_BEST_ROOT, GL2PS_OCCLUSION_CULL, GL2PS_NO_TEXT, GL2PS_LANDSCAPE, GL2PS_NO_PS3_SHADING, or GL2PS_NONE", argvList[i]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  /*
    PointSize
  */
  if ( Tcl_GetDouble(interp, argv[5], &num) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "wanted double but got %s", argv[5]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  } else {
    pointsize = (GLfloat) num;    
  }

  /*
    LineWidth
  */
  if ( Tcl_GetDouble(interp, argv[6], &num) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "wanted double but got %s", argv[6]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  } else {
    linewidth = (GLfloat) num;
  }

  if ( format == GL2PS_PDF ) {
    fp = fopen(argv[7], "wb");
  } else {
    fp = fopen(argv[7], "w");
  }

  if ( fp == NULL ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "couldn't open file %s", argv[7]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  fprintf(stderr, "\nCRY_gl2psPrintToglCb: printing to %s\n\n", argv[2]);

  glGetIntegerv(GL_VIEWPORT,viewport);

  while ( state == GL2PS_OVERFLOW ) {
    gl2psBeginPage ("XCRYSDEN image", "XCRYSDEN 1.x via GL2PS",
		    viewport,
		    format, sort, options,
		    GL_RGBA, 0, NULL, 
		    8, 8, 8,
		    buffersize, fp, argv[2]);  
    gl2psPointSize (pointsize);
    gl2psLineWidth (linewidth);  

    /**/
    if ( togl == mesa_togl ) 
      {
	(*xcDisplay)(togl);  
      }
    else 
      {
	NEW_WIN_CONTEXT *wc = FindWinContextByTogl( togl );
	(*wc->xcDisplay)(togl);
      }
    /**/

    state = gl2psEndPage();
    fprintf(stderr,"CRY_gl2psPrintToglCb$ state = %d\n", state);
    buffersize *= 2;
  }
  fclose(fp);

  return TCL_OK;
}

