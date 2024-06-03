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
* Source: $XCRYSDEN_TOPDIR/C/ppmPrintTogl.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <togl.h>
#include "struct.h"
#include "xcfunc.h"
#include "memory.h"

extern void (*xcDisplay)(struct Togl *togl);
extern struct Togl *mesa_togl;

realTimeMove makeMovie = { MOVIE_MODE_EVERY_SNAPSHOT, 0, 0, 0, (char*)NULL };
static char *filelist = (char*)NULL;
static int   filelist_size = 20000;

/*
 * this function takes care of real-time screen PPM capture for movie creation
 * 
 * Usage: 
 * 
 * xc_realtimemovie toglName begin mode dir
 * xc_realtimemovie toglName end
 * xc_realtimemovie toglName filelist
 * xc_realtimemovie toglName clear
 */
int 
CRY_RealTimeMovieCb(ClientData clientData, Tcl_Interp *interp,
		    int argc, const char *argv[])
{
  Togl *togl;

  if ( argc < 3 && argc > 5) {
    Tcl_SetResult(interp, "Usage: \"realtimemovie toglName begin mode dir\" or \"realtimemovie toglName end\"", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( togl != mesa_togl ) {
    Tcl_SetResult(interp, "so far toglName realtimemovie callback works only for .mesa togl", TCL_STATIC);
    return TCL_ERROR;
  }

  if (strncmp(argv[2],"beg",3) == 0) {
    if ( argc != 5) {
      Tcl_SetResult(interp, "Usage: toglName realtimemovie begin mode dir", TCL_STATIC);
      return TCL_ERROR;
    }
    /* BEGIN */
    if ( makeMovie.doit == 1 ) {
      Tcl_SetResult(interp, "can't call \"toglName realtimemovie begin\" twice ...", TCL_STATIC);
      return TCL_ERROR;
    }

    if ( filelist != NULL) {
      filelist[0] = '\0';
    }

    makeMovie.doit = 1;   
    makeMovie.nframe = 0;  

    /* mode */

    makeMovie.mode = MOVIE_MODE_EVERY_SNAPSHOT;
    if (strncmp(argv[3],"real",4) == 0) {
      makeMovie.mode = MOVIE_MODE_REALTIME_INTERVAL;
    }

    /* directory where to write files */
    if ( makeMovie.dir != NULL ) {
      Tcl_SetResult(interp, "\"toglName realtimemovie begin\" re-called without clearing ...", TCL_STATIC);
      return TCL_ERROR;
    }
    makeMovie.dir = (char*) xcCalloc( strlen(argv[4]) + 1, sizeof(char) );
    strcpy(makeMovie.dir, argv[4]);
  } 

  else if (strncmp(argv[2],"end",3) == 0) {
    /* END */
    if ( makeMovie.doit != 1 ) {
      Tcl_SetResult(interp, "can't call \"toglName realtimemovie end\" prior to \"toglName realtimemovie begin\"", TCL_STATIC);
      return TCL_ERROR;
    }

    xcFree(makeMovie.dir); 
    makeMovie.dir  = NULL;    
    makeMovie.doit = 0;    
  }

  else if (strncmp(argv[2],"filelist",4) == 0) 
    {
      if ( filelist != NULL) {
	/*char *result = Tcl_Alloc( sizeof(char) * strlen(filelist+1) );*/
	fprintf(stderr, "size: %ld\nfilelist: %s\n", strlen(filelist), filelist);
	Tcl_SetResult(interp,filelist,TCL_VOLATILE);    
      } else {
	Tcl_SetResult(interp,"",TCL_STATIC);
      }
    }

  else if (strncmp(argv[2],"clear",4) == 0) 
    {
      xcFree(makeMovie.dir); 
      makeMovie.dir    = NULL;    
      makeMovie.doit   = 0;  
      makeMovie.nframe = 0; 
      if ( filelist != NULL) {
	filelist[0] = '\0';
      }
    }

  return TCL_OK;
}

void
createMoviePPMFrame(struct Togl *togl)
{
  char *filename;

  if (makeMovie.dir != NULL) {
    filename = (char*) xcCalloc(strlen(makeMovie.dir) + 50, sizeof(char));
  } else {
    return;
  }

  fprintf(stderr,"makeMovie.printing=%d\n",makeMovie.printing);
  if (makeMovie.printing) return;

  makeMovie.printing = 1;

  if ( togl == mesa_togl ) {
    /* so far it works only for ".mesa" */

    /* take care of filelist */
    if ( filelist == NULL ) {    
      filelist = xcCalloc(filelist_size, sizeof(char));
    }
    sprintf(filename,"%s/frame-%05d.ppm", makeMovie.dir, makeMovie.nframe++);

    if (strlen(filelist) + strlen(filename) + 1 >= filelist_size) 
      {
	filelist_size *= 2;
	filelist       = xcRealloc(filelist, sizeof(char)*filelist_size);
      }
    strcat(filelist,filename);
    strcat(filelist," ");

    /* dump to PPM file */

    fprintf(stderr,"createMoviePPMFrame: making %s ... ", filename);

    Togl_DumpToPpmFile( togl, filename );
    fprintf(stderr,"ok\n");
  }

  makeMovie.printing = 0;

  xcFree((FREE_ARG)filename);
}
