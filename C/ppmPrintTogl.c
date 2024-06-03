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

extern struct Togl *mesa_togl;
extern const char *printImage;
static int FileWritePPM(Tcl_Interp *interp, const char *fileName, Tk_PhotoImageBlock *blockPtr);

/* --- xcDisplayFunc.c --- */
extern void (*xcDisplay)(struct Togl *togl);

/*
 * this function takes care of PPM
 * 
 * Usage: xc_dump2ppm toglName filename
 */
int
CRY_Dump2PpmCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
  Togl *togl;

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( argc != 3 ) {
    Tcl_SetResult(interp, "Usage: xc_dump2ppm toglName filename", TCL_STATIC);
    return TCL_ERROR;
  }

  /* /\* check if togl is .mesa or something else ... *\/ */
  /* if ( togl == mesa_togl ) { */
  /*   Togl_DumpToPpmFile( togl, argv[2] ); */
  /* } else { */
  /*   /\*NEW_WIN_CONTEXT *wc; */
  /*     wc = FindWinContextByTogl( togl );*\/ */
  /*   Togl_DumpToPpmFile( togl, argv[2] ); */
  /* } */

  Togl_DumpToPpmFile( togl, argv[2] );
  
  return TCL_OK;
}


int Togl_DumpToPpmFile(Togl *togl, const char *filename)  
{
  Tcl_Interp *interp = Togl_Interp(togl);
  /*ClientData clientData = Togl_GetClientData(togl);*/
  Tk_PhotoHandle photo;
  Tk_PhotoImageBlock blockPtr;
    
  photo = Tk_FindPhoto(interp, printImage);

  if (photo == NULL) {
    Tcl_AppendResult(interp, "image \"", printImage,
                     "\" doesn't exist or is not a photo image", NULL);
    return TCL_ERROR;
  }
  
  Togl_TakePhoto(togl, photo);  
  Tk_PhotoGetImage(photo, &blockPtr);
  FileWritePPM(interp, filename, &blockPtr);
  
  return TCL_OK;
}


/************************************************************************* 

   The below "FileWritePPM" routine is taken from tk8.6 sources, in
   particular from file "tk8.6.9/generic/tkImgPPM.c"

*************************************************************************/

/*
 *----------------------------------------------------------------------
 *
 * FileWritePPM --
 *
 *	This function is invoked to write image data to a file in PPM format
 *	(although we can read PGM files, we never write them).
 *
 * Results:
 *	A standard TCL completion code. If TCL_ERROR is returned then an error
 *	message is left in the interp's result.
 *
 * Side effects:
 *	Data is written to the file given by "fileName".
 *
 * Copyright (c) 1994 The Australian National University.
 * Copyright (c) 1994-1997 Sun Microsystems, Inc.
 *
 * See the file "$XCRYSDEN_TOPDIR/otherLICENSES/TclTk:LICENSE" for
 * information on usage and redistribution of this file, and for a
 * DISCLAIMER OF ALL WARRANTIES.
 *
 * Author: Paul Mackerras (paulus@cs.anu.edu.au),
 *      Department of Computer Science,
 *      Australian National University.
 */

static int
FileWritePPM(Tcl_Interp *interp,
             const char *fileName,
             Tk_PhotoImageBlock *blockPtr)
{
  Tcl_Channel chan;
  int w, h, greenOffset, blueOffset, nBytes;
  unsigned char *pixelPtr, *pixLinePtr;
  char header[16 + TCL_INTEGER_SPACE * 2];

  chan = Tcl_OpenFileChannel(interp, fileName, "w", 0666);
  if (chan == NULL) {
    return TCL_ERROR;
  }

  if (Tcl_SetChannelOption(interp, chan, "-translation", "binary")
      != TCL_OK) {
    Tcl_Close(NULL, chan);
    return TCL_ERROR;
  }
  if (Tcl_SetChannelOption(interp, chan, "-encoding", "binary")
      != TCL_OK) {
    Tcl_Close(NULL, chan);
    return TCL_ERROR;
  }

  sprintf(header, "P6\n%d %d\n255\n", blockPtr->width, blockPtr->height);
  Tcl_Write(chan, header, -1);

  pixLinePtr = blockPtr->pixelPtr + blockPtr->offset[0];
  greenOffset = blockPtr->offset[1] - blockPtr->offset[0];
  blueOffset = blockPtr->offset[2] - blockPtr->offset[0];

  if ((greenOffset == 1) && (blueOffset == 2) && (blockPtr->pixelSize == 3)
      && (blockPtr->pitch == (blockPtr->width * 3))) {
    nBytes = blockPtr->height * blockPtr->pitch;
    if (Tcl_Write(chan, (char *) pixLinePtr, nBytes) != nBytes) {
      goto writeerror;
    }
  } else {
    for (h = blockPtr->height; h > 0; h--) {
      pixelPtr = pixLinePtr;
      for (w = blockPtr->width; w > 0; w--) {
        if (    Tcl_Write(chan,(char *)&pixelPtr[0], 1) == -1 ||
                Tcl_Write(chan,(char *)&pixelPtr[greenOffset],1)==-1 ||
                Tcl_Write(chan,(char *)&pixelPtr[blueOffset],1) ==-1) {
          goto writeerror;
        }
        pixelPtr += blockPtr->pixelSize;
      }
      pixLinePtr += blockPtr->pitch;
    }
  }

  if (Tcl_Close(NULL, chan) == 0) {
    return TCL_OK;
  }
  chan = NULL;

 writeerror:
  Tcl_SetObjResult(interp, Tcl_ObjPrintf("error writing \"%s\": %s",
                                         fileName, Tcl_PosixError(interp)));
  if (chan != NULL) {
    Tcl_Close(NULL, chan);
  }
  return TCL_ERROR;
}
