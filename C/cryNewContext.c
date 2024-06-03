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
 * Source: $XCRYSDEN_TOPDIR/C/cryNewContext.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <togl.h>
#include <stdio.h> 
#include <stdlib.h> 
#include "struct.h"

NEW_WIN_CONTEXT *wc=NULL, *wcPtr=NULL;

NEW_WIN_CONTEXT *FindWinContextByTogl(struct Togl *togl) {
  NEW_WIN_CONTEXT *g = wcPtr;
  
  while (g) {
    if ( togl == g->togl )
      return g;
    g = g->ptr;
  }
  return NULL;
}

NEW_WIN_CONTEXT *FindWinContext(int index)
{
   NEW_WIN_CONTEXT *g = wcPtr;
   while (g) {
     if (index == g->index)
       return g;
     g = g->ptr;
   }
   return NULL;
}

void DestroyWinContext(struct Togl *togl)
{
  NEW_WIN_CONTEXT *prev = NULL;
  NEW_WIN_CONTEXT *pos = wcPtr;
  NEW_WIN_CONTEXT *g = FindWinContextByTogl( togl );
  
  while (pos) {
    if (pos == g ) {
      if (prev) {
	prev->ptr = pos->ptr;
      }
      else {
	wcPtr = pos->ptr;
      }
      return;
    }
    prev = pos;
    pos = pos->ptr;
  }
  free( (FREE_ARG) g );
}

void AddToWinContext(NEW_WIN_CONTEXT *g)
{
  g->ptr = wcPtr;
  wcPtr  = g;
}

int NewWinContext(void) 
{
  static int index = -1, i;
  wc = (NEW_WIN_CONTEXT *) malloc ( sizeof(NEW_WIN_CONTEXT) );
  wc->index=++index;
  /*---------------*/
  wc->VPf.stropened  = 0;
  wc->VPf.nsurface   = 0;
  wc->VPf.surface    = NULL; 
  wc->VPf.surfacePtr = NULL;
  wc->VPf.surfaceInd = NULL;
  wc->VPf.antialias  = GL_FALSE;
  wc->VPf.fog        = GL_FALSE;
    
  wc->tr.xtransl   = 0.0;
  wc->tr.ytransl   = 0.0;
  wc->tr.rotx      = 0.0;
  wc->tr.roty      = 0.0;
  wc->tr.zoom      = DEF_ZOOM;
  /* wc->tr.zoom = 1.0; */
  wc->tr.b1motion  = 0;
  wc->tr.b2motion  = 0;
  wc->recprim_cage = (void*)NULL;

  wc->fermiContext.toglVector = NULL;

  for(i=0; i<4; i++)
    wc->bg[i] = DefBg[i];
  AddToWinContext( wc );
  return index;
}

void
cryNewToglInit(struct Togl *togl) 
{
  wc = FindWinContext( NewWinContext() );
  wc->togl = togl;
}
  
