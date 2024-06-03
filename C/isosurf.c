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
 * Source: $XCRYSDEN_TOPDIR/C/isosurf.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <stdio.h> 
#include <stdlib.h> 
#include "isosurf.h"

ISOSURFACE *iso=NULL, *isoPtr=NULL;

ISOSURFACE *FindIsoSurf(int index)
{
   ISOSURFACE *g = isoPtr;
   while (g) {
     if (index == g->index)
       return g;
     g = g->ptr;
   }
   return NULL;
}
void AddToIsoSurfList(ISOSURFACE *g)
{
  g->ptr  = isoPtr;
  isoPtr  = g;
}
int NewIsoSurf(void) 
{
  static int index = -1;
  iso = (ISOSURFACE *) malloc ( sizeof(ISOSURFACE) );
  iso->index=++index;
  /*---------------*/
  iso->triangl_status = NULL;
  iso->vertex         = NULL;
  iso->vertex_orig    = NULL;
  iso->normal         = NULL;
  iso->color          = NULL;
  iso->tri2verIN      = NULL;
  iso->nver2triIN     = NULL;
  iso->ver2triIN      = NULL;
  AddToIsoSurfList( iso );
  return index;
}
