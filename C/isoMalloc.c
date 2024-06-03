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
 * Source: $XCRYSDEN_TOPDIR/C/isoMalloc.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "isosurf.h"
#include "memory.h"


void
xcFree_PLANEVERTEX(PLANEVERTEX **m)
{
  free( (FREE_ARG) m[0] );
  free( (FREE_ARG) m );
}


void
xcFree_GRIDVERTEX(GRIDVERTEX ***t)
{
  free( (FREE_ARG) t[0][0] );
  free( (FREE_ARG) t[0] );
  free( (FREE_ARG) t );
}


void
xcFree_LINE(LINE *t)
{
  free( (FREE_ARG) t );
}


PLANEVERTEX **xcMallocPLANEVERTEX(long nr, long nc)
{
  long i;
  PLANEVERTEX **m;
  /* allocate pointers to pointers to rows */
  m=(PLANEVERTEX **) malloc( (size_t) nr*sizeof(PLANEVERTEX *) );
  if (!m) xcError("allocation failure 1 in xcMallocPLANEVERTEX");

  /* allocate pointers to rows and set pointers to them */
  m[0]=(PLANEVERTEX *) malloc( (size_t) nr*nc*sizeof(PLANEVERTEX) );
  if (!m[0]) xcError("allocation failure 2 in xcMallocPLANEVERTEX");

  for(i=1; i<nr; i++)
    m[i] = m[i-1] + nc;
  
  return m;
}


GRIDVERTEX ***xcMallocGRIDVERTEX(long nr, long nc, long nd)
{  
  /* allocate a float 3tensor with range t[0..nr-1][0..nc-1][0..nd-1] */
  long i, j;
  GRIDVERTEX ***t;
  /* allocate pointers to pointers to rows */
  t=(GRIDVERTEX ***) malloc( (size_t) nr*sizeof(GRIDVERTEX **) );
  if (!t) xcError("allocation failure 1 in xcMallocGRIDVERTEX");

  /* allocate pointers to rows and set pointers to them */
  t[0]=(GRIDVERTEX **) malloc( (size_t) nr*nc*sizeof(GRIDVERTEX *) );
  if (!t[0]) xcError("allocation failure 2 in xcMallocGRIDVERTEX");

  /* allocate rows and set pointers to them */
  t[0][0]=(GRIDVERTEX *) malloc( (size_t) nr*nc*nd*sizeof(GRIDVERTEX) );
  if (!t[0][0]) xcError("allocation failure 3 in xcMallocGRIDVERTEX");
  
  for(i=0; i<nr; i++) {
    t[i] = t[0] + i*nc;
    for(j=0; j<nc; j++)
      t[i][j] = t[0][0] + nd*(i*nc + j);
  }
  return t;
}


LINE *xcMallocLINE(long nr)
{
  LINE *v;
  /* allocate pointers to pointers to rows */
  v=(LINE *) malloc( (size_t) nr*sizeof(LINE) );
  if (!v) xcError("allocation failure in xcMallocLINE");
  return v;
}


VERTEX *xcRealloc2xBiggerVertex(VERTEX *vec, int *size)
{
  register int i;
  VERTEX *oldvec = vec;

  vec = (VERTEX *) malloc( (size_t) sizeof(VERTEX) * 2 * *size);
  if (!vec) xcError("allocation failure in xcRealloc2xBiggerVertex");
  for(i=0; i<*size; i++) {
    vec[i].pos.x = oldvec[i].pos.x;
    vec[i].pos.y = oldvec[i].pos.y;
    vec[i].pos.z = oldvec[i].pos.z;
    vec[i].nml.x = oldvec[i].nml.x;
    vec[i].nml.y = oldvec[i].nml.y;
    vec[i].nml.z = oldvec[i].nml.z; /* tole bi moral biti constructor, pa verjetno ni */
  }
  
  free( (FREE_ARG) oldvec );
  *size *=2;

  return vec;
}


TRIANGLE *xcRealloc2xBiggerTriangl(TRIANGLE *vec, int *size)
{
  register int i;
  TRIANGLE *oldvec = vec;

  vec = (TRIANGLE *) malloc( (size_t) sizeof(TRIANGLE) * 2 * *size);
  if (!vec) xcError("allocation failure in xcRealloc2xBiggerTriangl");
  for(i=0; i<*size; i++) {
    vec[i].p[0] = oldvec[i].p[0];
    vec[i].p[1] = oldvec[i].p[1];
    vec[i].p[2] = oldvec[i].p[2]; /* tole bi moral biti constructor, pa verjetno ni */
  }
  free( (FREE_ARG) oldvec );
  *size *=2;

  return vec;
}


TRIG_INFO *xcRealloc2xBiggerTrigInfo(TRIG_INFO *vec, int *size)
{
  register int i;
  TRIG_INFO *oldvec = vec;

  vec = (TRIG_INFO *) malloc( (size_t) sizeof(TRIG_INFO) * 2 * *size);
  if (!vec) xcError("allocation failure in xcRealloc2xBiggerTrigInfo");
  for(i=0; i<*size; i++) {
    vec[i].i = oldvec[i].i;
    vec[i].j = oldvec[i].j;
    vec[i].k = oldvec[i].k; 

    vec[i].ver_ind[0] = oldvec[i].ver_ind[0];
    vec[i].ver_ind[1] = oldvec[i].ver_ind[1];
    vec[i].ver_ind[2] = oldvec[i].ver_ind[2];

    vec[i].ver_status[0] = oldvec[i].ver_status[0];
    vec[i].ver_status[1] = oldvec[i].ver_status[1];
    vec[i].ver_status[2] = oldvec[i].ver_status[2];
  }
  free( (FREE_ARG) oldvec );
  *size *=2;

  return vec;
}


LINE *xcRealloc2xBiggerLINE(LINE *vec, int *size)
{
  register int i;
  LINE *oldvec = vec;

  vec = (LINE *) malloc( (size_t) sizeof(LINE) * 2 * *size);
  if (!vec) xcError("allocation failure in xcRealloc2xBiggerLINE");
  for(i=0; i<*size; i++) {
    vec[i].p1.x = oldvec[i].p1.x;
    vec[i].p1.y = oldvec[i].p1.y;
    vec[i].p1.z = oldvec[i].p1.z;
    vec[i].p2.x = oldvec[i].p2.x;
    vec[i].p2.y = oldvec[i].p2.y;
    vec[i].p2.z = oldvec[i].p2.z;
  }    

  free( (FREE_ARG) oldvec );
  *size *= 2;

  return vec;
}
