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
 * Source: $XCRYSDEN_TOPDIR/C/memory.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "memory.h"

void *xcMalloc(size_t size) {
  void *ptr;
  if (size<=0) {
    xcError("can't locate zero memory in xcMalloc");
  }
  ptr = (void *) malloc(size);
  if (!ptr) xcError("out of memory in xcMalloc");
  return ptr;
}


void *xcCalloc(size_t nmemb, size_t size) {
  void *ptr;
  if (size<=0 || nmemb<=0) {
    xcError("can't locate zero memory in xcCalloc");  
  }
  ptr = (void *) calloc( nmemb, size );
  if (!ptr) xcError("out of memory in xcCalloc");
  return ptr;
}


void *xcRealloc(void *ptr, size_t size) {
  ptr = (void *) realloc( ptr, size );
  if (!ptr && size) xcError("out of memory in xcRealloc");
  return ptr;
}


void xcFree(FREE_ARG ptr) {
  if (ptr) free (ptr);
  /*fprintf(stderr,"xcFree();");*/
}

void *xcNullFree(FREE_ARG ptr) {  
  if (ptr) free (ptr);
  return (FREE_ARG)NULL;
}

void 
xcError(char error_text[])
{
  /* XCrySDen standard error handler */  
  fprintf(stderr,"XCrySDen run-time error:\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"Exit !!!\n");
  exit(1);
}

void
xcFree_Vectorf(float *v)
{
  free( (FREE_ARG) v );
}

void
xcFree_Matrixf(float **m)
{
  free( (FREE_ARG) m[0] );
  free( (FREE_ARG) m );
}

void
xcFree_Tensor3f(float ***t)
{
  free( (FREE_ARG) t[0][0] );
  free( (FREE_ARG) t[0] );
  free( (FREE_ARG) t );
}


void
xcFree_Tensor3i(int ***t)
{
  free( (FREE_ARG) t[0][0] );
  free( (FREE_ARG) t[0] );
  free( (FREE_ARG) t );
}


void
xcFree_ReallocatedTensor3f_00to11( float ***rt )
{
  free( (FREE_ARG) (rt[1] + 1) );
  free( (FREE_ARG) (rt + 1) );
}


float *xcMallocVectorf(long n)
{
  float *v;
  v=(float *) malloc( (size_t) n*sizeof(float) );
  if (!v) xcError("allocation failure in xcMallocVectorf");
  return v;
}


float **xcMallocMatrixf(long nr, long nc)
{
  long i;
  float **m;
  /* allocate pointers to pointers to rows */
  m=(float **) malloc( (size_t) nr*sizeof(float *) );
  if (!m) xcError("allocation failure 1 in xcMallocMatrixf");

  /* allocate pointers to rows and set pointers to them */
  m[0]=(float *) malloc( (size_t) nr*nc*sizeof(float) );
  if (!m[0]) xcError("allocation failure 2 in xcMallocMatrixf");

  for(i=1; i<nr; i++)
    m[i] = m[i-1] + nc;
  
  return m;
}


float ***xcMallocTensor3f(long nr, long nc, long nd)
{  
  /* allocate a float 3tensor with range t[0..nr-1][0..nc-1][0..nd-1] */
  long i, j;
  float ***t;
  /* allocate pointers to pointers to rows */
  t=(float ***) malloc( (size_t) nr*sizeof(float**) );
  if (!t) xcError("allocation failure 1 in xcMallocTensor3f");

  /* allocate pointers to rows and set pointers to them */
  t[0]=(float **) malloc( (size_t) nr*nc*sizeof(float*) );
  if (!t[0]) xcError("allocation failure 2 in xcMallocTensor3f");

  /* allocate rows and set pointers to them */
  t[0][0]=(float *) malloc( (size_t) nr*nc*nd*sizeof(float) );
  if (!t[0][0]) xcError("allocation failure 3 in xcMallocTensor3f");
  
  for(i=0; i<nr; i++) {
    t[i] = t[0] + i*nc;
    for(j=0; j<nc; j++)
      t[i][j] = t[0][0] + nd*(i*nc + j);
  }
  return t;
}

int ***xcMallocTensor3i(long nr, long nc, long nd)
{  
  /* allocate a int 3tensor with range t[0..nr-1][0..nc-1][0..nd-1] */
  long i, j;
  int ***t;
  /* allocate pointers to pointers to rows */
  t=(int ***) malloc( (size_t) nr*sizeof(int**) );
  if (!t) xcError("allocation failure 1 in xcMallocTensor3i");

  /* allocate pointers to rows and set pointers to them */
  t[0]=(int **) malloc( (size_t) nr*nc*sizeof(int*) );
  if (!t[0]) xcError("allocation failure 2 in xcMallocTensor3i");

  /* allocate rows and set pointers to them */
  t[0][0]=(int *) malloc( (size_t) nr*nc*nd*sizeof(int) );
  if (!t[0][0]) xcError("allocation failure 3 in xcMallocTensor3i");
  
  for(i=0; i<nr; i++) {
    t[i] = t[0] + i*nc;
    for(j=0; j<nc; j++)
      t[i][j] = t[0][0] + nd*(i*nc + j);
  }
  return t;
}

short int ***xcMallocTensor3si(long nr, long nc, long nd)
{  
  /* allocate a int 3tensor with range t[0..nr-1][0..nc-1][0..nd-1] */
  long i, j;
  short int ***t;
  /* allocate pointers to pointers to rows */
  t=(short int ***) malloc( (size_t) nr*sizeof(short int**) );
  if (!t) xcError("allocation failure 1 in xcMallocTensor3i");

  /* allocate pointers to rows and set pointers to them */
  t[0]=(short int **) malloc( (size_t) nr*nc*sizeof(short int*) );
  if (!t[0]) xcError("allocation failure 2 in xcMallocTensor3i");

  /* allocate rows and set pointers to them */
  t[0][0]=(short int *) malloc( (size_t) nr*nc*nd*sizeof(short int) );
  if (!t[0][0]) xcError("allocation failure 3 in xcMallocTensor3i");
  
  for(i=0; i<nr; i++) {
    t[i] = t[0] + i*nc;
    for(j=0; j<nc; j++)
      t[i][j] = t[0][0] + nd*(i*nc + j);
  }
  return t;
}

float ***xcReallocTensor3f_00to11(float ***t, long nr, long nc, long nd)
{
  /* convert from t[0..nr-1][0..nc-1][0..nd-1] to m[1..nr][1..nc][1..nd] */
  long i, j;
  float ***m;
  
  /* allocate pointers to pointers to rows */
  m=(float ***) malloc((size_t) (nr * sizeof(float **)));
  if (!m) xcError("allocation failure #1 in ***xcReallocTensor3f_00to11");
  
  m -= 1; /* now m[0] is undefined & m[1] == previous_m[0]; */
  m[1]=(float **) malloc((size_t) (nr * nc * sizeof(float *)));
  if (!m[1]) xcError("allocation failure #2 in ***xcReallocTensor3f_00to11");
  m[1] -= 1;  /* &(m[i][0]) -= 1 */
  for(i=1; i<=nr; i++) {
    m[i] = m[1] + (i-1)*nc;
    for(j=1; j<=nc; j++) {
      /* m[i][j] -= 1; */
      m[i][j] = t[i-1][j-1] - 1;
    }
  }
  return m;
}


void *cryMem_Malloc_W( size_t size, char *where ) {
  void *ptr = (void *) malloc( size );
  if (!ptr) {
    fprintf( stderr, 
	     "\n***\nMemory Error: can't allocate %ld bytes in %s\n", 
	     size, where );
    abort();
  }
  return ptr;
}
