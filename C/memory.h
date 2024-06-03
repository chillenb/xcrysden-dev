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
 * Source: $XCRYSDEN_TOPDIR/C/memory.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#ifndef FREE_ARG
#   define FREE_ARG char*
#endif

/* memory.c */
extern void *xcMalloc(size_t size);
extern void *xcCalloc(size_t nmemb, size_t size);
extern void *xcRealloc(void *ptr, size_t size);
extern void xcFree(FREE_ARG ptr);
extern void *xcNullFree(FREE_ARG ptr);
extern void xcError(char error_text[]);
extern void xcFree_Vectorf(float *v);
extern void xcFree_Matrixf(float **m);
extern void xcFree_Tensor3f(float ***t);
extern void xcFree_Tensor3i(int ***t);
extern void xcFree_ReallocatedTensor3f_00to11(float ***rt);
extern float *xcMallocVectorf(long n);
extern float **xcMallocMatrixf(long nr, long nc);
extern float ***xcMallocTensor3f(long nr, long nc, long nd);
extern int ***xcMallocTensor3i(long nr, long nc, long nd);
extern short int ***xcMallocTensor3si(long nr, long nc, long nd);
extern float ***xcReallocTensor3f_00to11(float ***t, long nr, long nc, long nd);
extern void *cryMem_Malloc_W(size_t size, char *where);

#ifdef XC_CPP_ISOSURF
   extern void xcFree_PLANEVERTEX(PLANEVERTEX **m);
   extern void xcFree_GRIDVERTEX(GRIDVERTEX ***t);
   extern void xcFree_LINE(LINE *t);
   extern PLANEVERTEX **xcMallocPLANEVERTEX(long nr, long nc);
   extern GRIDVERTEX ***xcMallocGRIDVERTEX(long nr, long nc, long nd);
   extern LINE *xcMallocLINE(long nr);
   extern VERTEX *xcRealloc2xBiggerVertex(VERTEX *vec, int *size);
   extern TRIANGLE *xcRealloc2xBiggerTriangl(TRIANGLE *vec, int *size);
   extern LINE *xcRealloc2xBiggerLINE(LINE *vec, int *size);
#endif
