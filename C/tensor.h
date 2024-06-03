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
 * Source: $XCRYSDEN_TOPDIR/C/tensor.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#define MALLOC_TENSOR2(type, p, n1, n2) \
do { \
register long _i, _n1 = (n1), _n2 = (n2); \
p = (type**) malloc((size_t) _n1 * sizeof(type*)); \
if (!p) xcError("memory exhausted"); \
p[0] = (type*) malloc((size_t) _n1*_n2 * sizeof(type)); \
if (!p[0]) xcError("memory exhausted"); \
 for (_i=1; _i<_n1; _i++) { \
   p[_i] = p[_i-1] + n2; \
 } \
} while(0)

#define FREE_TENSOR2(tensor) \
do { \
if (tensor[0]) free( (FREE_ARG) tensor[0] ); \
if (tensor)    free( (FREE_ARG) tensor ); \
} while(0)


#define MALLOC_TENSOR3(type, p, n1, n2, n3) \
do { \
register int _i, _j, _n1 = (n1), _n2 = (n2), _n3 = (n3); \
p = (type***) cryMem_Malloc_W((size_t) _n1 * sizeof(type**), \
            "MALLOC_TENSOR3: allocation failure 1"); \
p[0] = (type**) cryMem_Malloc_W((size_t) _n1*_n2 * sizeof(type*), \
            "MALLOC_TENSOR3: allocation failure 2"); \
p[0][0] = (type*) cryMem_Malloc_W((size_t) _n1*_n2*_n3 * sizeof(type), \
	    "MALLOC_TENSOR3: allocation failure 3"); \
 /*printf("MALLOC_TENSOR3: allocation completed !!!\n");*/ \
 for (_i=0; _i<_n1; _i++) { \
   p[_i] = p[0] + _i*_n2; \
  /* printf("MALLOC_TENSOR3: array assignment -> major index = %d/%d\n", _i,_n1);*/ \
   for (_j=0; _j<_n2; _j++) { \
     /*printf("MALLOC_TENSOR3: array assignment -> major - minor index = [%d/%d][%d/%d]\n", _i,_n1,_j,_n2);*/ \
     p[_i][_j] = p[0][0] + _n3*(_i*_n2 + _j); \
   } \
 } \
} while(0)
