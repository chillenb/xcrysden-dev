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
 * Source: $XCRYSDEN_TOPDIR/C/vector.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#define XC_CPP_VECTOR

#define VECTOR_ARROWSIZE 0.2
#define VECTOR_THICKF    1.3 /* basic vector thickness factor */
#define VECTOR_ARRTHICKF 1.5

#define VECTOR_PRIMITIVE    1
#define VECTOR_CONVENTIONAL 2
#define VECTOR_SUPERCELL    4

typedef struct {
  GLdouble coor[3][3];   /* vec[#].coor[x/y/z][begin/arrow/end] */
  GLdouble vecthick;
  GLdouble vecx; 
  GLdouble vecy; 
  GLdouble vecz;
  GLdouble vecfi; 
  GLdouble vecl;
  GLdouble arrthick;
  GLdouble arrx; 
  GLdouble arry; 
  GLdouble arrz;
  GLdouble arrfi; 
  GLdouble arrl;
} RenderVectors;
RenderVectors rnd_pvec[3];  /* primitive vectors */ 
RenderVectors rnd_cvec[3];  /* conventional vectors */
RenderVectors rnd_scvec[3]; /* supercell vectors */   
			    
typedef struct {	    
  GLfloat vertex[6][4][3];
  GLfloat normal[6][4][3];
} CellCage;
CellCage prim_cage;  /* primitive cell-cage */    
CellCage conv_cage;  /* conventional cell-cage */
CellCage sc_cage;    /* supercell cell-cage */   

