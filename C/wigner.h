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
 * Source: $XCRYSDEN_TOPDIR/C/wigner.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#define WIGNER_PRIM  0
#define WIGNER_CONV  1
#define WIGNER_CLEAR 2

#define WIGNER_NODEMODE_EVERY  0     /* every node is selected (DEFAULT) */
#define WIGNER_NODEMODE_SELECT 1     /* just selected nodes */

typedef struct {
  FILE  *fp;
  int   transparency;
  float color[4];
  int   is_initialized[2];
  int   nodemode[2];
  int   nnodes[2];
  float **nodes[2]; /* nodes[prim/conv][#][fx/fy/fz] */
} WIGNER_ATTRIB;
WIGNER_ATTRIB ws_attrib;



