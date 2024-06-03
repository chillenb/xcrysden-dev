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
 * Source: $XCRYSDEN_TOPDIR/C/remakestr.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <togl.h>
#include <math.h>
#include <stdio.h>
#include "struct.h"
#include "3D.h"
#include "memory.h"
#include "xcfunc.h"

extern Options3D is;

/* --- function protoypes ---*/
void ReMakeStr(void);
static void ReMakeBonds(void);
void ReDisplay(struct Togl *togl, GLuint remake3D);

/* --- extern function protoypes ---*/

/* --- xcDisplayFunc.c ---*/
/*extern void (*xcDisplay)(void);*/
/* extern void xcMakeBallLists(void); */
extern void xcMakeSticklLists(void);
extern void xcMakeSpaceFillLists(void);
/*  extern void xcMakeFrame3DLists(void); */
/*  extern void xcMakeBallLabel3DList(void); */
/*  extern void xcMakeSpaceLabel3DList(void); */

/* --- xcDisplayFunc2.c ---*/
extern void UpdateProjection(void);
extern void WhichDispFunc(void);

/* --- xcviewport.c --- */
extern void xcViewPort();

/* ---- readstrf.c ---- */
extern void MallocCoor(void);

/* ---- 3D.c ---- */
extern void MakeCylinderCoor(void);

/* ---- xcballstick.c ---- */
extern void MakeArcPoints(void);

/* ---- readstrf.c ---- */
extern void MakeBonds(void);


/* if "L_COV_SCALE", than chemical connectivity criteria is changed
 * this means that we must reassing bonds; 
 */
void
ReMakeStr(void)
{
  /* now radius has been changed so I must just remake bonds, that
   * is: first determine nbonds, than assign xbond, ..., assign
   * nobjects, MallocCoor & MakeCylinderCoor
   */
  
  /* specially care must be taken for "tmp_nobjects", because "nbonds" 
     parameter will be changed after ReMakeBonds */
  printf("ReMakeStr:: tmp_nobjects_BEFORE = %d\n",tmp_nobjects);
  printf("ReMakeStr:: nbonds_BEFORE = %d\n",nbonds);
  printf("ReMakeStr:: natoms_BEFORE = %d\n",natoms);
  printf("ReMakeStr:: nframes_BEFORE = %d\n",nframes);
  printf("ReMakeStr:: nobjects_BEFORE = %d\n",nobjects);
  tmp_nobjects = tmp_nobjects - nbonds;

  ReMakeBonds();
  
  tmp_nobjects = tmp_nobjects + nbonds;

  nobjects = natoms + nbonds + nframes;
  printf("ReMakeStr:: tmp_nobjects_AFTER = %d\n",tmp_nobjects);
  printf("ReMakeStr:: nbonds_AFTER = %d\n",nbonds);
  printf("ReMakeStr:: natoms_AFTER = %d\n",natoms);
  printf("ReMakeStr:: nframes_AFTER = %d\n",nframes);
  printf("ReMakeStr:: nobjects_AFTER = %d\n",nobjects);
  
  /* malloc "coor" & "coor3D" and stuff that goes nearby */
  MallocCoor();
  /* rewrite coor structure */  
  WhichDispFunc();
  
  /* --- now calculate coordinates for 3Dbonds & 3Dframes --- */
  MakeCylinderCoor();
}

 
void
ReDisplay(struct Togl *togl, GLuint remake3D)
{
  /* first make new projection  (glOrtho) */
  UpdateProjection();
  
  /*
  if (is.stickmode) xcMakeStick3DLists();
  if (is.ballmode)  xcMakeBall3DLists();
  if (is.spacefillmode) xcMakeSpaceFill3DLists();
  */

  /*    if (remake3D & REMAKE3D_FRAME ) xcMakeFrame3DLists(); */
  /*    if (remake3D & REMAKE3D_BALL_LABEL) xcMakeBallLabel3DList(); */
  /*    if (remake3D & REMAKE3D_SPACE_LABEL) xcMakeSpaceLabel3DList(); */
  
  /* this is safe enough, because if some ball radius has changed, than
     we check which balllist has to be "new-made" */

  if (dimType == XC_2D) {
    /* xcMakeBallLists();*/
    /* remake MakeArcPoints for ballstick #2 */
    MakeArcPoints();
  }
  /* if we changed radius we had to make new glOrtho */
  /* ############## TUKI ############### */
  /* now update a display */
  xcViewPort();
  Togl_PostRedisplay(togl);
}


static void 
ReMakeBonds(void)
{
  /* number of bonds will changed, so we must free, and than malloc 
     (in MakeBonds) the following */
  xcFree((void*) natbond);
  xcFree((void*) sqnbond); 
  xcFree((void*) bondend);
  xcFree((void*) xbond); 
  xcFree((void*) ybond); 
  xcFree((void*) zbond);  
  xcFree((void*) xbond2); 
  xcFree((void*) ybond2); 
  xcFree((void*) zbond2); 

  /* also free "coor" structure */
  xcFree((void*) coor);
  xcFree((void*) zorient);
  xcFree((void*) iwksp);

  MakeBonds();
  make_H_Bonds();
}
