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
 * Source: $XCRYSDEN_TOPDIR/C/sInfo.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <stdio.h>
#include <stdlib.h>
#include <tcl.h>
#include "struct.h"

#define SINFO_MAXARGS 50

void Set_sInfoArray( Tcl_Interp *interp );

void
Set_sInfoArray( Tcl_Interp *interp ) 
{
  char *name1 = "sInfo";
  char *name2[SINFO_MAXARGS];
  char *newValue[SINFO_MAXARGS];
  int i, flags = TCL_GLOBAL_ONLY;

  for(i=0; i<SINFO_MAXARGS; i++)
    newValue[i] = (char *) malloc( sizeof(char) * 256 );

  name2[0] = "dim";
  sprintf(newValue[0], "%d", xcr.dim);
  
  name2[1] = "groupn";
  sprintf(newValue[1], "%d", xcr.groupn);

  name2[2] = "ldimgroup";
  sprintf(newValue[2], "%d", xcr.ldimgroup);

  name2[3] = "lconvvec";
  sprintf(newValue[3], "%d", xcr.lconvvec);

  name2[4] = "lprimvec";
  sprintf(newValue[4], "%d", xcr.lprimvec);

  name2[5] = "lrecprimvec";
  sprintf(newValue[5], "%d", xcr.lrecprimvec);

  name2[6] = "lrecconvvec";
  sprintf(newValue[6], "%d", xcr.lrecconvvec);

  name2[7] = "lprimcoor";
  sprintf(newValue[7], "%d", xcr.lprimcoor);

  name2[8] = "lconvcoor";
  sprintf(newValue[8], "%d", xcr.lconvcoor);

  name2[9] = "lprimwigner";
  sprintf(newValue[9], "%d", xcr.lprimwigner);

  name2[10] = "lconvwigner";
  sprintf(newValue[10], "%d", xcr.lconvwigner);

  name2[11] = "lprimbz";
  sprintf(newValue[11], "%d", xcr.lprimbz);

  name2[12] = "lconvbz";
  sprintf(newValue[12], "%d", xcr.lconvbz);

  name2[13] = "ldatagrid2D";
  sprintf(newValue[13], "%d", xcr.ldatagrid2D);

  name2[14] = "ldatagrid3D";
  sprintf(newValue[14], "%d", xcr.ldatagrid3D);

  name2[15] = "natr";
  sprintf(newValue[15], "%d", xcr.natr);

  name2[16] = "ncell";
  sprintf(newValue[16], "%d", xcr.ncell);

  name2[17] = "nunit";
  sprintf(newValue[17], "%d %d %d", xcr.nunit[0], xcr.nunit[1], xcr.nunit[2]);
 
  name2[18] = "unit";
  sprintf(newValue[18], "%d", xcr.unit);

  name2[19] = "celltype";
  sprintf(newValue[19], "%d", xcr.celltype);

  name2[20] = "shape";
  sprintf(newValue[20], "%d", xcr.shape);

  name2[21] = "lforce";
  sprintf(newValue[21], "%d", xcr.lforce);


  name2[22] = "primvec";
  sprintf(newValue[22], "%15.10f %15.10f %15.10f   %15.10f %15.10f %15.10f   %15.10f %15.10f %15.10f", 
	  vec.prim[0][0], vec.prim[0][1], vec.prim[0][2],
	  vec.prim[1][0], vec.prim[1][1], vec.prim[1][2],
	  vec.prim[2][0], vec.prim[2][1], vec.prim[2][2]);

  name2[23] = "convvec";
  sprintf(newValue[23], "%15.10f %15.10f %15.10f   %15.10f %15.10f %15.10f   %15.10f %15.10f %15.10f", 
	  vec.conv[0][0], vec.conv[0][1], vec.conv[0][2],
	  vec.conv[1][0], vec.conv[1][1], vec.conv[1][2],
	  vec.conv[2][0], vec.conv[2][1], vec.conv[2][2]);

  for (i=0; i<24; i++)
    Tcl_SetVar2(interp, name1, name2[i], newValue[i], flags);
  
  /* free *newValue[] */
  for(i=0; i<SINFO_MAXARGS; i++)
    free((FREE_ARG) newValue[i]);
}
