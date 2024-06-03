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
 * Source: $XCRYSDEN_TOPDIR/C/molsurf.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/
#ifndef MOLSURF_H
#define MOLSURF_H

#define MOLS_GAUSSIAN     0
#define MOLS_EXP          1
#define MOLS_UNIGAUSS     2
#define MOLS_UNIEXP       3
#define MOLS_DISTFUNC     4
#define MOLS_FERMISURFACE 10
#define MOLS_COV          0
#define MOLS_VdW          1
#define MOLS_ATOMIC       0
#define MOLS_MONOCHROME    1
#define MOLS_DISTANCE     2
#define MOLS_MOLSURF      0
#define MOLS_GAP          1

typedef struct {
  short int gridindex, gridsubindex, bandindex, celltype, cropbz;
  float wirecellcolor[4], solidcellcolor[4];
} FS_DATA;


#ifndef XC_CPP_GLPARAM
#   include "xcGLparam.h"
#endif
#ifndef ISOSURF_H
#   include "isosurf.h"
#endif

struct MOL_SURF_struct {
  struct MOL_SURF_struct *ptr;
  int index;
  /* ----------------------- */  
  int isosurf_index; /* index for isusurface */
  /* ----------------------- */  
  short int  type, colorscheme, surftype;
  ISO_ATTRIB dispt;      /* display types */
  LIGHTMODEL lightm;     /* lightmodels */
  IsoExpand  isoexpand;   /* periodic attribute for molsurf */
  int     smooth_nstep;  /* # of smothing steps */
  float   smooth_weight; /* smoothing weight */
  int     frontface;
  int     revertnormals;
  double  *rad;
  float   isolevel;
  float   d;           /* resolution criteria */
  int     i,j,k;       /* grid size */
  float   lowcoor[3];  /* point with lowest coordinates */  
  float   size[3];     /* size of gridcube */
  float   cutoff;      /* cutoff radius (to determine how big the 
			  grid for triangulation must be) */
  float   r;           /* current cnl's radius */
  int     *atmindex;   /* here goes indices of atoms withing cutoff region   */
  int     natm;        /* number of atoms that are within cutoff region      */
  float   c[MAXNAT+1]; /* exponent for MolSurfFunc */
  float   monocolor[4], back_monocolor[4];
  float   vec[3][3];
  float   ***molgrid;
  float   (*MolSurfFunc)(struct MOL_SURF_struct *mols, int nat, float dis); /* pointer to MolSurfFunction */

  /*---- Fermi-surface ----*/
  FS_DATA fs;

  /*---- interpolation ----*/
  float   ***interp_molgrid;
  int     interp_n[3];
  int     tessellation_algorithm;
  int     normals_algorithm;
};
typedef struct MOL_SURF_struct MOL_SURF;

struct READSURF_OPT {
  char  *option;
  short int flag;
};

#endif
