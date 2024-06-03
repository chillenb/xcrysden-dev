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
 * Source: $XCRYSDEN_TOPDIR/C/struct.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 * ------                                                                    *
 * Modified by Eric Verfaillie ericverfaillie@yahoo.fr EV                    *
 * may 2004                                                                  *
 * modifcations are near EV comments                                         *
  *****************************************************************************

*/

/* 
   handle this near & far bug (needed for CYGWIN)
*/
#ifndef STRUCT_H_
#define STRUCT_H_

#ifdef NEAR_BUG
#ifdef near
#undef near
#endif
#ifdef far
#undef far
#endif
#endif /* NEAR_BUG */

/*****************************************************************************/
/* XC_CPP_* are pre-processor flags */
#define XC_CPP_STRUCT

#ifndef FREE_ARG
#  define FREE_ARG char*
#endif

#define XC_CPP_GLH
#include <GL/gl.h>

#ifndef MAX
#define MAX(a,b) ( (a)>(b)?(a):(b) )
#endif
#ifndef MIN
#define MIN(x,y) ( (x)<(y) ? (x) : (y) )
#endif
#ifndef ABS
#define ABS(x)   ( (x)>0?(x):-(x)  )
#endif

#define DOT_V3(a,b)  ( (a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2] )
#define NORM_V3(v)   ( sqrt( (v)[0]*(v)[0] + (v)[1]*(v)[1] + (v)[2]*(v)[2] ) )

#define LOAD_NULL_V(n, vec) \
do { \
  register int __iv; \
  for (__iv=0; __iv<n; __iv++) (vec)[__iv] = 0; \
} while(0)

#define COPY_V(n, dst, src) \
do { \
  register int __iv; \
  for (__iv=0; __iv<n; __iv++) (dst)[__iv] = (src)[__iv]; \
} while(0)

#define COPY_AND_SCALE_V(n, dst, scale, src) \
do { \
  register int __iv; \
  for (__iv=0; __iv<n; __iv++) (dst)[__iv] = scale * (src)[__iv]; \
} while(0)  

#define CLAMP_V(n, value, v) \
do { \
  register int __iv; \
  for (__iv=0; __iv<n; __iv++) \
    v[__iv] = v[__iv] > value ? value : v[__iv]; \
} while(0);

/*
  I think this is not needed ...

  #define COPY_V_(n, dst, src) \
  do { \
  int __iv; \
  for (__iv=0; __iv<n; __iv++) dst[__iv] = src[__iv]; \
  } while(0)

*/
#define COPY_M33(dst, src) \
do { \
  (dst)[0][0] = (src)[0][0]; (dst)[0][1] = (src)[0][1]; (dst)[0][2] = (src)[0][2]; \
  (dst)[1][0] = (src)[1][0]; (dst)[1][1] = (src)[1][1]; (dst)[1][2] = (src)[1][2]; \
  (dst)[2][0] = (src)[2][0]; (dst)[2][1] = (src)[2][1]; (dst)[2][2] = (src)[2][2]; \
} while(0)

#define INVERT_M33(type, inverse, direct) \
do { \
  type _d_e_t_, _m33_[3][3];                                           \
                                                                       \
  _m33_[0][0] = (direct)[1][1]*(direct)[2][2] - (direct)[1][2]*(direct)[2][1]; \
  _m33_[0][1] = (direct)[1][2]*(direct)[2][0] - (direct)[1][0]*(direct)[2][2]; \
  _m33_[0][2] = (direct)[1][0]*(direct)[2][1] - (direct)[2][0]*(direct)[1][1]; \
  _m33_[1][1] = (direct)[0][0]*(direct)[2][2] - (direct)[0][2]*(direct)[2][0]; \
  _m33_[1][2] = (direct)[2][0]*(direct)[0][1] - (direct)[0][0]*(direct)[2][1]; \
  _m33_[2][2] = (direct)[0][0]*(direct)[1][1] - (direct)[0][1]*(direct)[1][0]; \
  _m33_[1][0] = (direct)[2][1]*(direct)[0][2] - (direct)[0][1]*(direct)[2][2]; \
  _m33_[2][0] = (direct)[0][1]*(direct)[1][2] - (direct)[1][1]*(direct)[0][2]; \
  _m33_[2][1] = (direct)[0][2]*(direct)[1][0] - (direct)[0][0]*(direct)[1][2]; \
                                                                       \
  _d_e_t_ =                                                            \
    (direct)[0][0]*_m33_[0][0] +                                       \
    (direct)[0][1]*_m33_[0][1] +                                       \
    (direct)[0][2]*_m33_[0][2];                                        \
                                                                       \
  /* make a check for zero devision */                                 \
  /* ... here ... */                                                   \
                                                                       \
  (inverse)[0][0] = _m33_[0][0] / _d_e_t_;                             \
  (inverse)[1][1] = _m33_[1][1] / _d_e_t_;                             \
  (inverse)[2][2] = _m33_[2][2] / _d_e_t_;                             \
  (inverse)[0][1] = _m33_[1][0] / _d_e_t_;                             \
  (inverse)[0][2] = _m33_[2][0] / _d_e_t_;                             \
  (inverse)[1][2] = _m33_[2][1] / _d_e_t_;                             \
  (inverse)[1][0] = _m33_[0][1] / _d_e_t_;                             \
  (inverse)[2][0] = _m33_[0][2] / _d_e_t_;                             \
  (inverse)[2][1] = _m33_[1][2] / _d_e_t_;                             \
} while(0)

#define XC_OK        1
#define XC_ERROR     0
#define XC_TRUE      1
#define XC_FALSE     0
#define XC_LEFT      1
#define XC_RIGHT     0
#define XC_YES       1
#define XC_NO        0
#define XC_ALLOWED   1 
#define XC_FORBIDDEN 0
#define MAXNAT       100   /* maximum atomic number */
#define MINDIS       0.0001  /* minimum disitance between atoms (in A) allowed */
#define MINTOL       1.0e-15
#define MINTOL_F     1.0e-5

#ifndef PI
#  define PI           3.14159265358979323844
#endif
#define PI12         1.57079632679489661922
#define RAD2DEG      57.295779513084
#define OUTLINE      1.0   /* this is for linewidth when rendering balls */
#define CIRCLE       6     /* determine how precise/smooth a 2Dballs are drawn
		            * CIRCLE = 1 -> maximum smoothness
		            */
#define FRAME_PAT    0xAAAA /* this is stipple pattern for frames */

#define MAXSEL       25     /* maximum number of selected atoms */
#define N_SELECT_OBJ 2*MAXSEL     /* maximum number of selected objects (atoms+bonds) */

/*============================================================================
 * define flags for vorius objects 
 * this is used for glLists 
 */ 
#define N_LIST_PRIMITIVES 11

#define POINT           0
#define BALL            1
#define OUTLINEBALL     2

#define SPHEREWIRE	3
#define CYLINDERWIRE	4
#define CONEWIRE        5
#define BONDWIRE        6

#define SPHERESOLID	7
#define CYLINDERSOLID	8
#define CONESOLID       9
#define BONDSOLID      10

#define PIPEBALL       11

/*****************************************************************************
 * WARNING::                                                                 *
 *          here below are few "global" variables and this are the only      *
 *          variables that are not structures; all the others are structures *
 *****************************************************************************/

/*===============================================================*/
/* --- ATOMS & BONDS --- ATOMS & BONDS --- BEGIN --- BEGIN ---   */
/*===============================================================*/
#define MAX_BONDS 30

int nbonds, natoms, nframes; /* number of bonds, number of atoms,
				       number of "crystal frame-pieces" */

int tmp_nobjects;  /* number of objects to display: 
		  nobjects = if(nbonds) + if(natoms) + if(nframes) */
int nobjects; /* total number of objects: 
		 nobjects = nbonds + natoms + nframes */

double *xat, *yat, *zat;           /* XYZ coordinates of atoms */
double (*fv)[3];                     /* force vector */
double *xbond, *ybond, *zbond;     /* XYZ of first end of a bond 
				    * First End of a bond is always at atom
				    */
double *xbond2, *ybond2, *zbond2;  /* second end of a bond */

#define BOND_ATOM_TO_MIDBOND 1
#define BOND_MIDBOND_TO_ATOM 2
#define BOND_ATOM_TO_ATOM    4


#define VORONOI_DEFAULT     0
#define VORONOI_WIGNERSEITZ 1
#define VORONOI_BZ          2


#define CENTERED_LABEL -999.999



int *bondend; /* how the bond is oriented, ie. from atom to midbond or 
		 from midbond to atom !!! */

int *frametype;   /* flag for frame type (I/O connectivity) */ 
double *xframe, *yframe, *zframe;   /* XYZ of first end of a frame 
				     * First End of a frame is always at atom
				     */
double *xframe2, *yframe2, *zframe2; /* second end of a frame */


double *zbmid; /* center of a bond in Z-direction */

int *nat;  /* atomic number of an atom */
int *sqn;  /* sequential number of an atom */
int *natbond; /* to determine to what kind of an atom the bond bellong;
	       * all the rest of the attrib. will be derived from
	       * this parameter 
	       */
int *sqnbond; /* will replace natbond;
	       * to determine to what kind of an atom the bond bellong;
	       * all the rest of the attrib. will be derived from
	       * this parameter 
	       */
int *iwksp;   /* this is for Z-orientation (F77_INDEXX1) */
double mx, my, mz; /* geometrical centre of structure is always moved to
		    * (0,0,0); (mx,my,mz) is translation vector for that  
		    */

typedef struct {
  double mx, my, mz;
} MX_MY_MZ;
MX_MY_MZ mx_my_mz;
 
char displayMode2D[4];

/*GLuint *primitive_Lists[N_LIST_PRIMITIVES]; */ /* holds the glLists ID for 
					            various primitives */

/*=========================================================================*/
/*     above were 1st part of global variables that are not structures     */
/*=========================================================================*/

/* definitions for iso_objects:: */
#define ISOOBJ_BASE    0   /* base isosurface or colorplane/isolines */
#define ISOOBJ_PLANE1  1
#define ISOOBJ_PLANE2  2
#define ISOOBJ_PLANE3  3
#define MAX_ISOOBJECTS 4   /* for example: isosurface & three colorplanes */

/* -------PARAMETERS NEEDED FOR STRUCTURE DISPLAY -----------------*/
typedef struct { 
  double x;                 /* xcViewPort's x */
  double y;                 /* xcViewPort's y */
  double sizex;          /* xcViewPort's size when X-croped */
  double sizey;          /* xcViewPort's size when Y-croped */
  double size;              /* xcViewPort's size when not croped */
  int stropened;         /* when structure is opened this parameter is 
			    set to 1 
			  */
  int canvassize;       /* size of a canvas in screen units
		           if height < width ==> canvassize = height
		           if Ycnv > Xcnv ==> canvassize = width
		         */
  int width;            /* canvas width */
  int height;           /* canvas height */
  double VPfactor;      /* VPf.VPfactor = VPf.size / VPf.structsize; */
  float scrAnX;         /* "screen units to angstroms" factor for X-dir. */
  float scrAnY;         /* "screen units  to angstroms" factor for Y-dir. */
  int WFlinewidth;      /* linewidth for WireFrame */
  int WF3Dlinewidth;    /* linewidth for WireFrame */
  int PLlinewidth;      /* linewidth for PointLine */
  int OUTlinewidth;
  int PLradius;         /* Point radius for PointLine */

  int framewidth;       /* linewidth for crystal-frames */

  double rcovf;         /* rcov[] = rcovf * rcovdef[]; this is scale factor 
			   for chemical connectivity criteria */
  /* double scrf; */         /* scrf == how biger is canvas than structure */
  /* THIS IS FOR BALLSTICK MODE */
  double ballf;         /* scale for balls; rball = ballf * rcov */
  double rodf;          /* rod/bond scale factor; 
		           rrod = ballf * rcov[H] * rodf */
  double atradf;        /* scaled factor for balls and spacefills; 
		           it's howmany times are spacefills/balls greater
		           than rcov or vdw
		         */
  double tessFactor;    /* tessellation factor for xcPrimitives */
  double framef;        /* analogous to ballf, rodf, just for frames */
  GLboolean xyzOn;      /* true when coordinate sistem is displayed */
  GLboolean labelsOn;   /* true when atomic labels are displayed */
  GLboolean framesOn;   /* TRUE when "crystal" frames are displayed */
  GLboolean selection;  /* true when we are in "selection" */
  GLboolean atomadd;    /* true when we use "cell-adding selection" for 
                           cell-adding type of ATOMINSE crystal command */
  GLboolean projmade;
  GLboolean isosurf[MAX_ISOOBJECTS]; /* 
                              when isosurface is displayed, this is true */
  GLboolean colorplane[MAX_ISOOBJECTS]; /* 
                              when "colorplane is displayed, this is true */
  GLboolean isoline[MAX_ISOOBJECTS]; /* 
                              when "isoline-plane" is displayed this is true */
  GLboolean wignerseitz; /* render wigner-seitz cell */
  GLboolean supercell;   /* in supercell mode -> render vectors and 
			    cell-cage */
  GLboolean isospacesel2D; /* we are in IsoSpaceSel mode -> for isoplane case a
			      parallelogram is drawn; for isosurface case a
			      parallelepiped is drawn */
  GLboolean isospacesel3D;
  GLboolean force;     /* render the forces */
  GLboolean Hbond;     /* render H-bonds */
  GLboolean unibond;
  GLboolean fog;
  GLboolean antialias;
  GLboolean perspective;  
  double    perspective_fovy;
  double    perspective_far;
  double    perspective_size;

  /* NEW; version 0.4.x on */
  int       *surface;      /* new implementation of isosurfaces allow an 
			      arbitrary isosurfaces to be displayed;
			      hence we need a logical pointer */
  int       *surfaceInd;   /* mapping from mols->isosurf_index */
  void      **surfacePtr;  /* pointer to MOL_SURF structure */
  int       nsurface;			        
  int       dispLattice;   /* do we want to display lattice (generated 
			      separately, not read from FRAMES XSF section */
  int       dispLatType;   /* display type for dispLattice */
} StructPar;

StructPar VPf; /* all parameters for structure(molecule, crystal..) display 
		* + state parameters are here */

/*******************************/
/* this is for VPf.dispLatType */
#define CELL_WIRE           0
#define CELL_ROD            1
#define CELL_SOLID          2
#define CELL_SOLID_AND_WIRE 3
#define CELL_SOLID_AND_ROD  4
/*******************************/

typedef struct {
  double structsize;     /* each atomic structure has it's own size; 
		          * STRUCTSIZE is a parameter for this 
			  */
  double isosize;        /* size of isoXXX; xc_isoexpand is taken into account
                          */
  double wignerseitzsize; /* size of wigner-seitz cells */
  double isospaceselsize; /* size for iso_space_selelection 
			     parallelogram/parallelopiped */
  /* used for cryXXX stuff */
  float o_shift[3];  /* origin shift; geometrical center of picture should be 
			in the middle of MODELVIEW */
} ModelPar;
ModelPar MVf;
/*
 * this definitions here are used for ResetVar(int var)
 */
#define R_ATRAD                0
#define R_RCOV                 10000
#define R_RBALL                1
#define R_RROD                 2
#define R_ATCOL                3
#define R_WFLINEWIDTH          4
#define R_WF3DLINEWIDTH        10004
#define R_OUTLINEWIDTH         10005
#define R_PLLINEWIDTH          5
#define R_PLRADIUS             6
#define R_SCRF                 7
#define R_ALL                  8
#define R_FRAMECOL             9
#define R_FRAMELINEWIDTH       10
#define R_FRAMERODF            11
#define R_BACKGROUND           23
#define R_TESSELLATION         24
#define R_UNIBOND              25
#define R_UNIBONDCOLOR         26
#define R_PERSPECTIVE          27
#define R_PERSPECTIVEBACK      28
#define R_PERSPECTIVEFOVY      29
#define R_PERSPECTIVEFRONT     10029
#define R_FOG                  30
#define R_ANTIALIAS            31
#define R_AMBIENT_BY_DIFFUSE   32
#define R_CURRENTFILEFORMAT    33
#define R_HBOND                34
#define R_FORCE_RODTHICKF      35
#define R_FORCE_ARRTHICKF      36
#define R_FORCE_ARRLENF        37
#define R_FORCE_COLOR          38
#define R_XYZ_AXIS_COLOR       39
#define R_XYZ_XYPLANE_COLOR    40

/*
 * this definitions here are used for LoadNewValue(int var, double ....)
 */
#define L_SPACE_COV            0 
#define L_SPACE_VDW            1 
#define L_BALL_COV             2
#define L_BALL_VDW             3
#define L_RCOV_ONE             10004
#define L_ATRAD_ONE            4 
#define L_ATRAD_SCALE          5 
#define L_BALLF                6 
#define L_RODF                 7 
#define L_ATCOL_ONE            8 
#define L_WFLINEWIDTH          9 
#define L_WF3DLINEWIDTH        10009
#define L_OUTLINEWIDTH         10010
#define L_PLLINEWIDTH          10 
#define L_PLRADIUS             11
#define L_SCRF                 12
#define L_COV_SCALE            13
#define L_XYZ_ON               14  /* for coordinate sistem */ 
#define L_LABELS_ON            15  /* for atomic labels */
#define L_FRAME_ON             16
#define L_FRAMECOL             17
#define L_FRAMELINEWIDTH       18
#define L_FRAMERODF            19
#define L_LINEFRAME            20
#define L_TR_XTRANSL           21
#define L_TR_YTRANSL           22
#define L_BACKGROUND           23
#define L_TESSELLATION         24
#define L_UNIBOND              25
#define L_UNIBONDCOLOR         26
#define L_PERSPECTIVE          27
#define L_PERSPECTIVEBACK      28
#define L_PERSPECTIVEFOVY      29
#define L_PERSPECTIVEFRONT     10029
#define L_FOG                  30
#define L_ANTIALIAS            31
#define L_AMBIENT_BY_DIFFUSE   32
#define L_CURRENTFILEFORMAT    33
#define L_HBOND                34
#define L_FORCE_RODTHICKF      35
#define L_FORCE_ARRTHICKF      36
#define L_FORCE_ARRLENF        37
#define L_FORCE_COLOR          38
#define L_XYZ_AXIS_COLOR       39
#define L_XYZ_XYPLANE_COLOR    40


/*
 * this definitions here are used for GetDefault(int var)
 */
#define D_SCRF                 0
#define D_BALLF                1
#define D_RODF                 2
#define D_COVF                 3
#define D_ALL                  4
#define D_WFLINEWIDTH          5
#define D_WF3DLINEWIDTH        10005
#define D_OUTLINEWIDTH         10006
#define D_PLLINEWIDTH          6
#define D_PLRADIUS             7
#define D_ATCOL_ONE            8
#define D_ATRAD_SCALE          9
#define D_RCOV_ONE             10004
#define D_ATRAD_ONE            10
#define D_FRAMECOL             11
#define D_FRAMELINEWIDTH       12
#define D_FRAMERODF            13
#define D_MAXSTRUCTSIZE        14
#define D_BACKGROUND           23
#define D_TESSELLATION         24
#define D_UNIBOND              25
#define D_UNIBONDCOLOR         26
#define D_PERSPECTIVE          27
#define D_PERSPECTIVEBACK      28
#define D_PERSPECTIVEFOVY      29
#define D_PERSPECTIVEFRONT     10029
#define D_FOG                  30
#define D_ANTIALIAS            31
#define D_AMBIENT_BY_DIFFUSE   32
#define D_CURRENTFILEFORMAT    33
#define D_HBOND                34
#define D_FORCE_RODTHICKF      35
#define D_FORCE_ARRTHICKF      36
#define D_FORCE_ARRLENF        37
#define D_FORCE_COLOR          38
#define D_XYZ_AXIS_COLOR       39
#define D_XYZ_XYPLANE_COLOR    40



/* this is used just by GetValue */
#define GET_NATOMS       100
#define GET_NAT          101
#define GET_SS_MINX      115
#define GET_SS_MINY      116
#define GET_SS_MINZ      117
#define GET_SS_MAXX      118
#define GET_SS_MAXY      119
#define GET_SS_MAXZ      120
#define GET_AT_MINX      10115
#define GET_AT_MINY      10116
#define GET_AT_MINZ      10117
#define GET_AT_MAXX      10118
#define GET_AT_MAXY      10119
#define GET_AT_MAXZ      10120
#define GET_ATOMLABEL_LABEL       121
#define GET_ATOMLABEL_BRIGHTCOLOR 122
#define GET_ATOMLABEL_DARKCOLOR   123
#define GET_ATOMLABEL_DO_DISPLAY  124
#define GET_ATOMLABEL_ALL_ID      125 /* flag for the ID's of all atoms than have cutom-labels */
#define GET_GLOBALATOMLABEL_BRIGHTCOLOR 126
#define GET_GLOBALATOMLABEL_DARKCOLOR   127
#define GET_GLOBALATOMLABEL_DO_DISPLAY  128

#define GET_FOG_COLORMODE   130
#define GET_FOG_COLOR       131
#define GET_FOG_MODE        132
#define GET_FOG_DENSITY     133
#define GET_FOG_ORT_START_F 134
#define GET_FOG_ORT_END_F   135
#define GET_FOG_PERSP_F1    136
#define GET_FOG_PERSP_F2    137

#define GET_ANTIALIAS_DEGREE 140
#define GET_ANTIALIAS_OFFSET 141

#define GET_ALAT             150


#define SET_ATOMLABEL_DO_DISPLAY       200
#define SET_GLOBALATOMLABEL_DO_DISPLAY 201
#define SET_DO_NOT_DISPLAY_ATOMLABEL   202 /* apply to all atomic lables global and custom */

#define SET_FOG_COLORMODE   210
#define SET_FOG_COLOR       211
#define SET_FOG_MODE        212
#define SET_FOG_DENSITY     213
#define SET_FOG_ORT_START_F 214
#define SET_FOG_ORT_END_F   215
#define SET_FOG_PERSP_F1    216
#define SET_FOG_PERSP_F2    217

#define SET_ANTIALIAS_DEGREE 220
#define SET_ANTIALIAS_OFFSET 221

/* ========================================================================= */
/* "defines" for xcMaybeDelete3DLists & ReDisplay--------------------------- */
#define DELETE3D_NONE         0
#define DELETE3D_BALL         1
#define DELETE3D_SPACE        2
#define DELETE3D_BOND         4
#define DELETE3D_FRAME        8
#define DELETE3D_LINEFRAME   16
#define DELETE3D_BALL_LABEL  32
#define DELETE3D_SPACE_LABEL 64
#define DELETE3D_ALL        127

#define REMAKE3D_NONE         0
#define REMAKE3D_BALL         1
#define REMAKE3D_SPACE        2
#define REMAKE3D_BOND         4
#define REMAKE3D_FRAME        8
#define REMAKE3D_LINEFRAME   16
#define REMAKE3D_BALL_LABEL  32
#define REMAKE3D_SPACE_LABEL 64
#define REMAKE3D_ALL        127


/* --------------------------- */
/* --- some default values --- */
#define DEF_ZOOM  1.7
#define DEF_SCRF  1.5
#define DEF_BALLF 0.4
#define DEF_RODF  0.6
#define DEF_LineWidthWF   1
#define DEF_LineWidthWF3D 1
#define DEF_LineWidthOUT  1
#define DEF_LineWidthPL   1
#define DEF_RadiusPL      6
#define DEF_RCOVF         1.05     /* rcovdef is scaled by this factor and
				      this is criteria for chemmical 
				      connectivity */
#define DEF_ATRADF          1.40
#define DEF_FRAMEWIDTH      1
#define DEF_FRAMERODF       0.1
#define DEF_TESSELLATION    30.0
#define DEF_UNIBOND         0
#define DEF_PERSPECTIVE     0
#define DEF_PERSPECTIVEFOVY  2.5
#define DEF_PERSPECTIVEBACK  3.0 
#define DEF_PERSPECTIVEFRONT 1.0
#define DEF_FOG             0
#define DEF_ANTIALIAS       0
#define DEF_HBOND           0

/* ---------- TRANSLATION & ZOOM TRANSFORMATION ------- */
typedef struct {
  double xtransl;
  double ytransl;
  double zoom;
  float  rotx;
  float  roty;
  float  rotz;
  int    xrotold; /* this is for xc_rot2; B1-Motion rotation */
  int    yrotold; 
  int    trXold; /* this is for xc_transl2; B2-Motion translation */
  int    trYold;
  int    b1motion; /* do we have B1-Motion; boolean */
  int    b2motion; /* do we have B2-Motion; boolean */
  int    shiftB1motion;
} Transform;
Transform tr;


typedef struct {
  int *sqn;            /* sequential number of an atom (ex. atom n. 14) */
  int *sqnat;          /* sequential number of an atom of a kind (ex. atom n. O3) */
  GLfloat (*col)[3];   /* color for each atom; each displayed atom can have
			  its own color */
  int natomkind;
  int atomkind[MAXNAT+1];  /* array of elements that are in structure */
} AtomAtrib;
AtomAtrib atm;

/* this is for flags - to recognise objects */
#define ATOM    0 
#define SELATOM 1
#define BOND    2
#define SELLINE 3
#define SELBOND 4
#define FRAME   5

/* when we orientate objects some objects must be above the others in case,
 * that they are in the same plane
 */
#define Z_OFFSET(x)       ( (0.0001 * x) )

#define XC_NONE  0
#define XC_BOND  1
#define XC_ATOM  2
#define XC_LABEL 2   /* XC_ATOM and XC_LABEL are synonyms */
#define XC_FRAME 4
#define XC_HBOND 5

/* ========================================================================
 * --- AtomBond -- when we open structure, all data about structure goes in
 * AtomBond 
 */
typedef struct {
  int flag;    /* flag for ATOM or BOND or FRAME*/
  int nat;     /* corresponding atomic number 
		* from NAT all other atribbutes are assigned */
  int sqn;     /* sequential number of atom; from sqn other attribbutes are
		  assigned */
  int bondend;
  GLuint list1; /* identifier for OUTLINEBALL list to use when 
		   rendering balls */
  GLuint list2; /* identifier for BALL list to use when 
		   rendering balls */
  double x1;   
  double y1;   /* for ATOMS only X1 Y1 Z1 are needed */
  double z1;
  double x2;   /* for BONDS also X2 Y2 Z2 are neede */
  double y2;
  double z2;
} AtomBond;
AtomBond *coor;  /* coordinates in original units 
		  * on this coordinates I perform rotations,
		  * but never do Z-sorting;
		  * first goes nbonds, than natoms
		  */

/* for Z-orientation */
double *zorient;  /* for ATOMs z == z1
		   * for BONDs z == 0.5 * (z1 + z2) + little
		   * z is center of object for Z coord. */

/*****************************************************************************/
/* someday make the structure from that */
double rcov[MAXNAT + 1];  /* covalent radii */
extern double rvdw[MAXNAT + 1];  /* Van der Waals radii */
double atrad[MAXNAT + 1]; /* radii used to render atoms, balls; it's possible
			     to use rcov[], rvdw[], or some custom radii */
double rball[MAXNAT + 1]; /* balls radii: rball[]= VPf.ballf * atrad[]; */
double rrod;   /* rod's radius rrod= VPf.ballf * VPf.rodf * rcov[H] */
double rframe; /* crystal frame's radius */
extern GLfloat DefAtCol[MAXNAT + 1][3]; /* default atomic color */
float atcol[MAXNAT + 1][3];   /* atomic color used for rendering; it's possible
			      to define a custom colors or use default one */
extern GLfloat DefFrameCol[3];
float framecol[3];
extern GLclampf DefBg[4], DefUnibondCol[4];
GLclampf bg[4], unibondCol[4];


/*****************************************************************************/
/* THIS IS to manege 2D/3D */
#define XC_2D 0
#define XC_3D 1
GLboolean dimType;

/* ========================================================================= */
/* Options; various boolean flags, for managing display options ------------ */
typedef struct {
  GLboolean pipemode;
  GLboolean stickmode; /* TK_r */
  GLboolean ballmode; /* TK_b */
  GLboolean spacefillmode; /* TK_f */
  GLboolean smooth; /* if (smooth) -> GL_SMOOTH else GL_FLAT*/
  GLboolean solid; /* if (solidmode) -> SOLID else WIRE */
  GLboolean lineframe; /* if true LineFrameList is displyed, else
			  either SolidFrameList, either WireFrameList */
  GLboolean anaglyph; /* EV if (anglyphmode) -> ANAGLYPH else NULL  */
  GLboolean stereo;   /* GB do we use the SGI stereo mode */
} Options3D; 


/*****************************************************************************/
/* this is for max size in each directions -> so that we can efficiently 
 * determine glOrtho
 */

typedef struct {
  double x; /* maximum x */ 
  double y;
  double z;
  double r; /* maximum radius */
} MaxSize;

MaxSize max, min;

/* ========================================================================= */
/* this is base for atom-label lists */
GLuint atomlabelOffset; 
GLuint xyzlabelOffset; 

/* ========================================================================= 
 * definition of voriuos vectors; coordinate sistem base-vectors, direct
 * vectors, reciprocal vectors, ...
   ========================================================================= */
#define CRDS_SIZE 130
typedef struct {
  double crdmajor[4][4]; /* coor-sist major vec */
  double crdnew[4][4];   /* coor-sist new vec   */
  double crdold[4][4];   /* coor-sist old vec   */
  double crdvec[16];
  double prim[4][4];     /* primitiv-vectors; if we don't have crystal, then
                          * expample: for slabs, the third vector is 
                          * (0,0,max.z)
                          */
  double conv[4][4];     /* conventional-vectors; if not crystal do same as
                          * with prim[][]
			  */
  double recprim[4][4];  /* reciprocal primitive-vectors */
  double recconv[4][4];  /* reciprocal conventional-vectors */
} Vector;
Vector vec;


/* ========================================================================= */
/* structure for tracking the dimension of orthographic projection
 */
typedef struct {
  double size;
  double minx;
  double maxx;
  double miny;
  double maxy;
  double minz;
  double maxz;
  double fi, near, far; /* for PERSPECTIVE */
} OrthoProj;


/*****************************************************************************/
/* this is for file format !!!!! */
#define FORMAT_NONE       0
#define FORMAT_XSF        1
#define FORMAT_XYZ        2
#define FORMAT_PDB        3

int current_file_format;

/*****************************************************************************/
/* XC_CPP_* are pre-processor flags */
#ifndef XC_CPP_BOOLEAN
typedef int boolean;
#define XC_CPP_BOOLEAN
#endif


/*****************************************************************************/
/* this is structures that hold XCRYSDEN file format information             */
/* XCRYSDEN file-format is defined and made by GENGEOM program               */
#define XCR_CELL           0   /* cell is the unit of repetition; look below */
#define XCR_TR_ASYM        1   /* translational symetric part is unit of     
                                  repetition                                 */
#define XCR_NOCELL         0
#define XCR_PRIMCELL       1   /* primitiv cell                              */
#define XCR_CONVCELL       2   /* conventional cell                          */
#define XCR_PARAPIPEDAL    0   /* parapipedal shape of cell                  */
#define XCR_HEXAGONAL      1   /* hexagonal shape of cell                    */
#define XCR_BZ             2

/* this must be the same as in gengeom.f */
#define XCR_CELL_PC   1
#define XCR_CELL_AC   2
#define XCR_CELL_BC   3
#define XCR_CELL_CC   4
#define XCR_CELL_FC   5
#define XCR_CELL_IC   6
#define XCR_CELL_RC   7
#define XCR_CELL_HC   8
#define XCR_CELL_TNRC 9

#define XCR_NCELL_TYPES 10

#define XCR_CELL_IPC   1
#define XCR_CELL_IAC   2
#define XCR_CELL_IBC   2
#define XCR_CELL_ICC   2
#define XCR_CELL_IFC   4
#define XCR_CELL_IIC   2
#define XCR_CELL_IRC   3
#define XCR_CELL_IHC   3
#define XCR_CELL_ITNRC 3

typedef struct {
  int     dim;       /* dimensionality of sistem                             */
  int     groupn;    /* number of group; it's from c95's ftn34;
                        actually it's family not group                       */
  float   cellpos[XCR_NCELL_TYPES][4][3]; /* fractional positions within the cell */
  int     npos[XCR_NCELL_TYPES]; /* number of atoms(positions) within desired
				    cell type                                */
  boolean ldimgroup; /* true if group-number is present in xcr-file          */
  boolean lconvvec;  /* true if conv vectors are specified in xcr-file       */
  boolean lprimvec;  /* true if prim vectors are specified in xcr-file       */
  boolean lrecprimvec; /* true if prim reciprocal vectors are specified      */
  boolean lrecconvvec; /* true if conv reciprocal vectors are specified      */
  boolean lprimcoor; /* true if primcoord are specified in xcr-file          */
  boolean lconvcoor; /* true if convcoord are specified in xcr-file          */
  boolean lprimwigner; /* true if primitive Wigner-Seitz cell is specified   */
  boolean lconvwigner; /* true if conventional Wigner-Seitz cell is specified*/
  boolean lprimbz;   /* true if primitive Brillouin zone is specified        */
  boolean lconvbz;   /* true if conventional Brillouin zone is specified     */
  boolean ldatagrid2D; /* true if DATAGRID2D is specified                    */
  boolean ldatagrid3D; /* true if DATAGRID3D is specefied                    */
  boolean lbandgrid2D; /* true if BANDGRID2D is specified                    */
  boolean lbandgrid3D; /* true if BANDGRID3D is specefied                    */
  boolean lforce;      /* true if forces were specified                      */
  /*boolean lHbond;*/      /* true if hydrogen bonds are present                 */
  int natr;          /* number of atoms in primitiv cell                     */
  int ncell;         /* "atoms per cell", FC ->4, RC -> 3, etc.              */
  int     nunit[3];  /* number of units in each directions                   */
  int     unit;      /* what is the unit "cell" or translational asym. part  */
  int     celltype;  /* primitiv or convetional                              */
  int     shape;     /* parapipedal or hexagonal                             */
  int     *prim_nat; /* atomic numbers of atoms in primcell */
  int     *conv_nat; /* atomic numbers of atoms in primcell */
  double  (*prim_coor)[3];  /* Cartesian coordinates of atoms in primcell (in ANGSTROMS) */
  double  (*prim_fcoor)[3]; /* fractional coordinates of atoms in primcell  */
  double  (*prim_forc)[3];  /* Cartesian components of atomic forces in primcell */
  boolean prim_lforc;       /* are forces in PRIMCOORD present */
  double  (*conv_fcoor)[3];  /* fractional coordinates of atoms in convcell  */
  double  (*conv_coor)[3];  /* Cartesian coordinates of atoms in convcell (in ANGSTROMS) */  
} XcrInfo;
XcrInfo xcr;


/*****************************************************************************/
/* here are display-attributes for isosurface and colorplane                 */
typedef struct {
  int       drawstyle;   /* drawstyle for isosurface rendering
                       * ISOSURF_WIRE or ISOSURF_SOLID      */
  GLenum    shademodel;
  int       transparent;
  GLboolean lighting;    /* lighting NO/OFF */
} ISO_ATTRIB;


/*****************************************************************************/
/* this is structure used for evauating custom Tcl/tk command's options      */
#ifdef XC_CPP_GLPARAM  /* xcGLparam.h header file was included               */
#   define COM_ISOEXPAND_REPEAT_DEFAULT         0
#   define COM_ISOEXPAND_REPEAT_CONVCELL        1
#   define COM_ISOEXPAND_REPEAT_PRIMCELL        2
#   define COM_ISOEXPAND_SHAPE_DEFAULT          0
#   define COM_ISOEXPAND_SHAPE_PARAPIPEDAL      1
#   define COM_ISOEXPAND_SHAPE_HEXAGONAL        2
#   define COM_ISOEXPAND_EXPAND_WHOLE           0
#   define COM_ISOEXPAND_EXPAND_LIST            1
    typedef GetGlParam GetComOption;
#endif

/*****************************************************************************/
#define WIGNERSEITZ_MAXPOLY    26 /* hope this is enough */
#define WIGNERSEITZ_MAXVERTEX  15
typedef struct {
  int   npoly; /* number of polygons */
  int   nvert[WIGNERSEITZ_MAXPOLY];
  float max;   /* maximum value of coordinates (either X,Y,Z) */
  float poly[WIGNERSEITZ_MAXPOLY][WIGNERSEITZ_MAXVERTEX][3]; 
  float norm[WIGNERSEITZ_MAXPOLY][WIGNERSEITZ_MAXVERTEX][3]; 
} WIGNERSEITZ;
WIGNERSEITZ wsp, wsc;


/*****************************************************************************/
/*        LINESTIPPLE --- LINESTIPPLE --- LINESTIPPLE --- LINESTIPPLE        */
/*****************************************************************************/
#define LINESTIPPLE0      0x5555
#define LINESTIPPLE1   	  0x6b2d
#define LINESTIPPLE2   	  0x6735
#define LINESTIPPLE3   	  0xcf5e
#define LINESTIPPLE4   	  0xaf57
#define LINESTIPPLE5   	  0x3f49
#define LINESTIPPLE6   	  0xef4c
#define LINESTIPPLE7   	  0xe73a
#define LINESTIPPLE8   	  0x4f12
#define LINESTIPPLE9   	  0x8f32
#define LINESTIPPLE_SOLID 0xffff

#define LINESTIPPLE_FACTOR0  1
#define LINESTIPPLE_FACTOR1  2
#define LINESTIPPLE_FACTOR2  3
#define LINESTIPPLE_FACTOR3  4

/*************************/
/* DATAGRID --- DATAGRID */
/*************************/
#define DATAGRID_2D       2
#define DATAGRID_3D       3
#define RECIP_DATAGRID_3D 4
#define DATAGRID_MAXSUBINDEX 5000 /* do we need more that 5000 bands for Fermi surface ??? */

#ifndef XC_CPP_NO_STDIO /* load DATAGRID only if stdio.h was included */
struct DATAGRID {
  FILE   *fp;
  int    type;        /* type of DATAGRID, i.e. 2D or 3D */
  int    index;       /* identifier for *grid */
  char   *ident;     /* identifier for the user */
  int    n_of_subgrids; /* number of subgrids (up to DATAGRID_MAXSUBINDEX) */
  int    n[3];        /* number of points in two/three direction */
  float  orig[3];   /* xyz of origin */
  float  vec[3][3]; /* format of vec: [vec#][xyz] */
  struct DATAGRID *ptr;   /* pointer to previous struct */
  int    selected[DATAGRID_MAXSUBINDEX];
  long   fpos[DATAGRID_MAXSUBINDEX];       /* position if the file */
  float  signfactor[DATAGRID_MAXSUBINDEX]; /* factor for subindex's datagrid */
  char   *subident[DATAGRID_MAXSUBINDEX];
  int    lband, nband, *band_index[DATAGRID_MAXSUBINDEX];
  long   *band_fpos[DATAGRID_MAXSUBINDEX];
  float  minvalue[DATAGRID_MAXSUBINDEX]; /* minimum value in the grid */
  float  maxvalue[DATAGRID_MAXSUBINDEX]; /* maximum value in the grid */
};
#endif

/* here are some data about raster font stored in xcLabels.c */
typedef struct {
  GLsizei wid;
  GLsizei height;
  GLfloat w2;
  GLfloat h2;
  GLfloat wp2;
  GLfloat hp2;
} RasterFontSize;


#ifdef _TK
typedef struct {
  GLuint  base;
  GLsizei width;
  GLsizei height;
  GLfloat bright_color[3];
  GLfloat dark_color[3];  
  short   do_display;
  char    *label;
  Tk_Font tkfont;
} AtomicLabel;
#endif

/* minimum and maximum coorinates of structure atoms */
typedef struct {
  double minX;
  double minY;
  double minZ;
  double maxX;
  double maxY;
  double maxZ;
} StructSize;
  

/*
  Now, just ".mesa" togl is over. We can create an arbitrary nember of togls
  and display the structures into them. Everyting is dynmicaly allocated, and
  NEW_WIN_CONTEXT is a structure that holds everything !!!
*/
#ifdef TOGL_H

#define FS_SINGLE  0
#define FS_MULTI   1

struct FS_CONTEXT {
  int           ntogl;
  struct Togl **toglVector;
};

struct NEW_WIN_CONTEXTstr {
  struct NEW_WIN_CONTEXTstr *ptr;
  struct Togl               *togl;
  struct FS_CONTEXT         fermiContext;
  int           index;
  StructPar     VPf;
  ModelPar      MVf;
  Transform     tr;
  StructSize    ss;
  Vector        vec;
  void          *recprim_cage; /* pointer to CellCage recprim_cage; */
  GLclampf      bg[4];
  void          (*xcDisplay)(struct Togl *togl);
};

typedef struct NEW_WIN_CONTEXTstr NEW_WIN_CONTEXT;
#endif

#define XC_FORCE_TRESHHOLD    0.0005
#define XC_FORCE_LENGTHFACTOR 200

/* this is for rendering Force Vectors */
typedef struct {
  float  (*ScaleFunc)(float value);
  float  threshold;
  float  lengthfactor;

  float  rod_thickf;     /* vector's thickness factor */
  float  arr_thickf; /* thickness factor for the arrow */
  float  arr_lenf;   /* length of the arrow */

  float  color[4];

  /* 
     this is maximum size for atom+force (for making correct
     projection/viewport) 
  */
  float max_size; 
} ForceVector;


/* structure for gluPerspective */
typedef struct {
  GLdouble fovy;
  GLdouble aspect;
  GLdouble near;
  GLdouble far;
  GLdouble shiftZ;
} PERSPECTIVE;

/* fog (i.e. depth cuing) */
#define XC_FOG_BGCOLOR     0
#define XC_FOG_CUSTOMCOLOR 1
#define XC_FOG_LINEAR      2
#define XC_FOG_EXP         3
#define XC_FOG_EXP2        4

typedef struct {
  int     colormode; /* can be CUSTOMCOLOR or BGCOLOR */
  GLfloat color[4];
  GLint   mode;
  GLfloat density;
  GLfloat ort_start_f;
  GLfloat ort_end_f;
  GLfloat persp_f1;
  GLfloat persp_f2;
} XCfog;

typedef struct {
  int    degree;
  float  offset;
} XCantialias;


/* 3D vector:     double precision */
typedef double VEC3d[3];        

/* H-bonds */
typedef struct {
  int   n;                      /* number of H-bonds                         */
  int   max_n;                  /* allocated number of H-bonds               */
  int   *H_like_list;           /* lists of H-like atoms                     */
  int   *O_like_list;           /* lists of O-like atoms                     */
  float  color[4];              /* color of the H-bonds                      */
  double length_min;            /* minimum H-bond length                     */
  double length_max;            /* maximum H-bond length                     */  
  double angle_min;             /* minimum-angle; Namely the angle
				   between H-chemical bond and H-bond
				   must be close to 180 degrees
				   (i.e. greater than angle_min)             */
  unsigned short line_pattern;  /* line stipple pattern                      */
  int    line_patternsize;            
  float  line_width;            
  VEC3d  *start_coor;           /* start coordinates of the H-bond           */
  VEC3d  *end_coor;             /* end coordinates of the H-bond             */
} H_Bond;


int crystal_version; /* what is the version of CRYSTAL program we are using */

/*
  color for the Coordinate system display
*/
typedef struct {
  float  axis_color[4];
  float  xyplane_color[4];
} XYZ_Attrib;


/*
  movie making
*/

#define MOVIE_MODE_EVERY_SNAPSHOT 0
#define MOVIE_MODE_REALTIME_INTERVAL 1

typedef struct {
  int mode;
  int nframe;
  int doit;
  int printing;
  char *dir;
} realTimeMove;

extern const char *printImage;

#endif
