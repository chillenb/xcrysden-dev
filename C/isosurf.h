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
 * Source: $XCRYSDEN_TOPDIR/C/isosurf.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#ifndef ISOSURF_H
#define ISOSURF_H
 
#include <stdio.h>
#include <string.h>
#include <GL/gl.h>

#define XC_CPP_ISOSURF 

#define MAX_VERTEX_COORDINATION 20


#define MAX_ISOSURFLEVELS  2   /* positive & negative isolevel */
#define ISOLINE_MAXLEVEL   500

#define ISODATA_MAXSTACK    4
#define ISODATA_MAXFRAME    1000
#define ISODATA_MAXFILES    5     /* maximum number of files, that can be 
                                   * merged together; this is alo maximum 
                                   * number of frames in stack #3   */
#define ISODATA_MAXBLOCK    1000  /* block is something like frame in stack 
                                   * order 0 */

#define ISOSURF_SHADE_FLAT   GL_FLAT    /* flat-shading of isosurface */
#define ISOSURF_SHADE_SMOOTH GL_SMOOTH  /* smooth-shading of isosurface */
#define ISOSURF_WIRE         0     /* render isosurface as a wireframe     */
#define ISOSURF_SOLID        1     /* render isosurface in solid fashion   */
#define ISOSURF_DOT          2     /* render isosurface as doted           */
#define ISOSURF_CONTOUR      3     /* render isosurface in contour fashion */
#define ISOSURF_TRANSP_OFF   0     /* isosurface transparency is turnd off */
#define ISOSURF_TRANSP_ON    1     /* isosurface transparency is turnd on  */
#define ISOSURF_MARCHING_CUBES 0 /* use marching cubes algorithm */
#define ISOSURF_TETRAHEDRAL    1 /* use tetrahedral decomposition */
#define ISOSURF_NORMALS_GRADIENT 0 /* calculate normals by interpolating grid's gradient vectors */
#define ISOSURF_NORMALS_TRIANGLE 1 /* calaculate normals by averaging touching triangles normals */

#define COLORBASE_FIRST      0  /* to know what is the number of first item */
#define COLORBASE_MONO       0
#define COLORBASE_RAINBOW    1
#define COLORBASE_RGB        2
#define COLORBASE_GEOGRAPHIC 3
#define COLORBASE_BLUE_WHITE_RED 4
#define COLORBASE_BLACK_BROWN_WHITE 5
#define COLORBASE_LAST       5  /* to know what is the number of last item */

#define SCALE_FUNC_FIRST    0
#define SCALE_FUNC_LIN      0
#define SCALE_FUNC_LOG      1
#define SCALE_FUNC_LOG10    2
#define SCALE_FUNC_SQRT     3
#define SCALE_FUNC_ROOT3    4
#define SCALE_FUNC_GAUSS    5
#define SCALE_FUNC_SLATER   6
#define SCALE_FUNC_LAST     6

/*****************************************************************************/
/* constants for ReadBlock0 function                                         */
#define RB0_INIT   0      /* initialization */
#define RB0_FIND   1      /* find out the number of frames in stack #0 */
#define RB0_DUMMY  2      /* dummy reading  */
#define RB0_READ   3      /* read data      */
/****************************************************************************/

#define ISO_NULL        0
#define ISO_INIT        1
#define ISO_STACK       2
#define ISO_SIGN        4
#define ISO_FILES       8
#define ISO_POINTS      16
#define ISO_DATA        32
#define ISO_ISOLEVEL    64
#define ISO_POLYGONISE  128
#define ISO_COLORPLANE  256
#define ISO_COLORPLANE1 512
#define ISO_COLORPLANE2 1024
#define ISO_COLORPLANE3 2048
#define ISO_ISOLINE     4096
#define ISO_ISOLINE1    8192
#define ISO_ISOLINE2    16348
#define ISO_ISOLINE3    32768
#define ISO_EXPAND      65536

/*****************************************************************************/
/* this is for isosurf's Material properties                                 */
#define MAT_DEFAULT    1
#define MAT_CUSTOM     2
#define MAT_ONELEVEL   4
#define MAT_POSLEVEL   8
#define MAT_NEGLEVEL   16
#define MAT_ISOSURF    32
#define MAT_COLORPLANE 64
#define MAT_SOLID   128

/*****************************************************************************/
/* here is isosurface extraction specific stuff                              */
#define XC_CPP_XYZ
typedef struct {
  float x;
  float y;
  float z;
} XYZ;

#define XYZ_det3x3( v0, v1, v2 ) \
  ( ( v0.x*v1.y*v2.z + v0.y*v1.z*v2.x + v0.z*v1.x*v2.y - \
      v2.x*v1.y*v0.z - v2.y*v1.z*v0.x - v2.z*v1.x*v0.y ) )

typedef struct {
  XYZ pos;
  XYZ nml;
} VERTEX;

typedef struct {
  XYZ   p;
  XYZ   n;
  float val;
} GRIDVERTEX;

typedef struct {
  XYZ p[3];    /* Triangles Data */
} TRIANGLE;

typedef struct {
  int      numT[MAX_ISOSURFLEVELS];          /* Number of Triangles */
  int      numV[MAX_ISOSURFLEVELS];          /* Number of Vertices */
} SURFACE;

typedef struct {
  XYZ    p[8];
  float val[8];
} GRIDCELL;
typedef struct {
  float p[8][3];
  float val[8];
} GRIDCELL3;

typedef struct {
  int hexa;
  int deci;
} CASE;

typedef struct {
  int ix;
  int iy;
  int iz;
  int nx;
  int ny;
  int nz;
} GRID;

/* date Fri Apr  7 14:27:00 CEST 2000
   ----------------------------------
   MAJOR reprograming of polygonalization !!!!
   these are new tydefs:
*/
typedef struct {
  int nt;  /* number of triangles in grid-cube */
  int t0;  /* index of 0-th triangle in grid-cube */
} GRID_INFO;
typedef struct {
  int i,j,k;          /* cube indices */
  int ver_ind[3];     /* normal index */
  int ver_status[3];  /* status of individual vertices in triangle */
} TRIG_INFO;
typedef struct {
  int t; /* which triangle among all the triangles */
  int v; /* 1st/2nd/3rd vertex in a "t" triangle */
} VERT2TRIG_INFO;
typedef struct {
  unsigned int i,j,k;
} IJK;
struct ISOSURFACEstruct {
  struct ISOSURFACEstruct *ptr;
  int             index;
  /*-------------------------*/
  int             revertnormals;
  int             nvertex, ntriangl;
  int             smooth_nstep;
  float           smooth_weight;
  int             *triangl_status;
  XYZ             *vertex;
  XYZ             *vertex_orig; 
  XYZ             *normal;
  float           (*color)[4];
  TRIG_INFO       *tri2verIN;
  int             *nver2triIN;
  VERT2TRIG_INFO (*ver2triIN)[MAX_VERTEX_COORDINATION];
};
typedef struct ISOSURFACEstruct ISOSURFACE;

/*VERTEX     *vertex[MAX_ISOSURFLEVELS];*/
XYZ        *vertex[MAX_ISOSURFLEVELS];
XYZ        *normal[MAX_ISOSURFLEVELS];
TRIG_INFO  *tri2verIN[MAX_ISOSURFLEVELS];
TRIANGLE   *triangles[MAX_ISOSURFLEVELS];
SURFACE    surf;
CASE       cs;
GRID       grd, newgrd; /* grd is original as read from file, whereas
			   newgrd is used for interpolation --> thatwhy
			   newgrd is used in isorender.c && colorplane.c &&
			   polygonise.c */
GRIDVERTEX ***gridvertex;
/*****************************************************************************/


/*****************************************************************************/
/* here is stuff related to "colorplane"                                     */
typedef struct {
  int ncol;    /* up to 15 different basic colors */
  float basecol[15][4];
  float baseintv[15];
} BASECOLOR;
  

typedef struct {
  float p[3];
  float val;
  float col[4];
} PLANEVERTEX;



/*****************************************************************************
 * here is stuff related to isosurface && Tcl?tk implementation of xc_iso 
 * group of custom Tcl command
 */
#ifdef XC_CPP_STRUCT
typedef struct {
  int    nstack; /* n of stack used; 
                  * maximum nstack could be ISODATA_MAXSTACK */
  int    nframe[ISODATA_MAXFILES][ISODATA_MAXSTACK];    /* number of frames in
							 * each stack        */
  int    framesign[ISODATA_MAXSTACK][ISODATA_MAXFRAME]; /* sign of stack     */
  int    dim[MAX_ISOOBJECTS];             /* dimension of objects            */
  float  min;                             /* minimum vertex value            */
  float  max;                             /* maximim vertex value            */
  float  mintol;                          /* minimum tolerances used is 
					     VertexInterp */
  float  min_allowed[MAX_ISOOBJECTS];     /* minimum vertex value to be drawn*/
  float  max_allowed[MAX_ISOOBJECTS];     /* maximum vertex value to be drawn*/
  char   *isofiles[ISODATA_MAXFILES];     /* filenames                       */
  char   error[512];                      /* error message                   */
  FILE   *fp;                             /* filepointer                     */
  FILE   *bin_fp;                         /* filepointer for binary file     */
  FILE   *bin_vertex_fp;                  /* filepointer for binary file where 
					     current vertex / planevertex data
					     are stored                      */
  int    nlevel;                          /* number of levels for isosurface's
					     isolevel; max # of levels:
					     MAX_ISOSURFLEVELS               */
  float  isolevel[MAX_ISOSURFLEVELS];
  double points[MAX_ISOOBJECTS][4][3];
  XYZ    vec[MAX_ISOOBJECTS][3];          /* vectors that specify space of isosurface/isoplate */
  XYZ    colnml[MAX_ISOOBJECTS];          /* normal for colorplane           */
  int    smoothsteps;                     /* # of surface smoothing steps    */
  float  smoothweight;                    /* surface smoothing weight        */
  int    cell_orientation;                /* orientation of the unit-cell spanned by vecs */
} ISOSURF_DATA;
ISOSURF_DATA isodata;

typedef int logical;  /* synonym for logical type */
typedef struct {
  logical       gridvertex_malloc;
  logical       plvertex_malloc[MAX_ISOOBJECTS];
  logical       isoline2D_malloc[MAX_ISOOBJECTS];
  logical       vertex_malloc[MAX_ISOSURFLEVELS];
  logical       triangl_malloc[MAX_ISOSURFLEVELS];
  logical       tri2verIN_malloc[MAX_ISOSURFLEVELS];
  int           max_n_vertex[MAX_ISOSURFLEVELS];
  int           max_n_triangl[MAX_ISOSURFLEVELS];
  int           max_n_isoline2D[MAX_ISOOBJECTS][ISOLINE_MAXLEVEL];
  logical       bin_file_open;
  logical       bin_vertex_file_open;
  unsigned int  stateflag;
} ISOSTATE;

/* float (*datavalue)[MAXSIZE][MAXSIZE]; */

/****************************************************************************/
/* ---> this are poniters for storing ADDRESSES                             */
/* void   *gridvertex_address; */          /* memory address for gridvertex &
                                              plvertex                      */
/* void   *vertex_address; */              /* memory address for vertex     */
/* void   *triangl_address; */             /* memory address for triangles  */


/*****************************************************************************/
/* this is used for xc_isoexpand                                             */
typedef struct {
  double rep_vec[4][4]; /* here it's double, because it should be the same 
                           as in Vector structure                            */
  int    shape;
  int    expand;
  int    transl[3][3];  /* 3 translational vectors in fractional coordinates; 
                           used for isodata.dim = 2                          */
  int    irepvec[3];    /* how many times to repeat "iso*" in each direction */
} IsoExpand;

IsoExpand isoexpand[MAX_ISOOBJECTS];
#endif


/*****************************************************************************/
/* ISOLINES --- ISOLINES --- ISOLINES --- ISOLINES --- ISOLINES --- ISOLINES */
/*****************************************************************************/
#define ISOLINE_MONOCOLOR   0
#define ISOLINE_POLYCOLOR   1
#define ISOLINE_NODASH      2
#define ISOLINE_NEGDASH     3
#define ISOLINE_FULLDASH    4
typedef struct {
  XYZ p1;
  XYZ p2;
} LINE;

#ifdef XC_CPP_STRUCT
typedef struct {
  int    linecolor;                      /* code for linecolor type          */
  float  monocolor[4];                   /* monocolor                        */
  unsigned short linedash;               /* code for linedash  type          */
  int    nlevel;                         /* number of levels                 */
  int    iseg[ISOLINE_MAXLEVEL];         /* current # of segment for ilevel  */
  float  level[ISOLINE_MAXLEVEL];        /* values of levels                 */
  unsigned short dash[ISOLINE_MAXLEVEL]; /* types of line for each level     */
  int    dashfactor[ISOLINE_MAXLEVEL];   /* factor for linedash              */
  float  color[ISOLINE_MAXLEVEL][4];     /* color of ilevel contour          */
  float  linewidth;                      /* width of the isoline             */
  LINE   *segment[ISOLINE_MAXLEVEL];     /* segments data                    */
} ISOLINE2D;

ISOLINE2D isoline2D[MAX_ISOOBJECTS];
#endif

#endif
