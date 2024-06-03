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
 * Source: $XCRYSDEN_TOPDIR/C/polygonise.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

/***********************************************/
/* date Fri Apr  7 14:27:00 CEST 2000          */
/* MAJOR reprograming of polygonalization !!!! */
/***********************************************/

/* pri elektronski gostoti pozitivne normale kazejo v notranjost; zato sem dal
 * nml -= mnl_tmp;
 */
#include <stdio.h> 
#include <math.h>
#include <stdlib.h> 
#include "struct.h"
#include "isosurf.h"
#include "memory.h"
#include "xcfunc.h"

#define MINTOL_P 1e-12
#define FLOAT_MINTOL  3.0e-8
#define ASSIGN_GRID2IJK(a,b) a.i = b.ix; a.j = b.iy; a.k = b.iz
#define EQUAL_VERTICES(a,b)    ( (fabs(a.x-b.x) < FLOAT_MINTOL && fabs(a.y-b.y) < FLOAT_MINTOL && fabs(a.z-b.z) < FLOAT_MINTOL) )
#define GET_INDEX(i,j,k,N,M,L) ( (k + L*(j + M*i)) )

#define XYZ_A_SUB_B(dst,a,b) \
do { \
  dst.x = a.x - b.x; \
  dst.y = a.y - b.y; \
  dst.z = a.z - b.z; \
 }  while(0)

extern ISOSTATE isostate;

/*---------------------------------------------------------------------------*/
extern void
VertexRecognise(GRID          grd,
		int           ntriangl,
		TRIANGLE      *triangles,
		int           *triangl_status,
		TRIG_INFO     *tri2verIN,
		XYZ           *vertex,
		int           *vertex_status,
		VERT2TRIG_INFO (*ver2triIN)[MAX_VERTEX_COORDINATION],
		int           *nver2triIN,
		GRID_INFO     *grid_info,
		int           *nvertex);
extern void
SurfSmoothing(int            nstep, 
	      float          weight,
	      XYZ            *vertex,
	      int            *vertex_status,
	      VERT2TRIG_INFO (*ver2triIN)[MAX_VERTEX_COORDINATION],
	      TRIG_INFO      *tri2verIN,
	      int            *nver2triIN,
	      int            nvertex);
extern void
SurfNmlAver(int       ntriangl,
	    int       nvertex,
	    XYZ       *nml,
	    XYZ       *vertex,	    
	    TRIG_INFO *tri2verIN,
	    float     *sign,
	    int       is_tetrahedral);
/*---------------------------------------------------------------------------*/

void MarchingCubes(float isolevel, int ilevel, int smooth_nstep, float smooth_weight, int algorithm, int shade_model, int normals_model);

int Polygonise(GRIDCELL grid,
	       float isolevel,
	       GRID   ijk,
	       TRIANGLE *triangles,
	       TRIG_INFO *trig_info);
XYZ VertexInterp(float isolevel, XYZ p1, XYZ p2, float valp1, float valp2);
TRIG_INFO *xcRealloc2xBiggerTrigInfo(TRIG_INFO *vec, int *size);



static VERT2TRIG_INFO (*ver2triIN_P)[MAX_VERTEX_COORDINATION];
static int            *nver2triIN_P;
static int            *vertex_status_P, *triangl_status_P;
static float          *sign_P;
void
WrapperSurfSmoothing(int ilevel) {

  /* perform smoothing */
  SurfSmoothing(isodata.smoothsteps, isodata.smoothweight, vertex[ilevel], 
		vertex_status_P, ver2triIN_P, tri2verIN[ilevel],
		nver2triIN_P, surf.numV[ilevel]);
  SurfNmlAver(surf.numT[ilevel], surf.numV[ilevel], 
	      normal[ilevel], vertex[ilevel], tri2verIN[ilevel], sign_P, 0);
}


void
MarchingCubes(float isolevel, 
	      int   ilevel, 
	      int   smooth_nstep, 
	      float smooth_weight,
	      int   algorithm,
	      int   shade_model,
	      int   normals_model)
{
  register int   i, ix1, iy1, iz1;
  int            mem_size, ntriangl, index, new_maxelem, ivrt;
  GRIDCELL       grid;
  GRID_INFO      *grid_info;
  TRIANGLE       *triangles;
  VERT2TRIG_INFO (*ver2triIN)[MAX_VERTEX_COORDINATION];
  int            *nver2triIN, nvertex;
  int            *vertex_status, *triangl_status;
  float          *sign=NULL;
  XYZ            spanning_vec[3];

  mem_size = 
    MAX(newgrd.nx, newgrd.ny) * MAX(newgrd.ny, newgrd.nz);

  /* -----------------------------------------------------------------------
   * INITIALIZATIONS  
   * ----------------------------------------------------------------------- */
  surf.numT[ilevel] = 0;
  surf.numV[ilevel] = 0;

  isostate.max_n_triangl[ilevel] = mem_size;
  triangles = (TRIANGLE *) malloc( (size_t) sizeof(TRIANGLE) * isostate.max_n_triangl[ilevel] );
  if ( !triangles ) xcError("allocation error for triangles");
  isostate.triangl_malloc[ilevel] = 1; /* not used any more !!! */
  
  if ( isostate.tri2verIN_malloc[ilevel] ) free((FREE_ARG) tri2verIN[ilevel] );

  tri2verIN[ilevel] =  (TRIG_INFO *) malloc((size_t)mem_size * sizeof(TRIG_INFO));
  isostate.tri2verIN_malloc[ilevel] = 1;
  
  grid_info = (GRID_INFO *) malloc((size_t) (newgrd.nx*newgrd.ny*newgrd.nz) * sizeof(GRID_INFO));
  if ( !tri2verIN[ilevel] || !grid_info) 
    xcError("allocation error in Marching Cubes");

  ntriangl=0;
  ix1 = newgrd.nx - 1;
  iy1 = newgrd.ny - 1;
  iz1 = newgrd.nz - 1;

  /* calculate normals at the grid points */
  
  if ( normals_model == ISOSURF_NORMALS_GRADIENT ) gridNormals (gridvertex, newgrd);

  XYZ_A_SUB_B(spanning_vec[0],  gridvertex[ix1][0][0].p, gridvertex[0][0][0].p);
  XYZ_A_SUB_B(spanning_vec[1],  gridvertex[0][iy1][0].p, gridvertex[0][0][0].p);
  XYZ_A_SUB_B(spanning_vec[2],  gridvertex[0][0][iz1].p, gridvertex[0][0][0].p);
  
  /* do the marching ... */
  for (grd.ix = 0; grd.ix < ix1; grd.ix++) {
    for (grd.iy = 0; grd.iy < iy1; grd.iy++) { 
      for (grd.iz = 0; grd.iz < iz1; grd.iz++) {

	new_maxelem = ( algorithm == ISOSURF_MARCHING_CUBES ? 5 : 12 );

	/* take care of reallocation if necessary !!! */
	if ( ntriangl + new_maxelem >= isostate.max_n_triangl[ilevel] ) 
	  {
	    triangles         = xcRealloc2xBiggerTriangl( triangles, &isostate.max_n_triangl[ilevel] );
	    tri2verIN[ilevel] = xcRealloc2xBiggerTrigInfo( tri2verIN[ilevel], &mem_size );
	  }
	
	grid.p[0]   = gridvertex[grd.ix  ][grd.iy  ][grd.iz  ].p;
	grid.val[0] = gridvertex[grd.ix  ][grd.iy  ][grd.iz  ].val; 
	
	grid.p[1]   = gridvertex[grd.ix  ][grd.iy+1][grd.iz  ].p;
	grid.val[1] = gridvertex[grd.ix  ][grd.iy+1][grd.iz  ].val;

	grid.p[2]   = gridvertex[grd.ix+1][grd.iy+1][grd.iz  ].p;
	grid.val[2] = gridvertex[grd.ix+1][grd.iy+1][grd.iz  ].val;

	grid.p[3]   = gridvertex[grd.ix+1][grd.iy  ][grd.iz  ].p;
	grid.val[3] = gridvertex[grd.ix+1][grd.iy  ][grd.iz  ].val;

	grid.p[4]   = gridvertex[grd.ix  ][grd.iy  ][grd.iz+1].p;
	grid.val[4] = gridvertex[grd.ix  ][grd.iy  ][grd.iz+1].val;

	grid.p[5]   = gridvertex[grd.ix  ][grd.iy+1][grd.iz+1].p;
	grid.val[5] = gridvertex[grd.ix  ][grd.iy+1][grd.iz+1].val;

	grid.p[6]   = gridvertex[grd.ix+1][grd.iy+1][grd.iz+1].p;
	grid.val[6] = gridvertex[grd.ix+1][grd.iy+1][grd.iz+1].val;

	grid.p[7]   = gridvertex[grd.ix+1][grd.iy  ][grd.iz+1].p;
	grid.val[7] = gridvertex[grd.ix+1][grd.iy  ][grd.iz+1].val;

	index=GET_INDEX(grd.ix,grd.iy,grd.iz, ix1, iy1, iz1);

	grid_info[index].t0=ntriangl;

	if ( algorithm == ISOSURF_MARCHING_CUBES ) 
	  {
	    /* ===> POLYGONISE: MarchingCubes <=== */
	    ntriangl += Polygonise(grid, isolevel, grd, &triangles[ntriangl], 
				   (TRIG_INFO *) &tri2verIN[ilevel][ntriangl]);
	  }
	else
	  {
	    /* ===> POLYGONISE: TetrahedralDecomposition <=== */
	    ntriangl += PolygoniseTetrahedral(grid, isolevel, grd, &triangles[ntriangl], 
					      (TRIG_INFO *) &tri2verIN[ilevel][ntriangl]);
	  }
	
	/* TESTING ... */
 	grid_info[index].nt = ntriangl - grid_info[index].t0; 
 	/* grid_info[index].nt = ntriangl - grid_info[index].t0 + 1; */
      } /* iz */
    }
  }

  surf.numT[ilevel] = ntriangl;
  printf("NUMBER OF TRIANGLES: %d\n",surf.numT[ilevel]);

  /* malloc space for *vertex & *triangles */
  if ( isostate.vertex_malloc[ilevel]  &&  isostate.max_n_vertex[ilevel] < ntriangl*3 ) {
    free( (FREE_ARG) vertex[ilevel] );
    free( (FREE_ARG) normal[ilevel] );
  }
  
  isostate.max_n_vertex[ilevel] = ntriangl * 3;

  vertex[ilevel] = (XYZ *) malloc( (size_t) sizeof(VERTEX) * isostate.max_n_vertex[ilevel] );
  normal[ilevel] = (XYZ *) malloc( (size_t) sizeof(VERTEX) * isostate.max_n_vertex[ilevel] );

  triangl_status = (int *)       malloc(ntriangl   * sizeof(int));
  vertex_status  = (int *)       malloc(ntriangl*3 * sizeof(int));
  nver2triIN     = (int *)       malloc(ntriangl*3 * sizeof(int));
  /*sign           = (float *)    malloc(triangl    * sizeof(float));*/

  ver2triIN = (VERT2TRIG_INFO (*)[MAX_VERTEX_COORDINATION]) malloc((size_t )ntriangl*3 * sizeof(VERT2TRIG_INFO [MAX_VERTEX_COORDINATION]) );

  
  if ( !vertex[ilevel] || !normal[ilevel] || !triangl_status || !vertex_status || !nver2triIN )
    {
      xcError("allocation error in routine Marching Cubes !!!");
    }

  isostate.vertex_malloc[ilevel] = 1;

  if ( shade_model == ISOSURF_SHADE_SMOOTH ) 
    {      
      VertexRecognise(newgrd, ntriangl, triangles, triangl_status, 
		      tri2verIN[ilevel], vertex[ilevel], vertex_status, 
		      ver2triIN, nver2triIN, grid_info, &nvertex);
    } 
  else 
    {
      ivrt = 0; 
      nvertex = 3*ntriangl;

      for(i=0; i<ntriangl; i++)  
	{ 
	  vertex[ilevel][ivrt] = triangles[i].p[0]; tri2verIN[ilevel][i].ver_ind[0] = ivrt++; 
	  vertex[ilevel][ivrt] = triangles[i].p[1]; tri2verIN[ilevel][i].ver_ind[1] = ivrt++; 
	  vertex[ilevel][ivrt] = triangles[i].p[2]; tri2verIN[ilevel][i].ver_ind[2] = ivrt++; 	  
	}       
    }  
  
  SurfSmoothing(smooth_nstep, smooth_weight, vertex[ilevel], 
		vertex_status, ver2triIN, tri2verIN[ilevel], 
		nver2triIN, nvertex);

  if ( normals_model == ISOSURF_NORMALS_TRIANGLE ) 
    {
      SurfNmlAver(ntriangl, nvertex, normal[ilevel], vertex[ilevel], tri2verIN[ilevel], sign, 0);
    }
  else
    {
      /* TODO: for flat shading and gradient normals --> take normal from the triangle center */
      /* HERE */
      gradient_SurfNml(gridvertex, newgrd, spanning_vec, ntriangl, tri2verIN[ilevel], nvertex, vertex[ilevel], normal[ilevel]);
    }
  

  surf.numV[ilevel] = nvertex;

  /* assing *_P variables the right memmory locations */
  triangl_status_P  = triangl_status;
  vertex_status_P   = vertex_status;
  nver2triIN_P      = nver2triIN;
  sign_P            = sign;
  ver2triIN_P       = ver2triIN;

  /* t.k: there is a mass concerning the orientation of normals; try to remedy it */
  
  if ( normals_model == ISOSURF_NORMALS_TRIANGLE && isodata.cell_orientation == XC_LEFT ) {
    if ( isolevel < 0.0 ) {
      /* t.k: it seems that we must revert normals if isolevel is negative */
      for (i=0; i<surf.numV[ilevel]; i++) 
        RevertVectorfv( (float *) &(normal[ilevel][i]) );          
    }
  } else {
    if ( isolevel > 0.0 ) {
      /* t.k: it seems that we must revert normals if isolevel is positive */
      for (i=0; i<surf.numV[ilevel]; i++) 
        RevertVectorfv( (float *) &(normal[ilevel][i]) );          
    }
  }
}

/*

   ************************************************************************
   NOTE:

   This is a modified routine of the original one due to Paul Bourke
   (used with permission). It is enhanced for the present purpose. For
   the original routine ee "Modelling" on:
   http://astronomy.swin.edu.au/~pbourke/modelling/
   ************************************************************************

   Given a grid cell and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the array "triangles"
   will be loaded up with the vertices at most 5 triangular facets.
        0 will be returned if the grid cell is either totally above
   of totally below the isolevel.
*/
int Polygonise(GRIDCELL grid,
	       float isolevel,
	       GRID   ijk,
	       TRIANGLE *triangles,
	       TRIG_INFO *trig_info)
{
  int i,in,ntriang;
  int cubeindex;
  XYZ vertlist[12], ver0, ver1, ver2;

  int edgeTable[256]={
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };
  
  int triTable[256][16] =
{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
 {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
 {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
 {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
 {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
 {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
 {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
 {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
 {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
 {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
 {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
 {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
 {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
 {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
 {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
 {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
 {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
 {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
 {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
 {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
 {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
 {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
 {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
 {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
 {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
 {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
 {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
 {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
 {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
 {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
 {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
 {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
 {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
 {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
 {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
 {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
 {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
 {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
 {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
 {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
 {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
 {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
 {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
 {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
 {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
 {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
 {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
 {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
 {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
 {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
 {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
 {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
 {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
 {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
 {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
 {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
 {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
 {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
 {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
 {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
 {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
 {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
 {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
 {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
 {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
 {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
 {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
 {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
 {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
 {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
 {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
 {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
 {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
 {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
 {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
 {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
 {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
 {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
 {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
 {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
 {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
 {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
 {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
 {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
 {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
 {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
 {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
 {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
 {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
 {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
 {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
 {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
 {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
 {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
 {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
 {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
 {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
 {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
 {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
 {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
 {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
 {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
 {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
 {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
 {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
 {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
 {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
 {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
 {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
 {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
 {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
 {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
 {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
 {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
 {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
 {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
 {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
 {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
 {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
 {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
 {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
 {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
 {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
 {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
 {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
 {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
 {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
 {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
 {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
 {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
 {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
 {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
 {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
 {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
 {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
 {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
 {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
 {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
 {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
 {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
 {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
 {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
 {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
 {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
 {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
 {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
 {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
 {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
 {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
 {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
 {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
 {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
 {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
 {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
 {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
 {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
 {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
 {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
 {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
 {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
 {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
 {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
 {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
 {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
 {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
 {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
 {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
 {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
 {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
 {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
 {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
 {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
 {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
 {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
 {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
 {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
 {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
 {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
 {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
 {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
 {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

  /*
    Determine the index into the edge table which
    tells us which vertices are inside of the surface
  */
  cubeindex = 0;
  if (grid.val[0] < isolevel) cubeindex |= 1;
  if (grid.val[1] < isolevel) cubeindex |= 2;
  if (grid.val[2] < isolevel) cubeindex |= 4;
  if (grid.val[3] < isolevel) cubeindex |= 8;
  if (grid.val[4] < isolevel) cubeindex |= 16;
  if (grid.val[5] < isolevel) cubeindex |= 32;
  if (grid.val[6] < isolevel) cubeindex |= 64;
  if (grid.val[7] < isolevel) cubeindex |= 128;

  /* Cube is entirely in/out of the surface */
  if (edgeTable[cubeindex] == 0)
    return(0);

  /* Find the vertices where the surface intersects the cube */
  if (edgeTable[cubeindex] & 1)
    vertlist[0] =
      VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
  if (edgeTable[cubeindex] & 2)
    vertlist[1] =
      VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
  if (edgeTable[cubeindex] & 4)
    vertlist[2] =
      VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
  if (edgeTable[cubeindex] & 8)
    vertlist[3] =
      VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
  if (edgeTable[cubeindex] & 16)
    vertlist[4] =
      VertexInterp(isolevel,grid.p[4],grid.p[5],grid.val[4],grid.val[5]);
  if (edgeTable[cubeindex] & 32)
    vertlist[5] =
      VertexInterp(isolevel,grid.p[5],grid.p[6],grid.val[5],grid.val[6]);
  if (edgeTable[cubeindex] & 64)
    vertlist[6] =
      VertexInterp(isolevel,grid.p[6],grid.p[7],grid.val[6],grid.val[7]);
  if (edgeTable[cubeindex] & 128)
    vertlist[7] =
      VertexInterp(isolevel,grid.p[7],grid.p[4],grid.val[7],grid.val[4]);
  if (edgeTable[cubeindex] & 256)
    vertlist[8] =
      VertexInterp(isolevel,grid.p[0],grid.p[4],grid.val[0],grid.val[4]);
  if (edgeTable[cubeindex] & 512)
    vertlist[9] =
      VertexInterp(isolevel,grid.p[1],grid.p[5],grid.val[1],grid.val[5]);
  if (edgeTable[cubeindex] & 1024)
    vertlist[10] =
      VertexInterp(isolevel,grid.p[2],grid.p[6],grid.val[2],grid.val[6]);
  if (edgeTable[cubeindex] & 2048)
    vertlist[11] =
      VertexInterp(isolevel,grid.p[3],grid.p[7],grid.val[3],grid.val[7]);
  
  /* Create the triangle */
  ntriang = 0;
  for (i=0;triTable[cubeindex][i]!=-1;i+=3) {
    triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i  ]];
    triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i+1]];
    triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i+2]];

    ASSIGN_GRID2IJK( trig_info[ntriang], ijk);
    trig_info[ntriang].ver_status[0] = XC_TRUE;
    trig_info[ntriang].ver_status[1] = XC_TRUE;
    trig_info[ntriang].ver_status[2] = XC_TRUE;

    /* check for bad triangles !!! */
    in=0;
    ver0=triangles[ntriang].p[0];
    ver1=triangles[ntriang].p[1];
    ver2=triangles[ntriang].p[2];
    if (EQUAL_VERTICES(ver0,ver1)) in++;
    if (EQUAL_VERTICES(ver0,ver2)) in++;
    if (EQUAL_VERTICES(ver1,ver2)) in++;
    /* if triangle is good, count it */
    if(!in) ntriang++;
  }
  
   return(ntriang);
}

/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/

XYZ VertexInterp(float isolevel, XYZ p1, XYZ p2, float valp1, float valp2)
{
  float mu;
  XYZ p;
  
/*   if (ABS(valp1-valp2) < FLOAT_MINTOL) { */
/*     p.x = 0.5 * (p2.x + p1.x); */
/*     p.y = 0.5 * (p2.y + p1.y); */
/*     p.z = 0.5 * (p2.z + p1.z); */
/*     return(p); */
/*   } else if (ABS(isolevel-valp1) < FLOAT_MINTOL) { */
/*     return(p1); */
/*   } else if (ABS(isolevel-valp2) < FLOAT_MINTOL) { */
/*     return(p2); */
/*   } */

  if (ABS(valp1-valp2) < MINTOL_P) {
    /*fprintf(stderr, "polygonise_vertex_interp: valp1 == valp2 occurred !!!\n");*/
    p.x = 0.5 * (p2.x + p1.x);
    p.y = 0.5 * (p2.y + p1.y);
    p.z = 0.5 * (p2.z + p1.z);
    return(p);
  } else if (ABS(isolevel-valp1) < MINTOL_P) {
    /*fprintf(stderr, "polygonise_vertex_interp: valp1 == isolevel occurred !!!\n");*/
    return(p1);
  } else if (ABS(isolevel-valp2) < MINTOL_P) {
    /*fprintf(stderr, "polygonise_vertex_interp: valp2 == isolevel occurred !!!\n");*/
    return(p2);
  }

  mu = (isolevel - valp1) / (valp2 - valp1);
  
  /* testing */
  /*if (mu < 0.1)      mu = 0.0;
    else if (mu > 0.9) mu = 1.0;*/
  /* end */

  p.x = p1.x + mu * (p2.x - p1.x);
  p.y = p1.y + mu * (p2.y - p1.y);
  p.z = p1.z + mu * (p2.z - p1.z);
  
  return(p);
}





/*
   ************************************************************************
   NOTE:

   Below are modified routines of the original ones due to Paul Bourke
   (used with permission). It is enhanced for the present purpose. For
   the original routine ee "Modelling" on:
   http://astronomy.swin.edu.au/~pbourke/modelling/
   ************************************************************************

  Polygonise a tetrahedron given its vertices within a cube
  This is an alternative algorithm to polygonisegrid.
  It results in a smoother surface but more triangular facets.

         + 0
        /|\
       / | \
      /  |  \
  3  +-------+ 1
     \   |   /
      \  |  /
       \ | /
        \|/
         + 2

  It's main purpose is still to polygonise a gridded dataset and
  would normally be called 6 times, one for each tetrahedron making
  up the grid cell.
*/

int PolygoniseTetrahedral(GRIDCELL grid, 
			  float     iso, 
			  GRID      ijk, 
			  TRIANGLE  *triangl, 
			  TRIG_INFO *trig_info);

static 
int PolygoniseTri(GRIDCELL  g, 
		  float     iso, 
		  GRID      ijk, 
		  TRIANGLE  *triangl, 
		  TRIG_INFO *trig_info,
		  int v0, int v1, int v2, int v3);


int PolygoniseTetrahedral(GRIDCELL  grid, 
			  float     iso, 
			  GRID      ijk, 
			  TRIANGLE  *triangl, 
			  TRIG_INFO *trig_info)
{
  int ntriangl = 0;

  /* 6-tetrah. decomposition */
  ntriangl += PolygoniseTri(grid,iso,ijk,&triangl[ntriangl],&trig_info[ntriangl],0,2,3,7);
  ntriangl += PolygoniseTri(grid,iso,ijk,&triangl[ntriangl],&trig_info[ntriangl],0,7,6,2); /*0,7,6,2*/
  ntriangl += PolygoniseTri(grid,iso,ijk,&triangl[ntriangl],&trig_info[ntriangl],0,7,4,6); /*0,7,4,6*/
  ntriangl += PolygoniseTri(grid,iso,ijk,&triangl[ntriangl],&trig_info[ntriangl],0,6,1,2);
  ntriangl += PolygoniseTri(grid,iso,ijk,&triangl[ntriangl],&trig_info[ntriangl],0,4,1,6); /* 0,4,1,6*/
  ntriangl += PolygoniseTri(grid,iso,ijk,&triangl[ntriangl],&trig_info[ntriangl],5,6,1,4); /* reserve: 4.5.6.1*/

  /* 5-tetra decomposition: this doesn't work */
  /* ntriangl += PolygoniseTri(grid,iso,ijk,&triangl[ntriangl],&trig_info[ntriangl],0,5,1,2); */
  /* ntriangl += PolygoniseTri(grid,iso,ijk,&triangl[ntriangl],&trig_info[ntriangl],0,7,4,5);  */
  /* ntriangl += PolygoniseTri(grid,iso,ijk,&triangl[ntriangl],&trig_info[ntriangl],7,6,5,2);  */
  /* ntriangl += PolygoniseTri(grid,iso,ijk,&triangl[ntriangl],&trig_info[ntriangl],0,2,3,7); */
  /* ntriangl += PolygoniseTri(grid,iso,ijk,&triangl[ntriangl],&trig_info[ntriangl],0,5,2,7);  */
 
  return ntriangl;
}


static
int PolygoniseTri(GRIDCELL  g,
		  float     iso,
		  GRID      ijk, 
		  TRIANGLE  *triangl,
		  TRIG_INFO *trig_info,
		  int v0,int v1,int v2,int v3)
{
  int ntriangl, aux_ntriangl = 0;
  int i, in, triindex;
  TRIANGLE aux_triangl[2];

   /*
      Determine which of the 16 cases we have given which vertices
      are above or below the isosurface
   */
   triindex = 0;
   if (g.val[v0] > iso) triindex |= 1;
   if (g.val[v1] > iso) triindex |= 2;
   if (g.val[v2] > iso) triindex |= 4;
   if (g.val[v3] > iso) triindex |= 8;

   /* Form the vertices of the triangles for each case */

   switch (triindex) {

   case 0x00:
   case 0x0F:
      break;

   case 0x0E: 
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      aux_ntriangl++;
      break;
   case 0x01:
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      aux_ntriangl++;
      break;

   case 0x0D: 
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v1],g.p[v0],g.val[v1],g.val[v0]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      aux_ntriangl++;
      break;
   case 0x02:
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v1],g.p[v0],g.val[v1],g.val[v0]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      aux_ntriangl++;
      break;

   case 0x0C: 
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      aux_ntriangl++;
      aux_triangl[1].p[0] = aux_triangl[0].p[2];
      aux_triangl[1].p[1] = aux_triangl[0].p[1];
      aux_triangl[1].p[2] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      aux_ntriangl++;
      break;
   case 0x03:
     /* correcting */
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      aux_ntriangl++;
      aux_triangl[1].p[0] = aux_triangl[0].p[2];
      aux_triangl[1].p[1] = aux_triangl[0].p[1];
      aux_triangl[1].p[2] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      aux_ntriangl++;
      break;

   case 0x0B: 
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v2],g.p[v0],g.val[v2],g.val[v0]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v2],g.p[v1],g.val[v2],g.val[v1]);
      aux_ntriangl++;
      break;
   case 0x04:
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v2],g.p[v0],g.val[v2],g.val[v0]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v2],g.p[v1],g.val[v2],g.val[v1]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      aux_ntriangl++;
      break;

   case 0x0A: 
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      aux_ntriangl++;
      aux_triangl[1].p[0] = aux_triangl[0].p[0];
      aux_triangl[1].p[1] = aux_triangl[0].p[2];
      aux_triangl[1].p[2] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      aux_ntriangl++;
      break;
   case 0x05:
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v0],g.p[v3],g.val[v0],g.val[v3]);
      aux_ntriangl++;
      aux_triangl[1].p[0] = aux_triangl[0].p[0];
      aux_triangl[1].p[1] = VertexInterp(iso,g.p[v1],g.p[v2],g.val[v1],g.val[v2]);
      aux_triangl[1].p[2] = aux_triangl[0].p[1];
      aux_ntriangl++;
      break;

   case 0x09: 
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      aux_ntriangl++;
      aux_triangl[1].p[0] = aux_triangl[0].p[0];
      aux_triangl[1].p[1] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      aux_triangl[1].p[2] = aux_triangl[0].p[1];
      aux_ntriangl++;
      break;
   case 0x06:
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v0],g.p[v1],g.val[v0],g.val[v1]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v1],g.p[v3],g.val[v1],g.val[v3]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v2],g.p[v3],g.val[v2],g.val[v3]);
      aux_ntriangl++;
      aux_triangl[1].p[0] = aux_triangl[0].p[0];
      aux_triangl[1].p[1] = aux_triangl[0].p[2];
      aux_triangl[1].p[2] = VertexInterp(iso,g.p[v0],g.p[v2],g.val[v0],g.val[v2]);
      aux_ntriangl++;
      break;

   case 0x07: 
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v3],g.p[v0],g.val[v3],g.val[v0]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v3],g.p[v1],g.val[v3],g.val[v1]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v3],g.p[v2],g.val[v3],g.val[v2]);
      aux_ntriangl++;
      break;
   case 0x08:
      aux_triangl[0].p[0] = VertexInterp(iso,g.p[v3],g.p[v0],g.val[v3],g.val[v0]);
      aux_triangl[0].p[1] = VertexInterp(iso,g.p[v3],g.p[v2],g.val[v3],g.val[v2]);
      aux_triangl[0].p[2] = VertexInterp(iso,g.p[v3],g.p[v1],g.val[v3],g.val[v1]);
      aux_ntriangl++;
      break;
   }

   /* eliminate bad triangles, i.e. those that contain identical points */

   ntriangl = 0;
   for (i=0; i<aux_ntriangl; i++)
     {
       in=0;
       if (EQUAL_VERTICES(aux_triangl[i].p[0], aux_triangl[i].p[1])) in++;
       if (EQUAL_VERTICES(aux_triangl[i].p[0], aux_triangl[i].p[2])) in++;
       if (EQUAL_VERTICES(aux_triangl[i].p[1], aux_triangl[i].p[2])) in++;
       
      if(!in) {
	/* triangle is good, count it */
	
	triangl[ntriangl].p[0] = aux_triangl[i].p[0];
	triangl[ntriangl].p[1] = aux_triangl[i].p[1];
	triangl[ntriangl].p[2] = aux_triangl[i].p[2];
	
	ASSIGN_GRID2IJK( trig_info[ntriangl], ijk);
	trig_info[ntriangl].ver_status[0] = XC_TRUE;
	trig_info[ntriangl].ver_status[1] = XC_TRUE;
	trig_info[ntriangl].ver_status[2] = XC_TRUE;

	ntriangl++;
      }
    }
      
   return(ntriangl);
}

