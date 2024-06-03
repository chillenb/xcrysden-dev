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
 * Source: $XCRYSDEN_TOPDIR/C/MarchCubes.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <stdio.h> 
#include <math.h>
#include <stdlib.h> 
#include "struct.h"
#include "isosurf.h"
#include "memory.h"
#include "xcfunc.h"

#define ASSIGN_GRID2IJK(a,b) a.i = b.ix; a.j = b.iy; a.k = b.iz
#define EQUAL_VERTICES(a,b)    ( fabs(a.x-b.x) < MINTOL_F && fabs(a.y-b.y) < MINTOL_F && fabs(a.z-b.z) < MINTOL_F )
#define GET_INDEX(i,j,k,N,M,L) ( (k + L*(j + M*i)) )

extern int Polygonise(GRIDCELL grid,
		      float isolevel,
		      GRID   ijk,
		      TRIANGLE *triangles,
		      TRIG_INFO *trig_info);
extern void VertexRecognise(GRID          grd,
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
extern void SurfSmoothing(int            nstep, 
			  float          weight,
			  XYZ            *vertex,
			  int            *vertex_status,
			  VERT2TRIG_INFO (*ver2triIN)[MAX_VERTEX_COORDINATION],
			  TRIG_INFO      *tri2verIN,
			  int            *nver2triIN,
			  int            nvertex);
extern void SurfNmlAver(int       ntriangl,
			int       nvertex,
			XYZ       *nml,
			XYZ       *vertex,	    
			TRIG_INFO *tri2verIN,
			float     *sign,
			int       is_tetrahedral);

TRIANGLE  *xcRealloc2xBiggerTriangl(TRIANGLE *vec, int *size);
TRIG_INFO *xcRealloc2xBiggerTrigInfo(TRIG_INFO *vec, int *size);


void
MarchCubeNew(ISOSURFACE *iso,
	     float ***gridvalue,
	     float vec[][3],
	     float origin[],
	     int ix, int iy, int iz,
	     float isolevel,
	     int   algorithm,
	     int   shade_model,
	     int   normals_model)
	      
{
  register int i, l, iix, iiy, iiz;
  register float di, dj, dk, di0, dj0, dk0, di1, dj1, dk1;
  float incr_vec[3][3];
  int mem_size, ntriangl, ntr, nvrt, index, new_maxelem;
  GRID           gr;
  GRIDCELL       grid;
  GRID_INFO      *grid_info;
  TRIANGLE       *triangles;
  int            nvertex;
  int            *vertex_status;
  float          *sign=NULL;

  /* DEBUGGING_PURPOSES_BEGIN */
/*   if (ABS(ix) > 100) ix=8; */
/*   if (ABS(iy) > 100) iy=12; */
/*   if (ABS(iz) > 100) iz=12; */
  /* END */

  gr.nx=ix;
  gr.ny=iy;
  gr.nz=iz;

  iix = ix - 1;
  iiy = iy - 1;
  iiz = iz - 1;
  for (i=0; i<3; i++) {
    incr_vec[0][i] = vec[0][i]/(float)iix;
    incr_vec[1][i] = vec[1][i]/(float)iiy;
    incr_vec[2][i] = vec[2][i]/(float)iiz;
  }

  /*
    if (iso->vertex)     free( (FREE_ARG) iso->vertex );
    if (iso->normal)     free( (FREE_ARG) iso->normal );
    if (iso->nver2triIN) free( (FREE_ARG) iso->nver2triIN );
    if (iso->tri2verIN)  free( (FREE_ARG) iso->tri2verIN);
  */
  mem_size = MAX(gr.nx, gr.ny) * MAX(gr.ny, gr.nz);

  triangles = (TRIANGLE *) malloc( (size_t) sizeof(TRIANGLE) * mem_size );
  if ( !triangles ) 
    {
      xcError("allocation error for triangles");
    }

  iso->tri2verIN = (TRIG_INFO *) realloc((void *) iso->tri2verIN, (size_t)mem_size * sizeof(TRIG_INFO));
  grid_info      = (GRID_INFO *) malloc((size_t) (gr.nx*gr.ny*gr.nz) * sizeof(GRID_INFO));
  if ( !iso->tri2verIN || !grid_info) 
    {
      xcError("allocation error in Marching Cubes");
    }

  ntriangl=0;
  for (gr.ix = 0; gr.ix < iix; gr.ix++) 
    {
      di = (float) gr.ix;

      for (gr.iy = 0; gr.iy < iiy; gr.iy++) 
	{
	  dj = (float) gr.iy;

	  for (gr.iz = 0; gr.iz < iiz; gr.iz++) 
	    {
	      dk = (float) gr.iz;
	      
	      for (l=0; l<3; l++) 
		{
		  di0 =       di*incr_vec[0][l];
		  di1 = (di+1.0)*incr_vec[0][l];
		  
		  dj0 =       dj*incr_vec[1][l];
		  dj1 = (dj+1.0)*incr_vec[1][l];
		  
		  dk0 =       dk*incr_vec[2][l];
		  dk1 = (dk+1.0)*incr_vec[2][l];
		  
		  if(l==0) {
		    grid.p[0].x = origin[l] + di0 + dj0 + dk0;
		    grid.p[1].x = origin[l] + di0 + dj1 + dk0;
		    grid.p[2].x = origin[l] + di1 + dj1 + dk0;
		    grid.p[3].x = origin[l] + di1 + dj0 + dk0;
		    grid.p[4].x = origin[l] + di0 + dj0 + dk1;
		    grid.p[5].x = origin[l] + di0 + dj1 + dk1;
		    grid.p[6].x = origin[l] + di1 + dj1 + dk1;
		    grid.p[7].x = origin[l] + di1 + dj0 + dk1;
		  } else if (l==1) {
		    grid.p[0].y	= origin[l] + di0 + dj0 + dk0;
		    grid.p[1].y	= origin[l] + di0 + dj1 + dk0;
		    grid.p[2].y	= origin[l] + di1 + dj1 + dk0;
		    grid.p[3].y	= origin[l] + di1 + dj0 + dk0;
		    grid.p[4].y	= origin[l] + di0 + dj0 + dk1;
		    grid.p[5].y	= origin[l] + di0 + dj1 + dk1;
		    grid.p[6].y	= origin[l] + di1 + dj1 + dk1;
		    grid.p[7].y	= origin[l] + di1 + dj0 + dk1;
		  } else {
		    grid.p[0].z = origin[l] + di0 + dj0 + dk0;
		    grid.p[1].z = origin[l] + di0 + dj1 + dk0;
		    grid.p[2].z = origin[l] + di1 + dj1 + dk0;
		    grid.p[3].z = origin[l] + di1 + dj0 + dk0;
		    grid.p[4].z = origin[l] + di0 + dj0 + dk1;
		    grid.p[5].z = origin[l] + di0 + dj1 + dk1;
		    grid.p[6].z = origin[l] + di1 + dj1 + dk1;
		    grid.p[7].z = origin[l] + di1 + dj0 + dk1;
		  }
		}

	      new_maxelem = ( algorithm == ISOSURF_MARCHING_CUBES ? 5 : 12 );
	      
	      /* take care of reallocation if necessary !!! */
	      if ( ntriangl + new_maxelem >= mem_size ) 
		{
		  int mem_size1 = mem_size;
		  triangles      = xcRealloc2xBiggerTriangl ( triangles, &mem_size );
		  iso->tri2verIN = xcRealloc2xBiggerTrigInfo( iso->tri2verIN, &mem_size1 );
		}
	      
	      grid.val[0] = gridvalue[gr.ix  ][gr.iy  ][gr.iz  ]; 
	      grid.val[1] = gridvalue[gr.ix  ][gr.iy+1][gr.iz  ];
	      grid.val[2] = gridvalue[gr.ix+1][gr.iy+1][gr.iz  ];
	      grid.val[3] = gridvalue[gr.ix+1][gr.iy  ][gr.iz  ];
	      grid.val[4] = gridvalue[gr.ix  ][gr.iy  ][gr.iz+1];
	      grid.val[5] = gridvalue[gr.ix  ][gr.iy+1][gr.iz+1];
	      grid.val[6] = gridvalue[gr.ix+1][gr.iy+1][gr.iz+1];
	      grid.val[7] = gridvalue[gr.ix+1][gr.iy  ][gr.iz+1];

	      /* index=GET_INDEX(gr.ix,gr.iy,gr.iz, gr.nx, gr.ny, gr.nz); */
	      index=GET_INDEX(gr.ix,gr.iy,gr.iz, iix, iiy, iiz);	
	      grid_info[index].t0=ntriangl;

	      if ( algorithm == ISOSURF_MARCHING_CUBES ) 
		{
		  /* ===> POLYGONISE: MarchingCubes <=== */
		  ntriangl += Polygonise(grid, isolevel, gr, &triangles[ntriangl], 
					 (TRIG_INFO *) &iso->tri2verIN[ntriangl]);
		}
	      else 
		{
		  /* ===> POLYGONISE: TetrahedralDecomposition <=== */
		  ntriangl += PolygoniseTetrahedral(grid, isolevel, gr,  &triangles[ntriangl], 
						    (TRIG_INFO *) &iso->tri2verIN[ntriangl]);
		}

	      grid_info[index].nt = ntriangl - grid_info[index].t0 + 1;
	    } /* iz */
	}
    }
  
  nvertex = ntriangl * 3;
  
  /* just in case the ntriangl is 0, then we would have allocation error */
  ntr=ntriangl;
  nvrt=nvertex;
  if (ntriangl == 0) 
    {
      ntr=1;
      nvrt=1;
    }

  vertex_status = (int *) malloc(ntr*3 * sizeof(int));

  iso->triangl_status = (int *) realloc((void*) iso->triangl_status, ntr * sizeof(int));
  iso->vertex         = (XYZ *) realloc((void*) iso->vertex, (size_t) nvrt * sizeof(VERTEX) );
  iso->vertex_orig    = (XYZ *) realloc((void*) iso->vertex_orig, (size_t) nvrt * sizeof(VERTEX) );
  iso->normal         = (XYZ *) realloc((void*) iso->normal, (size_t) nvrt * sizeof(VERTEX) );
  iso->nver2triIN     = (int *) realloc((void*) iso->nver2triIN, ntr*3 * sizeof(int));
  iso->ver2triIN      = (VERT2TRIG_INFO (*)[MAX_VERTEX_COORDINATION]) realloc((void*) iso->ver2triIN, (size_t ) ntr*3 * sizeof(VERT2TRIG_INFO [MAX_VERTEX_COORDINATION]) );
  
  if ( !iso->triangl_status || !vertex_status || !iso->vertex || !iso->vertex_orig
       || !iso->normal || !iso->nver2triIN || !iso->ver2triIN )
    {
      xcError("allocation error in routine Marching Cubes !!!");
    }

  VertexRecognise(gr, ntriangl, triangles, iso->triangl_status, 
		  iso->tri2verIN, iso->vertex, vertex_status, 
		  iso->ver2triIN, iso->nver2triIN, grid_info, &nvertex);

  /* if ( shade_model == ISOSURF_SHADE_SMOOTH )  */
  /*   {             */
  /*     VertexRecognise(gr, ntriangl, triangles, iso->triangl_status,  */
  /* 		      iso->tri2verIN, iso->vertex, vertex_status,  */
  /* 		      iso->ver2triIN, iso->nver2triIN, grid_info, &nvertex); */
  /*   } */
  /* else */
  /*   { */
  /*     int ivrt = 0;  */
  /*     nvertex = 3*ntriangl; */
      
  /*     for(i=0; i<ntriangl; i++)   */
  /* 	{  */
  /* 	  iso->vertex[ivrt] = triangles[i].p[0]; iso->tri2verIN[i].ver_ind[0] = ivrt++;  */
  /* 	  iso->vertex[ivrt] = triangles[i].p[1]; iso->tri2verIN[i].ver_ind[1] = ivrt++;  */
  /* 	  iso->vertex[ivrt] = triangles[i].p[2]; iso->tri2verIN[i].ver_ind[2] = ivrt++; 	   */
  /* 	}        */
  /*   }   */

  /* copy vertices to vertex_orig !!! */

  for(i=0; i<nvertex; i++)
    {
      iso->vertex_orig[i] = iso->vertex[i];
    }

  /*for(i=0; i<ntriangl; i++)
    iso->triangl_status[i] = 1;*/
  
  SurfSmoothing(iso->smooth_nstep, iso->smooth_weight, iso->vertex, 
		vertex_status, iso->ver2triIN, iso->tri2verIN, 
		iso->nver2triIN, nvertex);
  
  if ( normals_model == ISOSURF_NORMALS_TRIANGLE ) 
    {
      SurfNmlAver(ntriangl, nvertex, iso->normal, iso->vertex, iso->tri2verIN, sign, 0);
    }
  else
    {
      /* TODO: for flat shading and gradient normals --> take normal from the triangle center */
      /* HERE */
      gradient_SurfNmlNew(origin, incr_vec, gridvalue, gr, ntriangl, iso->tri2verIN, nvertex, iso->vertex, iso->normal);
    }
  
  /* t.k: it seems that we must revert normals if isolevel is positive */
  /*
  for (i=0; i<nvertex; i++) {
    if ( isolevel > 0.0 ) 
      RevertVectorfv( (float *) &(iso->normal[i]) );
  }
  */

  iso->nvertex       = nvertex;
  iso->ntriangl      = ntriangl;
  iso->revertnormals = 0;
  
  free( (FREE_ARG) triangles );
  free( (FREE_ARG) grid_info );
  free( (FREE_ARG) vertex_status );
}

void
xcSurfSmoothing( ISOSURFACE *iso, float isolevel)
{
  register int i;
  int *vertex_status;
  float *sign = NULL;

  if (iso->ntriangl == 0) return;
  if (iso->smooth_nstep == 0) return;
  
  vertex_status = (int *) xcCalloc(3*iso->ntriangl, sizeof(int)); 
  
  for(i=0; i<iso->nvertex; i++)
    iso->vertex[i] = iso->vertex_orig[i];

  SurfSmoothing(iso->smooth_nstep, iso->smooth_weight, iso->vertex, 
		vertex_status, iso->ver2triIN, iso->tri2verIN, 
		iso->nver2triIN, iso->nvertex);

  if (1)
    {
      /* normals_algorithm == ISOSURF_NORMALS_TRIANGLE */
      SurfNmlAver(iso->ntriangl, iso->nvertex, iso->normal, iso->vertex, iso->tri2verIN, sign, 0);
    } 
  else 
    {
      /* enter code here */
    }
  
  xcFree((void*) vertex_status);
}
