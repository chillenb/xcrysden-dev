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
 * Source: $XCRYSDEN_TOPDIR/C/isoInterpolate.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <stdio.h>
#include "struct.h"
#include "isosurf.h"
#include "memory.h"

extern void xcRegularSplineInt2(float gridX[], float gridY[], float **Fn,
				int nX, int nY, float ndegree, float **newFn);
extern void xcRegularSplineInt3(float gridX[], float gridY[], float gridZ[], 
				float ***Fn, int nX, int nY, int nZ,
				float ndegree, float ***newFn);

extern ISOSTATE isostate;
extern PLANEVERTEX ***plvertex;

void
isoInterpolate( int degree )
{  
  register int i, j, k;

  if ( degree == 0 ) { /* I can not imagine what degree 0 would stand for */
    degree = 1;
  }  

  if ( degree < 0 ) { /* so far lowering of resolution is not supported */
    degree = 1;
  }

  if (isodata.dim[ISOOBJ_BASE] == 3 ) {
    /*
     * if grd == newgrd && degree == 1 --> we already have what we want
     */
    if ( grd.nx == newgrd.nx && grd.ny == newgrd.ny && 
	 grd.nz == newgrd.nz && degree == 1) {
      return;
    } else {
      float ***newGv;
      float ***gv = xcMallocTensor3f( grd.nx, grd.ny, grd.nz );
      
      /*
       * forget about gridvertex and get new one
       */
      if ( isostate.gridvertex_malloc ) {
	xcFree_GRIDVERTEX( gridvertex );
	isostate.gridvertex_malloc = 0;
      }

      /*
       * rewind bin_vertex_fp && read the data
       */
      fflush( isodata.bin_vertex_fp );
      rewind( isodata.bin_vertex_fp );
      for(i=0; i<grd.nz; i++)
	for(j=0; j<grd.ny; j++)
	  for(k=0; k<grd.nx; k++)
	    fread(&gv[k][j][i], sizeof(float), 1, isodata.bin_vertex_fp);
             
      /*
       * malloc new gridvertex
       */
      if ( degree < 1 ) {
	/*
	 * we want to lower the resolution of the grid
	 */

	/* NOT SUPPORTED YET */
	newGv = NULL;
	xcError("interpolation degree < 1 not supported yet !!!");
      } else if ( degree == 1 ) {
	/*
	 * just copy gv -> newGv
	 */
	newgrd.nx = grd.nx;
	newgrd.ny = grd.ny;
	newgrd.nz = grd.nz;
	newGv = gv;
      } else {
	/*
	 * perform the interpolation
	 */
	int max;
	float *x;
	
	/* assing x[] */
	max = (grd.nx > grd.ny) ? grd.nx : grd.ny;
	max = (max > grd.nz) ? max : grd.nz;
	x = xcMallocVectorf( max );
	for(i=0; i<max; i++)
	  x[i] = (float) i;
	
	newgrd.nx = degree * (grd.nx-1) + 1;
	newgrd.ny = degree * (grd.ny-1) + 1;
	newgrd.nz = degree * (grd.nz-1) + 1;
	newGv = xcMallocTensor3f( newgrd.nx, newgrd.ny, newgrd.nz );
	xcRegularSplineInt3(x, x, x, gv, grd.nx, grd.ny, grd.nz, 
			    degree, newGv);
	xcFree_Vectorf( x );
	xcFree_Tensor3f( gv );
      }

      /* malloc new gridvertex */
      gridvertex = xcMallocGRIDVERTEX( newgrd.nx, newgrd.ny, newgrd.nz );
      isostate.gridvertex_malloc = 1;
      /*
       * now copy newGv to gridvertex[][][].val && assing gridvertex[][][].p
       */
      for (i=0; i<newgrd.nx; i++)
	for (j=0; j<newgrd.ny; j++)
	  for (k=0; k<newgrd.nz; k++) {
      	    gridvertex[i][j][k].val = newGv[i][j][k];
	    gridvertex[i][j][k].p.x =
	      ((float) i / (newgrd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].x +
	      ((float) j / (newgrd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].x +
	      ((float) k / (newgrd.nz - 1)) * isodata.vec[ISOOBJ_BASE][2].x;  
	    gridvertex[i][j][k].p.y =   
	      ((float) i / (newgrd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].y + 
	      ((float) j / (newgrd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].y +
	      ((float) k / (newgrd.nz - 1)) * isodata.vec[ISOOBJ_BASE][2].y; 
	    gridvertex[i][j][k].p.z =   
	      ((float) i / (newgrd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].z +  
	      ((float) j / (newgrd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].z +  
	      ((float) k / (newgrd.nz - 1)) * isodata.vec[ISOOBJ_BASE][2].z;   
	  }
      xcFree_Tensor3f( newGv );
    }
  } /* isodata.dim[obj] */
  else if ( isodata.dim[ISOOBJ_BASE] == 2 ) {
    /*
     * if grd == newgrd && degree == 1 --> we already have what we want
     */
    if ( grd.nx == newgrd.nx && grd.ny == newgrd.ny && degree == 1) {
      return;
    } else {
      float **newGv;
      float **gv = xcMallocMatrixf( grd.nx, grd.ny );
      
      /*
       * forget about plvertex[ISOOBJ_BASE] and get new one
       */
      if ( isostate.plvertex_malloc[ISOOBJ_BASE] ) {
	xcFree_PLANEVERTEX( plvertex[ISOOBJ_BASE] );
	isostate.plvertex_malloc[ISOOBJ_BASE] = 0;
      }

      /*
       * rewind bin_vertex_fp && read the data
       */
      fflush( isodata.bin_vertex_fp );
      rewind( isodata.bin_vertex_fp );
      for(j=0; j<grd.ny; j++)
	for(k=0; k<grd.nx; k++)
	  fread(&gv[k][j], sizeof(float), 1, isodata.bin_vertex_fp);
             
      /*
       * malloc new gridvertex
       */
      if ( degree < 1 ) {
	/*
	 * we want to lower the resolution of the grid
	 */

	/* NOT SUPPORTED YET */
	newGv = NULL;
	xcError("interpolation degree < 1 not supported yet !!!");
      } else if ( degree == 1 ) {
	/*
	 * just copy gv -> newGv
	 */
	newgrd.nx = grd.nx;
	newgrd.ny = grd.ny;
	newGv = gv;
      } else {
	/*
	 * perform the interpolation
	 */
	int max;
	float *x;
	
	/* assing x[] */
	max = (grd.nx > grd.ny) ? grd.nx : grd.ny;
	x = xcMallocVectorf( max );
	for(i=0; i<max; i++)
	  x[i] = (float) i;
	
	newgrd.nx = degree * (grd.nx-1) + 1;
	newgrd.ny = degree * (grd.ny-1) + 1;
	newGv = xcMallocMatrixf( newgrd.nx, newgrd.ny );
	xcRegularSplineInt2(x, x, gv, grd.nx, grd.ny, degree, newGv);
	xcFree_Vectorf( x );
	xcFree_Matrixf( gv );
      }

      /* malloc new gridvertex */
      plvertex[ISOOBJ_BASE] = xcMallocPLANEVERTEX( newgrd.nx, newgrd.ny );
      isostate.plvertex_malloc[ISOOBJ_BASE] = 1;
      /*
       * now copy newGv to gridvertex[][][].val && assing gridvertex[][][].p
       */
      for (i=0; i<newgrd.nx; i++)
	for (j=0; j<newgrd.ny; j++) {
	  plvertex[ISOOBJ_BASE][i][j].val = newGv[i][j];
	  plvertex[ISOOBJ_BASE][i][j].p[0] =
	    ((float) i / (newgrd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].x +
	    ((float) j / (newgrd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].x;
	  plvertex[ISOOBJ_BASE][i][j].p[1] =   
	    ((float) i / (newgrd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].y + 
	    ((float) j / (newgrd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].y;
	  plvertex[ISOOBJ_BASE][i][j].p[2] =   
	    ((float) i / (newgrd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].z +  
	    ((float) j / (newgrd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].z;  
	}
      xcFree_Matrixf( newGv );
    }
  }
}
	 




