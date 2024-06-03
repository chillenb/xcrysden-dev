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
 * Source: $XCRYSDEN_TOPDIR/C/isoline.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#define XC_CPP_NO_STDIO
#include "struct.h"
#include "isosurf.h"
#include "memory.h"

extern void GetIsoLine2D_Attributes(int obj, int cb, int fn);
extern XYZ VertexInterp(float isolevel, XYZ p1, XYZ p2, 
			float valp1, float valp2);

extern PLANEVERTEX ***plvertex;

typedef struct {
  XYZ   p[3];
  float val[3];
} TRIANGL_GRID;

int IsoLine2D(int obj, int cb, int fn, int nx, int ny);
static void Segmentise( int obj, TRIANGL_GRID triangl, int nlevel, 
			float level[ISOLINE_MAXLEVEL] );

extern ISOSTATE isostate;
int
IsoLine2D(int obj, int cb, int fn, int nx, int ny)
{
  register int i, j;
  TRIANGL_GRID triangl1;
  TRIANGL_GRID triangl2;
  int nnx, nny;

  GetIsoLine2D_Attributes(obj, cb, fn);
  /*  reset the seqment counter */
  for (i=0; i<ISOLINE_MAXLEVEL; i++)
    isoline2D[obj].iseg[i] = 0;

  nnx = nx - 1;
  nny = ny - 1;
  for (i=0; i<nnx; i++)
    for (j=0; j<nny; j++) 
      {
	/* first triangle */
	triangl1.p[0].x   = plvertex[obj][i][j].p[0];
	triangl1.p[0].y   = plvertex[obj][i][j].p[1];
	triangl1.p[0].z   = plvertex[obj][i][j].p[2];
	triangl1.val[0]   = plvertex[obj][i][j].val;

	triangl1.p[1].x   = plvertex[obj][i][j+1].p[0];
	triangl1.p[1].y   = plvertex[obj][i][j+1].p[1];
	triangl1.p[1].z   = plvertex[obj][i][j+1].p[2];
	triangl1.val[1]   = plvertex[obj][i][j+1].val;

	triangl1.p[2].x   = plvertex[obj][i+1][j].p[0];
	triangl1.p[2].y   = plvertex[obj][i+1][j].p[1];
	triangl1.p[2].z   = plvertex[obj][i+1][j].p[2];
	triangl1.val[2]   = plvertex[obj][i+1][j].val;

	Segmentise( obj, triangl1, isoline2D[obj].nlevel, 
		    isoline2D[obj].level );
 
	/* second triangle */
	triangl2.p[0].x   = plvertex[obj][i][j+1].p[0];
	triangl2.p[0].y   = plvertex[obj][i][j+1].p[1];
	triangl2.p[0].z   = plvertex[obj][i][j+1].p[2];
        triangl2.val[0]   = plvertex[obj][i][j+1].val;

        triangl2.p[1].x   = plvertex[obj][i+1][j+1].p[0];
        triangl2.p[1].y   = plvertex[obj][i+1][j+1].p[1];
        triangl2.p[1].z   = plvertex[obj][i+1][j+1].p[2];
        triangl2.val[1]   = plvertex[obj][i+1][j+1].val;

        triangl2.p[2].x   = plvertex[obj][i+1][j].p[0];
        triangl2.p[2].y   = plvertex[obj][i+1][j].p[1];
        triangl2.p[2].z   = plvertex[obj][i+1][j].p[2];
        triangl2.val[2]   = plvertex[obj][i+1][j].val;

        Segmentise( obj, triangl2, isoline2D[obj].nlevel, 
		    isoline2D[obj].level );
      }
  return XC_OK;
}


static void
Segmentise( int obj, TRIANGL_GRID triangl, int nlevel, 
	    float level[ISOLINE_MAXLEVEL] )
{
  int   il, nseg;
  float v[3];
  XYZ   p[3];
  
  for (il=0; il<nlevel; il++) {
    /* if the malloc space get to small for isoline2D[obj], 
     * malloc bigger vector */
    if ( isoline2D[obj].iseg[il] + 2 > isostate.max_n_isoline2D[obj][il] )
      isoline2D[obj].segment[il] = 
	xcRealloc2xBiggerLINE( &(isoline2D[obj].segment[il][0]), 
			       &isostate.max_n_isoline2D[obj][il] );

    nseg=0;
    v[0] = triangl.val[0] - level[il];
    v[1] = triangl.val[1] - level[il];
    v[2] = triangl.val[2] - level[il];
    
    if ( v[0]*v[1] < 0.0 ) {
      p[nseg] = VertexInterp(level[il], 
			     triangl.p[0], triangl.p[1], 
			     triangl.val[0], triangl.val[1]);
      nseg++;
    } 
    if ( v[1]*v[2] < 0.0 ) {
      p[nseg] = VertexInterp(level[il],
                             triangl.p[1], triangl.p[2], 
			     triangl.val[1], triangl.val[2]);
      nseg++;
    }
    if ( v[2]*v[0] < 0.0 ) {
      p[nseg] = VertexInterp(level[il],
                             triangl.p[2], triangl.p[0], 
			     triangl.val[2], triangl.val[0]);
      nseg++;
    }

    if (nseg == 2) {
      isoline2D[obj].segment[il][isoline2D[obj].iseg[il]].p1 = p[0];
      isoline2D[obj].segment[il][isoline2D[obj].iseg[il]].p2 = p[1];
      isoline2D[obj].iseg[il]++;
    }
  }
}
