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
 * Source: $XCRYSDEN_TOPDIR/C/cryClip.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "isosurf.h"
#include "memory.h"
#include "xcfunc.h"

/* cryClip.c */
int ClipSimple(XYZ *p, XYZ n, XYZ p0);
int ClipSimple3(XYZ *p, float n[3], float p0[3]);
ISOSURFACE *cryIsosurfClip(ISOSURFACE *in, XYZ n, XYZ p0);

/* 
   Simple Cliping; either IN or OUT
   triangle is IN if at least one vertex is "inside the plane"

   The normal to the plane is n
   A point on the plane is p0

*/
int ClipSimple(XYZ *p, XYZ n, XYZ p0)
{
   float A,B,C,D;
   float l;
   float side[3];

   /*
      Determine the equation of the plane as
      Ax + By + Cz + D = 0
   */
   l = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);
   A = n.x / l;
   B = n.y / l;
   C = n.z / l;
   D = -(n.x*p0.x + n.y*p0.y + n.z*p0.z);

   /*
      Evaluate the equation of the plane for each vertex
      If side < 0 then it is on the side to be retained
      else it is to be clipped
   */
   side[0] = A*p[0].x + B*p[0].y + C*p[0].z + D;
   side[1] = A*p[1].x + B*p[1].y + C*p[1].z + D;
   side[2] = A*p[2].x + B*p[2].y + C*p[2].z + D;

   /* Are all the vertices on the clipped side */
   if (side[0] >= 0 && side[1] >= 0 && side[2] >= 0)
      return(0);

   /* Are all the vertices on the non-clipped side */
   if (side[0] <= 0 || side[1] <= 0 || side[2] <= 0)
      return(1);

   return(0);
}

int ClipSimple3(XYZ *p, float n[3], float p0[3])
{
   float A,B,C,D;
   float l;
   float side[3];

   /*
      Determine the equation of the plane as
      Ax + By + Cz + D = 0
   */
   l = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
   A = n[0] / l;
   B = n[1] / l;
   C = n[2] / l;
   D = -(n[0]*p0[0] + n[1]*p0[1] + n[2]*p0[2]);

   /*
      Evaluate the equation of the plane for each vertex
      If side < 0 then it is on the side to be retained
      else it is to be clipped
   */
   side[0] = A*p[0].x + B*p[0].y + C*p[0].z + D;
   side[1] = A*p[1].x + B*p[1].y + C*p[1].z + D;
   side[2] = A*p[2].x + B*p[2].y + C*p[2].z + D;

   /* Are all the vertices on the clipped side */
   if (side[0] >= 0 && side[1] >= 0 && side[2] >= 0)
      return(0);

   /* Are all the vertices on the non-clipped side */
   if (side[0] <= 0 || side[1] <= 0 || side[2] <= 0)
      return(1);

   return(0);
}


/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ */

#define COPY_XYZ(dst, src) \
do { \
  dst.x = src.x; \
  dst.y = src.y; \
  dst.z = src.z; \
} while(0)

#define COPY_V3(dst, src) \
do { \
  (dst)[0] = (src)[0]; \
  (dst)[1] = (src)[1]; \
  (dst)[2] = (src)[2]; \
} while(0)

#define COPY_M33(dst, src) \
do { \
  (dst)[0][0] = (src)[0][0]; (dst)[0][1] = (src)[0][1]; (dst)[0][2] = (src)[0][2]; \
  (dst)[1][0] = (src)[1][0]; (dst)[1][1] = (src)[1][1]; (dst)[1][2] = (src)[1][2]; \
  (dst)[2][0] = (src)[2][0]; (dst)[2][1] = (src)[2][1]; (dst)[2][2] = (src)[2][2]; \
} while(0)

/* 
   WARNING - WARNING - WARNING -> extremely dengerous MACRO !!!
   ------------------------------------------------------------

   1.) It is meant only for cryIsosurfClip function !!!
   2.) Any change of variables names or structure REQUIRES an update 
       of this Macro !!!
*/
#define NEW_TRIANGLE(i0, i1, i2) \
tc++; \
/*cryAux_DEBUG("NEW_TRIANGLE");*/ \
\
out->triangl_status[nt] = 1; \
\
out->tri2verIN[nt].ver_status[0] = 1; \
out->tri2verIN[nt].ver_status[1] = 1; \
out->tri2verIN[nt].ver_status[2] = 1; \
\
out->tri2verIN[nt].ver_ind[0] = nv; \
out->tri2verIN[nt].ver_ind[1] = nv+1; \
out->tri2verIN[nt].ver_ind[2] = nv+2; \
nt++; \
\
COPY_XYZ(out->vertex[nv], p[i0]); \
COPY_XYZ(out->vertex_orig[nv], p[i0]); \
COPY_XYZ(out->normal[nv], nml[i0]); \
nv++;\
COPY_XYZ(out->vertex[nv], p[i1]); \
COPY_XYZ(out->vertex_orig[nv], p[i1]); \
COPY_XYZ(out->normal[nv], nml[i1]); \
nv++;\
COPY_XYZ(out->vertex[nv], p[i2]); \
COPY_XYZ(out->vertex_orig[nv], p[i2]); \
COPY_XYZ(out->normal[nv], nml[i2]); \
nv++

/*
   ************************************************************************
   NOTE about cryIsosurfClip:

   This is a modified routine (adapted for present pupose) of the 
   original one which is due to P. Bourke; see Modelling on
   http://www.swin.edu.au/astronomy/pbourke/modelling/
   ************************************************************************
*/

/*!
 *-----------------------------------------------------------------------------
 *
 * cryIsosurfClip -- \brief clip an isosurface by desired plane
 *
 * \description  
 *    This function clips an isosurface by desired plane. The memory
 *    of input isosurface is kept intact, hence take care of possible
 *    memory deallocation in the calling routine.
 *
 * \warning
 *     This routine is very memory allocation demanding, since for
 *     every new clipping plane new isosurface is allocated. In the
 *     future this can be done as: create a driver routine for all
 *     clipping planes, then make a vertex and triangle variables and
 *     associated "exist/non-exist" flag. Then clip all the planes,
 *     and remove all the non-existent entities. The mapping between
 *     new and old vertex ID's can be made by kind of lookup table.
 *
 *     Otherwise, this can be beautifully done with linked list, if
 *     triangle is out, simply skip it !!!
 *
 * \return
 *    Newly created clipped isosurface.
 *
 * \sideeffect
 *    None.
 *
 *----------------------------------------------------------------------------- */

ISOSURFACE *cryIsosurfClip(ISOSURFACE *in,  /*!< isosurface to clip */
			   XYZ        n,  /*!< normal of clipping plane */
			   XYZ        p0) /*!< point on clipping plane  */
{
  register int it, nt=0, nv=0, tc=0;
  int nv_max = (int) (1.1 * (float) in->nvertex);
  int nt_max = (int) (1.1 * (float) in->ntriangl);    
  ISOSURFACE *out;
  XYZ p[4], nml[4], q, qn;
  float A,B,C,D;
  float side[3];
  
  /*  allocate out ISOSURFACE and initialize it !!! */
  out =  (ISOSURFACE *) calloc (1, sizeof(ISOSURFACE));
  if (!out) xcError("allocation error for ISOSURFACE");

  /*
  fprintf(stderr,"Number of IN triangles:  %d\n", in->ntriangl);
  fprintf(stderr,"Number of IN vertices:   %d\n", in->nvertex);
  */

  out->index          = in->index;
  out->smooth_nstep   = in->smooth_nstep;
  out->smooth_weight  = in->smooth_weight;
  out->triangl_status = (int *)        xcRealloc((void*) out->triangl_status, (size_t) nt_max * sizeof(int) );
  out->vertex         = (XYZ *)        xcRealloc((void*) out->vertex,         (size_t) nv_max * sizeof(VERTEX) );
  out->vertex_orig    = (XYZ *)        xcRealloc((void*) out->vertex_orig,    (size_t) nv_max * sizeof(VERTEX) );
  out->color          = (float (*)[4]) xcRealloc((void*) out->color,          (size_t) nv_max * sizeof(float [4]) );
  out->normal         = (XYZ *)        xcRealloc((void*) out->normal,         (size_t) nv_max * sizeof(VERTEX) );
  out->tri2verIN      = (TRIG_INFO *)  xcRealloc((void*) out->tri2verIN,      (size_t) nt_max * sizeof(TRIG_INFO) );
  out->nver2triIN     = (int *)        xcRealloc((void*) out->nver2triIN,     (size_t) nt_max*3 * sizeof(int) );
  out->ver2triIN      = (VERT2TRIG_INFO (*)[MAX_VERTEX_COORDINATION]) 
    xcRealloc((void*) out->ver2triIN, (size_t )  nt_max*3 * sizeof(VERT2TRIG_INFO [MAX_VERTEX_COORDINATION]) );

  /*
    Determine the equation of the plane as:
    Ax + By + Cz + D = 0

    l = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);
    A = n.x / l;
    B = n.y / l;
    C = n.z / l;
  */

  A = n.x;
  B = n.y;
  C = n.z;
  D = -(n.x*p0.x + n.y*p0.y + n.z*p0.z);

  /* MAPPING:
     
  cryIso->triangl[it][ipol] = ISO->tri2verIN[it].ver_ind[ipol]
  
  */

  for (it=0; it<in->ntriangl; it++) {
    /*
      check for reallocation
    */
    if ( (nt+2) >= nt_max ) {
      nt_max = (int) (1.5 * (float) nt_max);
      out->tri2verIN = (TRIG_INFO *) xcRealloc((void*) out->tri2verIN, (size_t) nt_max * sizeof(TRIG_INFO) );
    }
    if ( (nv+6) >= nv_max ) {
      nv_max = (int) (1.5 * (float) nv_max);
      out->vertex      = (XYZ *)        xcRealloc((void*) out->vertex,      (size_t) nv_max * sizeof(VERTEX) );
      out->vertex_orig = (XYZ *)        xcRealloc((void*) out->vertex_orig, (size_t) nv_max * sizeof(VERTEX) );
      out->color       = (float (*)[4]) xcRealloc((void*) out->color,       (size_t) nv_max * sizeof(float [4]) );
      out->normal      = (XYZ *)        xcRealloc((void*) out->normal,      (size_t) nv_max * sizeof(VERTEX) );
    }

    /*      
      Clip a 3 vertex facet in place. The 3 point facet is defined by
      vertices p[0],p[1],p[2]. There must be a fourth point ("p[3]")
      as a 4 point facet may result. The normal to the plane is "n". A
      point on the plane is "p0". The side of the plane containing the
      normal is clipped away. Return the number of vertices in the
      clipped polygon 
    */

    COPY_XYZ (p[0], in->vertex[ in->tri2verIN[it].ver_ind[0] ]);
    COPY_XYZ (p[1], in->vertex[ in->tri2verIN[it].ver_ind[1] ]);
    COPY_XYZ (p[2], in->vertex[ in->tri2verIN[it].ver_ind[2] ]);

    COPY_XYZ (nml[0], in->normal[ in->tri2verIN[it].ver_ind[0] ]);
    COPY_XYZ (nml[1], in->normal[ in->tri2verIN[it].ver_ind[1] ]);
    COPY_XYZ (nml[2], in->normal[ in->tri2verIN[it].ver_ind[2] ]);

    /*
      Evaluate the equation of the plane for each vertex. If side < 0
      then it is on the side to be retained else it is to be clipped 
    */
    side[0] = A*p[0].x + B*p[0].y + C*p[0].z + D;
    side[1] = A*p[1].x + B*p[1].y + C*p[1].z + D;
    side[2] = A*p[2].x + B*p[2].y + C*p[2].z + D;
    /*cryAux_DEBUG("Clipping plane: (%f, %f, %f)\n", side[0], side[1], side[2]);*/

    /*
      now find-out the CLIPPING situation !!!
    */

    if (side[0] >= 0 && side[1] >= 0 && side[2] >= 0) {
      /* All the vertices are on the clipped side */
      continue;
    } 
    else if (side[0] <= 0 && side[1] <= 0 && side[2] <= 0) {
      /* All the vertices are on the not-clipped side */
      /* new triangle is 0,1,2 */
      NEW_TRIANGLE(0,1,2);
    } 
    else if (side[0] > 0 && side[1] < 0 && side[2] < 0) {
      /* p0 is the only point on the clipped side */
      /* new triangles are 0,1,2 & 0,2,3 */
      float s020 = side[0] / (side[2] - side[0]);
      float s010 = side[0] / (side[1] - side[0]);
      q.x = p[0].x - s020 * (p[2].x - p[0].x);
      q.y = p[0].y - s020 * (p[2].y - p[0].y);
      q.z = p[0].z - s020 * (p[2].z - p[0].z);
      qn.x = nml[0].x - s020 * (nml[2].x - nml[0].x);
      qn.y = nml[0].y - s020 * (nml[2].y - nml[0].y);
      qn.z = nml[0].z - s020 * (nml[2].z - nml[0].z);
      p[3] = q;
      nml[3] = qn;
      
      q.x = p[0].x - s010 * (p[1].x - p[0].x);
      q.y = p[0].y - s010 * (p[1].y - p[0].y);
      q.z = p[0].z - s010 * (p[1].z - p[0].z);
      qn.x = nml[0].x - s010 * (nml[1].x - nml[0].x);
      qn.y = nml[0].y - s010 * (nml[1].y - nml[0].y);
      qn.z = nml[0].z - s010 * (nml[1].z - nml[0].z);
      p[0] = q;
      nml[0] = qn;

      NEW_TRIANGLE(0,1,2);
      NEW_TRIANGLE(0,2,3);
    }
    else if (side[1] > 0 && side[0] < 0 && side[2] < 0) {
      /* p1 is the only point on the clipped side */
      /* new triangles are 0,1,2 & 0,2,3 */
      float s121 = side[1] / (side[2] - side[1]);
      float s101 = side[1] / (side[0] - side[1]);
      p[3] = p[2];
      nml[3] = nml[2];
      q.x = p[1].x - s121 * (p[2].x - p[1].x);
      q.y = p[1].y - s121 * (p[2].y - p[1].y);
      q.z = p[1].z - s121 * (p[2].z - p[1].z);
      qn.x = nml[1].x - s121 * (nml[2].x - nml[1].x);
      qn.y = nml[1].y - s121 * (nml[2].y - nml[1].y);
      qn.z = nml[1].z - s121 * (nml[2].z - nml[1].z);
      p[2] = q;
      nml[2] = qn;
      q.x = p[1].x - s101 * (p[0].x - p[1].x);
      q.y = p[1].y - s101 * (p[0].y - p[1].y);
      q.z = p[1].z - s101 * (p[0].z - p[1].z);
      qn.x = nml[1].x - s101 * (nml[0].x - nml[1].x);
      qn.y = nml[1].y - s101 * (nml[0].y - nml[1].y);
      qn.z = nml[1].z - s101 * (nml[0].z - nml[1].z);
      p[1] = q;
      nml[1] = qn;

      NEW_TRIANGLE(0,1,2);
      NEW_TRIANGLE(0,2,3);      
    }
    else if (side[2] > 0 && side[0] < 0 && side[1] < 0) {
      /* p2 is the only point on the clipped side */
      /* new triangles are 0,1,2 & 0,2,3 */      
      float s202 = side[2] / (side[0] - side[2]);
      float s212 = side[2] / (side[1] - side[2]);
      q.x = p[2].x - s202 * (p[0].x - p[2].x);
      q.y = p[2].y - s202 * (p[0].y - p[2].y);
      q.z = p[2].z - s202 * (p[0].z - p[2].z);
      qn.x = nml[2].x - s202 * (nml[0].x - nml[2].x);
      qn.y = nml[2].y - s202 * (nml[0].y - nml[2].y);
      qn.z = nml[2].z - s202 * (nml[0].z - nml[2].z);
      p[3] = q;
      nml[3] = qn;
      q.x = p[2].x - s212 * (p[1].x - p[2].x);
      q.y = p[2].y - s212 * (p[1].y - p[2].y);
      q.z = p[2].z - s212 * (p[1].z - p[2].z);
      qn.x = nml[2].x - s212 * (nml[1].x - nml[2].x);
      qn.y = nml[2].y - s212 * (nml[1].y - nml[2].y);
      qn.z = nml[2].z - s212 * (nml[1].z - nml[2].z);
      p[2] = q;
      nml[2] = qn;

      NEW_TRIANGLE(0,1,2);
      NEW_TRIANGLE(0,2,3);
    }
    else if (side[0] < 0 && side[1] > 0 && side[2] > 0) {
      /* p0 is the only point on the not-clipped side */
      /* new triangle is 0,1,2 */
      float s010 = side[0] / (side[1] - side[0]);
      float s020 = side[0] / (side[2] - side[0]);
      q.x = p[0].x - s010 * (p[1].x - p[0].x);
      q.y = p[0].y - s010 * (p[1].y - p[0].y);
      q.z = p[0].z - s010 * (p[1].z - p[0].z);
      qn.x = nml[0].x - s010 * (nml[1].x - nml[0].x);
      qn.y = nml[0].y - s010 * (nml[1].y - nml[0].y);
      qn.z = nml[0].z - s010 * (nml[1].z - nml[0].z);
      p[1] = q;
      nml[1] = qn;
      q.x = p[0].x - s020 * (p[2].x - p[0].x);
      q.y = p[0].y - s020 * (p[2].y - p[0].y);
      q.z = p[0].z - s020 * (p[2].z - p[0].z);
      qn.x = nml[0].x - s020 * (nml[2].x - nml[0].x);
      qn.y = nml[0].y - s020 * (nml[2].y - nml[0].y);
      qn.z = nml[0].z - s020 * (nml[2].z - nml[0].z);
      p[2] = q;
      nml[2] = qn;

      NEW_TRIANGLE(0,1,2);
    }    
    else if (side[1] < 0 && side[0] > 0 && side[2] > 0) {
      /* p1 is the only point on the not-clipped side */
      /* new triangle is 0,1,2 */
      float s101 = side[1] / (side[0] - side[1]);
      float s121 = side[1] / (side[2] - side[1]);
      q.x = p[1].x - s101 * (p[0].x - p[1].x);
      q.y = p[1].y - s101 * (p[0].y - p[1].y);
      q.z = p[1].z - s101 * (p[0].z - p[1].z);
      qn.x = nml[1].x - s101 * (nml[0].x - nml[1].x);
      qn.y = nml[1].y - s101 * (nml[0].y - nml[1].y);
      qn.z = nml[1].z - s101 * (nml[0].z - nml[1].z);
      p[0] = q;
      nml[0] = qn;
      q.x = p[1].x - s121 * (p[2].x - p[1].x);
      q.y = p[1].y - s121 * (p[2].y - p[1].y);
      q.z = p[1].z - s121 * (p[2].z - p[1].z);
      qn.x = nml[1].x - s121 * (nml[2].x - nml[1].x);
      qn.y = nml[1].y - s121 * (nml[2].y - nml[1].y);
      qn.z = nml[1].z - s121 * (nml[2].z - nml[1].z);
      p[2] = q;
      nml[2] = qn;

      NEW_TRIANGLE(0,1,2);
    }
    else if (side[2] < 0 && side[0] > 0 && side[1] > 0) {
      /* p2 is the only point on the not-clipped side */
      /* new triangle is 0,1,2 */
      float s212 = side[2] / (side[1] - side[2]);
      float s202 = side[2] / (side[0] - side[2]);
      q.x = p[2].x - s212 * (p[1].x - p[2].x);
      q.y = p[2].y - s212 * (p[1].y - p[2].y);
      q.z = p[2].z - s212 * (p[1].z - p[2].z);
      qn.x = nml[2].x - s212 * (nml[1].x - nml[2].x);
      qn.y = nml[2].y - s212 * (nml[1].y - nml[2].y);
      qn.z = nml[2].z - s212 * (nml[1].z - nml[2].z);
      p[1] = q;
      nml[1] = qn;
      q.x = p[2].x - s202 * (p[0].x - p[2].x);
      q.y = p[2].y - s202 * (p[0].y - p[2].y);
      q.z = p[2].z - s202 * (p[0].z - p[2].z);
      qn.x = nml[2].x - s202 * (nml[0].x - nml[2].x);
      qn.y = nml[2].y - s202 * (nml[0].y - nml[2].y);
      qn.z = nml[2].z - s202 * (nml[0].z - nml[2].z);
      p[0] = q;
      nml[0] = qn;

      NEW_TRIANGLE(0,1,2);
    }
  }

  out->ntriangl = nt;
  out->nvertex  = nv;

  /*
  fprintf(stderr,"Number of OUT triangles: %d\n", out->ntriangl);
  fprintf(stderr,"Number of OUT vertices:  %d\n\n", out->nvertex);
  fflush(stderr);
  */

  /*
    out->tri2verIN  = (TRIG_INFO *)  xcRealloc((void*) out->tri2verIN,      (size_t) nt * sizeof(TRIG_INFO) );
    out->vertex     = (XYZ *)        xcRealloc((void*) out->vertex,         (size_t) nv * sizeof(VERTEX) );
  */

  return out;
}
