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
 * Source: $XCRYSDEN_TOPDIR/C/polygonise_auxil.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <stdio.h> 
#include <math.h>
#include <stdlib.h> 
#include "struct.h"
#include "isosurf.h"
#include "xcfunc.h"

/* 
   BEWARE: MINTOL_F should not be lower that 1e-6, otherwise it won't work
*/
#define EQUAL_VERTICES(a,b)    ( fabs(a.x-b.x) < MINTOL_F && fabs(a.y-b.y) < MINTOL_F && fabs(a.z-b.z) < MINTOL_F )
#define GET_INDEX(i,j,k,N,M,L) ( (k + L*(j + M*i)) )

XYZ TriangleWeightNormalS(XYZ p0, XYZ p1, XYZ p2, float sign);
XYZ TriangleWeightNormal(XYZ p0, XYZ p1, XYZ p2);
float distXYZ(XYZ v);

void
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
		int           *nvertex)
{
  register int nt, i3, ii, jj, kk, l, j3, ver_ind, nvt, ix1, iy1, iz1;
  int ijk, t0l;
  IJK ind;
		
  ix1=grd.nx-1;
  iy1=grd.ny-1;
  iz1=grd.nz-1;

  /*
    nt  ... number of all triangles
    nvt ... number of triangles that belong to a given vertex 
  */

  /* VERTEX RECOGNITION */
  ver_ind = 0;

  for(nt=0; nt<ntriangl; nt++) {

    triangl_status[nt] = XC_TRUE;
    ind.i=tri2verIN[nt].i;
    ind.j=tri2verIN[nt].j;
    ind.k=tri2verIN[nt].k;

    for(i3=0; i3<3; i3++) {

      nvt = 0;

      if (!tri2verIN[nt].ver_status[i3]) continue;
      
      tri2verIN[nt].ver_ind[i3]    = ver_ind;
      tri2verIN[nt].ver_status[i3] = XC_FALSE;
      
      vertex[ver_ind]        = triangles[nt].p[i3];
      vertex_status[ver_ind] = XC_TRUE;

      ver2triIN[ver_ind][nvt  ].t = nt;
      ver2triIN[ver_ind][nvt++].v = i3;
      
      for(ii=ind.i-1; ii<=ind.i+1; ii++) {
	if (ii<0 || ii>=ix1) continue;  
	
	for(jj=ind.j-1; jj<=ind.j+1; jj++) {
	  if (jj<0 || jj>=iy1) continue;
	  
	  for(kk=ind.k-1; kk<=ind.k+1; kk++) {
	    if (kk<0 || kk>=iz1) continue;
	  
	    ijk=GET_INDEX(ii,jj,kk, ix1, iy1, iz1);
      
	    for(l=0; l<grid_info[ijk].nt; l++) 
	      {
		t0l=grid_info[ijk].t0+l;

		for(j3=0; j3<3;j3++) 
		  {
		    if( EQUAL_VERTICES(triangles[nt].p[i3],triangles[t0l].p[j3]) ) 
		      {
			if ( nt != t0l || j3 != i3 ) 
			  {
			    ver2triIN[ver_ind][nvt  ].t = t0l;
			    ver2triIN[ver_ind][nvt++].v = j3;
			    tri2verIN[t0l].ver_status[j3] = XC_FALSE;
			    tri2verIN[t0l].ver_ind[j3]    = ver_ind;
			  } 
		      }

		    /* if we are at boundary the EQUAL_VERTICES should consider
		       the translational periodicity !!!
		       ....INSERT CODE HERE.... */
		    
		  } /* j3 */
	      } /* l */
	  } /* kk */
	} /* jj */
      } /* ii */
      nver2triIN[ver_ind] = nvt;
      ver_ind++;
    } /* i3 */
  } /* nt */
  
  *nvertex=ver_ind;
}


void
SurfSmoothing(int            nstep, 
	      float          weight,
	      XYZ            *vertex,
	      int            *vertex_status,
	      VERT2TRIG_INFO (*ver2triIN)[MAX_VERTEX_COORDINATION],
	      TRIG_INFO      *tri2verIN,
	      int            *nver2triIN,
	      int            nvertex)
{
  register int is, i, j;
  int i0, i1;    

  /*
    X(n+1) = (1 + w*L) X(n)

    L = 1 / E sum{vi - vj/|vi - vj|}
  */
  float E, u, w;
  XYZ L, dv;
  XYZ *aux; 

  if (nstep==0) return;

  aux = (XYZ*) malloc(sizeof(XYZ)*nvertex);

  /* 
     k = 1/w + 1/u 

     0.1 - 1/w = 1/u --> u = 1 / (0.1 - 1/w)
  */
  u = 1.0 / (0.1 - 1.0 / weight);

  weight *= 0.5; /* because we count each vertex twice */
  u      *= 0.5;

  /* SURFACE SMOOTHING */
  nstep *= 2; /* shrinking + expanding steps */
  for(is=0; is<nstep; is++)
    {      
      for (i=0; i<nvertex; i++) aux[i] = vertex[i];

      for (i=0; i<nvertex; i++) {
	/* if someday surface-simplification will be implemented !!!*/
	if(!vertex_status[i]) continue;
	
	E = L.x = L.y = L.z = 0.0;

	for(j=0; j<nver2triIN[i]; j++) {
	  
	  if (ver2triIN[i][j].v == 0) {
	    i0 = tri2verIN[ ver2triIN[i][j].t ].ver_ind[1];
	    i1 = tri2verIN[ ver2triIN[i][j].t ].ver_ind[2];
	  } else if (ver2triIN[i][j].v == 1) {
	    i0 = tri2verIN[ ver2triIN[i][j].t ].ver_ind[0];
	    i1 = tri2verIN[ ver2triIN[i][j].t ].ver_ind[2];
	  } else {
	    /* if (ver2triIN[i][j].v == 2) */
	    i0 = tri2verIN[ ver2triIN[i][j].t ].ver_ind[0];
	    i1 = tri2verIN[ ver2triIN[i][j].t ].ver_ind[1];
	  } 
	  
	  dv.x = aux[i0].x - aux[i].x;
	  dv.y = aux[i0].y - aux[i].y;
	  dv.z = aux[i0].z - aux[i].z;

	  /*
	  e = MAX(sqrt(dv.x*dv.x + dv.y*dv.y + dv.z*dv.z),1e-3);
	  E   += e;
	  L.x += dv.x / e;
	  L.y += dv.y / e;
	  L.z += dv.z / e;
	  */
	  
	  L.x += dv.x;
	  L.y += dv.y;
	  L.z += dv.z;	  
	  
	  dv.x = aux[i1].x - aux[i].x;
	  dv.y = aux[i1].y - aux[i].y;
	  dv.z = aux[i1].z - aux[i].z;
	  
	  /*
	  e = MAX(sqrt(dv.x*dv.x + dv.y*dv.y + dv.z*dv.z),1e-3);
	  E   += e;
	  L.x += dv.x / e;
	  L.y += dv.y / e;
	  L.z += dv.z / e;
	  */

	    L.x += dv.x;
	    L.y += dv.y;
	    L.z += dv.z;
	}
	
	/* insert here some condition, 
	   that prevents modifying the border atoms */

	E=(float)nver2triIN[i];
	if (E > MINTOL) {
	  /*fprintf(stderr, "Vertex smoothing (ntriangl = %2d) : %5.6f %5.6f %5.6f  --> ",
		  nver2triIN[i],  vertex[i].x, vertex[i].y, vertex[i].y);
	  */
	  
	  /*if (E<1.0) E = (float)nver2triIN[i];*/
	  
	  if ( is%2 == 0 ) {
	    w = weight;  /* shrinking step: use "weight" */
	  } else {
	    w = u; /* expansion step */
	  }
	  vertex[i].x += w * (1.0/E * L.x);  
	  vertex[i].y += w * (1.0/E * L.y);  
	  vertex[i].z += w * (1.0/E * L.z);      
	  
	  /*
	  fprintf(stderr, "%5.2f %5.2f %5.2f ; (E=%3f)\n",
		  vertex[i].x, vertex[i].y, vertex[i].y, E);
	  fprintf(stderr, "Vertex smoothing (E=%f, L: %f %f %f)\n", E, 
		  L.x, L.y, L.z);
	  */
	}
      }
    }
  free(aux);
}


void
SurfNmlAver(int       ntriangl,
	    int       nvertex,
	    XYZ       *nml,
	    XYZ       *vertex,	    
	    TRIG_INFO *tri2verIN,
	    float     *sign,
	    int       is_tetrahedral)
{
  register int i, j3;
  XYZ norm;

  /*fprintf(stderr,"        SurfNmlAver used,             normals address=%ld\n", nml);*/

  /* NORMAL AVERAGING */
  for(i=0; i<nvertex; i++)
    nml[i].x = nml[i].y = nml[i].z = 0.0;
  
  if(is_tetrahedral) {

    /* for tetrahedral the sign should be allocated ! */

    if ( sign == NULL ) abort();

    for(i=0; i<ntriangl; i++) {
      /* next line for the sake of surface simplification */
      /*if(!triangl_status[i]) continue;*/
      
      norm = 
	TriangleWeightNormalS(vertex[ tri2verIN[i].ver_ind[0] ], 
			      vertex[ tri2verIN[i].ver_ind[1] ], 
			      vertex[ tri2verIN[i].ver_ind[2] ], 
			      sign[i] );      
      for(j3=0; j3<3; j3++) {
	nml[ tri2verIN[i].ver_ind[j3] ].x += norm.x;
	nml[ tri2verIN[i].ver_ind[j3] ].y += norm.y;
	nml[ tri2verIN[i].ver_ind[j3] ].z += norm.z;
      }
    }
  } else {
    for(i=0; i<ntriangl; i++) {
      /* next line for the sake of surface simplification */
      /*if(!triangl_status[i]) continue;*/
      
      norm = 
	TriangleWeightNormal(vertex[ tri2verIN[i].ver_ind[0] ], 
			     vertex[ tri2verIN[i].ver_ind[1] ], 
			     vertex[ tri2verIN[i].ver_ind[2] ]);      
      for(j3=0; j3<3; j3++) {
	nml[ tri2verIN[i].ver_ind[j3] ].x += norm.x;
	nml[ tri2verIN[i].ver_ind[j3] ].y += norm.y;
	nml[ tri2verIN[i].ver_ind[j3] ].z += norm.z;
      }
    }
  }

  for(i=0; i<nvertex; i++)
    normalizeXYZ( &nml[i] );
}



/* Calculate normal; take vector v1 & v2 and return it's vector product 
 *
 * returned normal is normalized
 * ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 */
XYZ
TriangleWeightNormalS(XYZ p0, XYZ p1, XYZ p2, float sign)
{
  float dis12;
  /* Normals; take cross product of [i+1,i] & [i+2,i+1] */

  XYZ v1, v2, nml;

  v1.x = p1.x - p0.x;
  v1.y = p1.y - p0.y;
  v1.z = p1.z - p0.z;

  v2.x = p2.x - p0.x;
  v2.y = p2.y - p0.y;
  v2.z = p2.z - p0.z;

  nml.x = sign * (v1.y * v2.z - v2.y * v1.z);
  nml.y = sign * (v1.z * v2.x - v2.z * v1.x);
  nml.z = sign * (v1.x * v2.y - v2.x * v1.y);


  dis12=distXYZ(nml);
  if (dis12 > MINTOL) {
    nml.x /= dis12;
    nml.y /= dis12;
    nml.z /= dis12;
  } else { 
    nml.x = .0;
    nml.y = .0;
    nml.z = .0;
  }
  
  return(nml);
}
XYZ
TriangleWeightNormal(XYZ p0, XYZ p1, XYZ p2)
{
  float dis12;
  /* Normals; take cross product of [i+1,i] & [i+2,i+1] */

  XYZ v1, v2, nml;

  v1.x = p1.x - p0.x;
  v1.y = p1.y - p0.y;
  v1.z = p1.z - p0.z;

  v2.x = p2.x - p0.x;
  v2.y = p2.y - p0.y;
  v2.z = p2.z - p0.z;

  nml.x = v1.y * v2.z - v2.y * v1.z;
  nml.y = v1.z * v2.x - v2.z * v1.x;
  nml.z = v1.x * v2.y - v2.x * v1.y;


  dis12=distXYZ(nml);
  if (dis12 > MINTOL) {    
    nml.x /= (dis12*dis12);
    nml.y /= (dis12*dis12);
    nml.z /= (dis12*dis12);    
  } else { 
    nml.x = 0.0;
    nml.y = 0.0;
    nml.z = 0.0;
  }
  
  return(nml);
}
XYZ
VertexNormal(XYZ v1, XYZ v2)
{
  XYZ nml;

  nml.x = v1.y * v2.z - v2.y * v1.z;
  nml.y = v1.z * v2.x - v2.z * v1.x;
  nml.z = v1.x * v2.y - v2.x * v1.y;

  normalizepvfv( (float *) &nml );
  return(nml);
}

float 
distXYZ(XYZ v)
{
  float dis;
  dis = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
  /*printf("dis=%f\n",dis);*/
  return(dis);
}

int 
normalizeXYZ(XYZ *v)
{
  float dis;

  dis = distXYZ(*v);
  if ( dis < 1.e-12 ) {
    (*v).x = 0.0;
    (*v).y = 0.0;
    (*v).z = 0.0;
    return 0;
  }
  (*v).x /= dis;
  (*v).y /= dis;
  (*v).z /= dis;
  
  return 1;
}

