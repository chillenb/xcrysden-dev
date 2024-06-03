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
 * Source: $XCRYSDEN_TOPDIR/C/gridNormals.c
 * ------                                                                    *
 * Copyright (c) 2009 by Anton Kokalj                                        *
 *****************************************************************************

*/

#include <stdio.h>
#include <math.h>
#include "struct.h"
#include "isosurf.h"
#include "xcfunc.h"

/* note that M = N-1, where N is the number of points in the grid direction */
#define sign(x) ( (x) >= 0 ? +1 : -1 )
#define PBC_INDEX( i, M ) ( ((i)>=(M) || (i) < 0) ? (i)-sign(i)*(M) : (i) )

#define XYZ_ApTBmA(dst,a,t,b) \
do { \
  dst.x = a.x + t*(b.x-a.x); \
  dst.y = a.y + t*(b.y-a.y); \
  dst.z = a.z + t*(b.z-a.z); \
 }  while(0)

#define DOT_V3(a,b)  ( (a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2] )

/* static void  */
/* metricTensor33f(const float VEC_in[3][3], float T_out[3][3])  */
/* { */
/*   T_out[0][0] = DOT_V3 (VEC_in[0], VEC_in[0]); */
/*   T_out[1][1] = DOT_V3 (VEC_in[1], VEC_in[1]); */
/*   T_out[2][2] = DOT_V3 (VEC_in[2], VEC_in[2]); */
  
/*   T_out[0][1] = T_out[1][0] = DOT_V3 (VEC_in[0], VEC_in[1]); */
/*   T_out[0][2] = T_out[2][0] = DOT_V3 (VEC_in[0], VEC_in[2]); */
/*   T_out[1][2] = T_out[2][1] = DOT_V3 (VEC_in[1], VEC_in[2]); */
/* } */


void
gridNormals(GRIDVERTEX ***g, GRID grd)
{
  /* we implicitly assume that grid describes a periodic object, but
     we take care of reducible points at the right-cell boundaries */

  register int ix, iy, iz, ixp, ixm, iyp, iym, izp, izm;
  int ixN = grd.nx - 1;
  int iyN = grd.ny - 1;
  int izN = grd.nz - 1;

  for (ix = 0; ix < ixN; ix++) {
    for (iy = 0; iy < iyN; iy++) {
      for (iz = 0; iz < izN; iz++)
	{
	  ixp = ix+1;
	  ixm = ix-1;
	  
	  iyp = iy+1;
	  iym = iy-1;
	  
	  izp = iz+1;
	  izm = iz-1;
	  
	  /* 
	     gradient (normal) in crystal coordinates:
	  
	     grad_{a,b,c} g = (dg/da, dg/db, dg/dc) 
	  */
	  
	  g[ix][iy][iz].n.x = 0.5 * (g[ PBC_INDEX(ixp,ixN) ][iy][iz].val - g[ PBC_INDEX(ixm,ixN) ][iy][iz].val);
	  g[ix][iy][iz].n.y = 0.5 * (g[ix][ PBC_INDEX(iyp,iyN) ][iz].val - g[ix][ PBC_INDEX(iym,iyN) ][iz].val);
	  g[ix][iy][iz].n.z = 0.5 * (g[ix][iy][ PBC_INDEX(izp,izN) ].val - g[ix][iy][ PBC_INDEX(izm,izN) ].val);
	}
    }
  }
}


void gradient_SurfNml(GRIDVERTEX ***g, GRID grd, XYZ spanning_vec[3], int ntriangl, TRIG_INFO *tri2verIN, int nvertex, XYZ *vertex, XYZ *normal)
{

  int it, iv, ig, i, j, k, i1, j1, k1;
  float dx, dy, dz, dis;
  XYZ gp[8], gn[8], gn_xyz[8]; /* vertices and normals of a given grid cube */
  XYZ h;
  int iN = grd.nx - 1;
  int jN = grd.ny - 1;
  int kN = grd.nz - 1;
  float lat[3][3], invM[3][3];

  for (i=0; i<3; i++) 
    {
      lat[i][0] = spanning_vec[i].x / (float) iN;
      lat[i][1] = spanning_vec[i].y / (float) jN;
      lat[i][2] = spanning_vec[i].z / (float) kN;
    }

  INVERT_M33(float, invM, lat);

  for(iv=0; iv<nvertex; iv++) 
    {
      normal[iv].x = normal[iv].y = normal[iv].z = 0.0;
    }

  for(it=0; it<ntriangl; it++) 
    {   
      i = tri2verIN[it].i;
      j = tri2verIN[it].j;
      k = tri2verIN[it].k;
      
      i1 = PBC_INDEX(i+1,iN);
      j1 = PBC_INDEX(j+1,jN);
      k1 = PBC_INDEX(k+1,kN);
      
      gp[0] = g[i ][j ][k ].p;
      gp[1] = g[i ][j1][k ].p;
      gp[2] = g[i1][j1][k ].p;
      gp[3] = g[i1][j ][k ].p;
      gp[4] = g[i ][j ][k1].p;
      gp[5] = g[i ][j1][k1].p;
      gp[6] = g[i1][j1][k1].p;
      gp[7] = g[i1][j ][k1].p;

      gn[0] = g[i ][j ][k ].n;
      gn[1] = g[i ][j1][k ].n;
      gn[2] = g[i1][j1][k ].n;
      gn[3] = g[i1][j ][k ].n;
      gn[4] = g[i ][j ][k1].n;
      gn[5] = g[i ][j1][k1].n;
      gn[6] = g[i1][j1][k1].n;
      gn[7] = g[i1][j ][k1].n;      

      for (iv=0; iv<3; iv++)
	{
	  h.x = h.y = h.z = 0.0;

	  for (ig=0; ig<8; ig++)
	    {
	      /* 
		 transform grad_{a,b,c} to grad_{x,y,z}, i.e.:
		 
		 grad_{x,y,z} = lat^-1 grad_{a,b,c}, where lat is a Bravais matrix
	      */

	      gn_xyz[ig].x = invM[0][0]*gn[ig].x + invM[0][1]*gn[ig].y + invM[0][2]*gn[ig].z;
	      gn_xyz[ig].y = invM[1][0]*gn[ig].x + invM[1][1]*gn[ig].y + invM[1][2]*gn[ig].z;
	      gn_xyz[ig].z = invM[2][0]*gn[ig].x + invM[2][1]*gn[ig].y + invM[2][2]*gn[ig].z;

	      /*
		cartesian distance between vertex and grid point "ig"
	      */
	      
	      dx = vertex[ tri2verIN[it].ver_ind[iv] ].x - gp[ig].x;
	      dy = vertex[ tri2verIN[it].ver_ind[iv] ].y - gp[ig].y;
	      dz = vertex[ tri2verIN[it].ver_ind[iv] ].z - gp[ig].z;
	      
	      dis = sqrt(dx*dx + dy*dy + dz*dz);
	      
	      if ( dis < MINTOL_F ) {
		h = gn_xyz[ig];
		break;
	      }
	      
	      /* don't need the normalization, only direction matters */
	      
	      h.x += gn_xyz[ig].x / dis;
	      h.y += gn_xyz[ig].y / dis;
	      h.z += gn_xyz[ig].z / dis;
	    }
	  
	  normalizeXYZ ( &h );
	  normal[ tri2verIN[it].ver_ind[iv] ].x += h.x;
	  normal[ tri2verIN[it].ver_ind[iv] ].y += h.y;
	  normal[ tri2verIN[it].ver_ind[iv] ].z += h.z;
	}
    }

  for(iv=0; iv<nvertex; iv++)
    {
      normalizeXYZ( &normal[iv] );
    }
}

void
gradient_SurfNmlNew(float origin[], 
		    float incr_vec[3][3], 
		    float ***g, 
		    GRID  gr, 
		    int   ntriangl, 
		    TRIG_INFO *tri2verIN, 
		    int  nvertex, 
		    XYZ  vertex[], 
		    XYZ  nml[])
{
  float invM[3][3], org[3], rver[3], iver[3];
  XYZ   hc, gc[8], grad[8];
  XYZ   u0, u1, u2, u3, v0, v1;
  int   i, j, k;
  int   it, ic, iv, ix, iy, iz, ixm, ixp, iym, iyp, izm, izp;
  int   ixN = gr.nx - 1;
  int   iyN = gr.ny - 1;
  int   izN = gr.nz - 1;
  int   cube[8][3] = {
    {0,0,0},
    {0,1,0},
    {1,1,0},
    {1,0,0},
    {0,0,1},
    {0,1,1},
    {1,1,1},
    {1,0,1}
  };

  /*fprintf(stderr,"gradient_SurfNmlNew used,             normals address=%ld\n", nml);*/
  
  for(iv=0; iv<nvertex; iv++) 
    {
      nml[iv].x = nml[iv].y = nml[iv].z = 0.0;
    }

  INVERT_M33(float, invM, incr_vec);

  for(it=0; it<ntriangl; it++) 
    {   
      i = tri2verIN[it].i;
      j = tri2verIN[it].j;
      k = tri2verIN[it].k;

      /* org^T = origin^T + (i,j,k)*incr_vec */

      org[0] = origin[0] + (float)i * incr_vec[0][0] + (float)j * incr_vec[1][0] + (float)k * incr_vec[2][0];
      org[1] = origin[1] + (float)i * incr_vec[0][1] + (float)j * incr_vec[1][1] + (float)k * incr_vec[2][1];
      org[2] = origin[2] + (float)i * incr_vec[0][2] + (float)j * incr_vec[1][2] + (float)k * incr_vec[2][2];

      for(ic=0; ic<8; ic++)
	{
	  ix = i + cube[ic][0];
	  iy = j + cube[ic][1];
	  iz = k + cube[ic][2];
	  
	  ixp = ix+1;
	  ixm = ix-1;
	  
	  iyp = iy+1;
	  iym = iy-1;
	  
	  izp = iz+1;
	  izm = iz-1;

	  /* grad_{a,b,c} g = (dg/da, dg/db, dg/dc) */
	  
	  grad[ic].x = 0.5 * (g[ PBC_INDEX(ixp,ixN) ][iy][iz] - g[ PBC_INDEX(ixm,ixN) ][iy][iz]);
	  grad[ic].y = 0.5 * (g[ix][ PBC_INDEX(iyp,iyN) ][iz] - g[ix][ PBC_INDEX(iym,iyN) ][iz]);
	  grad[ic].z = 0.5 * (g[ix][iy][ PBC_INDEX(izp,izN) ] - g[ix][iy][ PBC_INDEX(izm,izN) ]);

	  /* 
	     transform grad_{a,b,c} to grad_{x,y,z}, i.e.:
	     
	     grad_{x,y,z} = invM * grad_{a,b,c}, where invM is the inverse Bravais matrix
	  */

	  gc[ic].x =  invM[0][0]*grad[ic].x + invM[0][1]*grad[ic].y + invM[0][2]*grad[ic].z;
	  gc[ic].y =  invM[1][0]*grad[ic].x + invM[1][1]*grad[ic].y + invM[1][2]*grad[ic].z;
	  gc[ic].z =  invM[2][0]*grad[ic].x + invM[2][1]*grad[ic].y + invM[2][2]*grad[ic].z;
	}
	
      for (iv=0; iv<3; iv++)
	{
	  rver[0] = vertex[ tri2verIN[it].ver_ind[iv] ].x - org[0];
	  rver[1] = vertex[ tri2verIN[it].ver_ind[iv] ].y - org[1];
	  rver[2] = vertex[ tri2verIN[it].ver_ind[iv] ].z - org[2];
	  
	  /* iver^T = rver^T * invM */

	  iver[0] = rver[0] * invM[0][0] + rver[1] * invM[1][0] + rver[2] * invM[2][0];
	  iver[1] = rver[0] * invM[0][1] + rver[1] * invM[1][1] + rver[2] * invM[2][1];
	  iver[2] = rver[0] * invM[0][2] + rver[1] * invM[1][2] + rver[2] * invM[2][2];
	  /*fprintf(stderr,"iver[%d] = (%f,%f,%f)\n", iv, iver[0], iver[1], iver[2]);*/

	  for (ic=0; ic<3; ic++)
	    {
	      /* take care of round-off errors: we only use floats */
	      if ( iver[ic] < 0.0 ) iver[ic] = 0.0;
	      if ( iver[ic] > 1.0 ) iver[ic] = 1.0;
	    }
	  
	  /* 
	     TRI-LINEAR INTERPOLATION 
	  */
	  
	  /* z-dir */
	  XYZ_ApTBmA(u0, gc[0], iver[2], gc[4]); /* u0 = gc[0] + iver[2]*(gc[4]-gc[0]); */
	  XYZ_ApTBmA(u1, gc[1], iver[2], gc[5]); /* u1 = gc[1] + iver[2]*(gc[5]-gc[1]); */
	  XYZ_ApTBmA(u2, gc[2], iver[2], gc[6]); /* u2 = gc[2] + iver[2]*(gc[6]-gc[2]); */
	  XYZ_ApTBmA(u3, gc[3], iver[2], gc[7]); /* u3 = gc[3] + iver[2]*(gc[7]-gc[3]); */
	  
	  /* y-dir */
	  XYZ_ApTBmA(v0, u0, iver[1], u1); /* v0 = u0 + iver[1]*(u1-u0); */
	  XYZ_ApTBmA(v1, u3, iver[1], u2); /* v1 = u3 + iver[1]*(u2-u3); */
	  
	  /* x-dir */
	  XYZ_ApTBmA(hc, v0, iver[0], v1); /*  h = v0 + iver[0]*(v1-v0); */

	  /* 
	     WEIGHTED  NEAREST-NEIGHBOR INTERPOLATION 
	  */
	  
	  /* hc.x = hc.y = hc.z = 0.0; */
	  
	  /* for (ic=0; ic<8; ic++) */
	  /*   { */
	  /*     da = iver[0] - (float)cube[ic][0]; */
	  /*     db = iver[1] - (float)cube[ic][1]; */
	  /*     dc = iver[2] - (float)cube[ic][2]; */

	  /*     /\*  */
	  /* 	 transform distance-vector to Cartesian coordinate system: */

	  /* 	 (dx,dy,dz) = (da,db,dc) * incr_vec */
	  /*      *\/ */
	      
	  /*     dx = da * incr_vec[0][0] + db * incr_vec[1][0] + dc * incr_vec[2][0]; */
	  /*     dy = da * incr_vec[0][1] + db * incr_vec[1][1] + dc * incr_vec[2][1]; */
	  /*     dz = da * incr_vec[0][2] + db * incr_vec[1][2] + dc * incr_vec[2][2]; */

	  /*     dis = sqrt(dx*dx + dy*dy + dz*dz); */
	      
	  /*     if ( dis < MINTOL_F ) { */
	  /* 	hc = gc[ic]; */
	  /* 	break; */
	  /*     } */
	      
	  /*     /\* don't need the normalization, only direction matters *\/ */
	  /*     /\* w   += 1.0 / dis; *\/ */
	           
	  /*     hc.x += gc[ic].x / dis; */
	  /*     hc.y += gc[ic].y / dis; */
	  /*     hc.z += gc[ic].z / dis; */
	  /*   } */

	  /* fprintf(stderr, "nml[%d]=(%f,%f,%f); ver=(%f, %f, %f)\n", tri2verIN[it].ver_ind[iv], hc.x, hc.y, hc.z,  */
	  /* 	  vertex[ tri2verIN[it].ver_ind[iv] ].x, vertex[ tri2verIN[it].ver_ind[iv] ].y, vertex[ tri2verIN[it].ver_ind[iv] ].z); */
	  normalizeXYZ ( &hc );
	  nml[ tri2verIN[it].ver_ind[iv] ].x += hc.x;
	  nml[ tri2verIN[it].ver_ind[iv] ].y += hc.y;
	  nml[ tri2verIN[it].ver_ind[iv] ].z += hc.z;
	}
    }
  
  for(iv=0; iv<nvertex; iv++)
    {
      normalizeXYZ( &nml[iv] );
    }

  /* for(it=0; it<ntriangl; it++)  */
  /*   {    */
  /*     for(iv=0; iv<3; iv++) */
  /* 	{ */
  /* 	  ic=tri2verIN[it].ver_ind[iv]; */
  /* 	  fprintf(stderr,"ver(%d)=(%f,%f,%f); nml=(%f,%f,%f)\n", ic, */
  /* 		  vertex[ ic ].x, vertex[ ic ].y, vertex[ ic ].z, */
  /* 		  nml[ ic ].x, nml[ ic ].y, nml[ ic ].z); */
  /* 	} */
  /*   } */
}
    
/* #include <stdio.h> */
/* #include <math.h> */
/* int main(int ac, char** av) */
/* { */
/*   int N=10; */
  
/*   printf("N=%d\n", N); */
/*   printf("PBC_INDEX(-3,%d) = %d\n", N, PBC_INDEX(-3,N)); */
/*   printf("PBC_INDEX(-2,%d) = %d\n", N, PBC_INDEX(-2,N)); */
/*   printf("PBC_INDEX(-1,%d) = %d\n", N, PBC_INDEX(-1,N)); */
/*   printf("PBC_INDEX(0 ,%d) = %d\n", N, PBC_INDEX( 0,N)); */
/*   printf("PBC_INDEX(1 ,%d) = %d\n", N, PBC_INDEX( 1,N)); */
/*   printf("PBC_INDEX(9 ,%d) = %d\n", N, PBC_INDEX( 9,N)); */
/*   printf("PBC_INDEX(10,%d) = %d\n", N, PBC_INDEX(10,N)); */
/*   printf("PBC_INDEX(11,%d) = %d\n", N, PBC_INDEX(11,N)); */
/*   printf("PBC_INDEX(12,%d) = %d\n", N, PBC_INDEX(12,N)); */
/* } */
