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
 * Source: $XCRYSDEN_TOPDIR/C/fs.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <togl.h>
#include <stdio.h>
#include <stdlib.h>
#include "struct.h"
#include "bz.h"
#include "xcGLparam.h"
#include "isosurf.h"
#include "molsurf.h"
#include "memory.h"

extern int iroundf(float x);

/* cryClip.c */
extern int ClipSimple(XYZ *p, XYZ n, XYZ p0);
extern int ClipSimple3(XYZ *p, float n[3], float p0[3]);
extern ISOSURFACE *cryIsosurfClip(ISOSURFACE *in, XYZ n, XYZ p0);

/* fs.c */
int fsReadBand(struct DATAGRID *grid, MOL_SURF *mols);
void CropBz(NEW_WIN_CONTEXT *wc, ISOSURFACE *iso);
void CropBz_New(NEW_WIN_CONTEXT *wc, ISOSURFACE *iso);

int
fsReadBand(struct DATAGRID *grid, MOL_SURF *mols)
{
  register int i, j, k, ii, jj, kk, i1, j1, k1;
  float ***bagrid;
  float frmin[3] = {1.0, 1.0, 1.0}, frmax[3] = {-1.0, -1.0, -1.0}, a;
  int N[3], Imin[3], Imax[3];
  float frImin[3], frImax[3];

  bagrid = xcMallocTensor3f(grid->n[0], grid->n[1], grid->n[2]);
  
  fseek(grid->fp, 
	grid->band_fpos[ mols->fs.gridsubindex ][ mols->fs.bandindex ], 
	SEEK_SET);

  for (i=0; i<grid->n[0]; i++)
    for (j=0; j<grid->n[1]; j++)
      for (k=0; k<grid->n[2]; k++)
	fread(&(bagrid[i][j][k]), sizeof(float), 1, grid->fp);

  for (k=0; k<3; k++)
    N[k]=grid->n[k]-1;

  fprintf(stderr,"FS:----------------\n");
  fprintf(stderr,"FS: N == %d, %d, %d\n", N[0], N[1], N[2]);

  if ( mols->fs.celltype == XCR_BZ ) {
    /**********/
    /* XCR_BZ */
    /**********/
    /* 
       query for BZ size then make a box just big enough to hold the 
       whole BZ and then reorder the grid so that it will match the BZ;
       use E(k) = E(K+k) relation; 
       NO:  use vectors from bandgrid not from RECIP-!!! 
       YES: use vectors from RECIP-!!!
    */
    for (i=0; i<bz[BZ_PRIMCELL].npoly; i++) 
      for (j=0; j<bz[BZ_PRIMCELL].nvert[i]; j++) 
	for (k=0; k<3; k++) {
	  /* use transposed direct primitive vectors */
	  a = 
	    bz[BZ_PRIMCELL].poly[i][j][0]*vec.prim[k][0] + 
	    bz[BZ_PRIMCELL].poly[i][j][1]*vec.prim[k][1] + 
	    bz[BZ_PRIMCELL].poly[i][j][2]*vec.prim[k][2];
	  if ( a < frmin[k] ) frmin[k] = a;
	  if ( a > frmax[k] ) frmax[k] = a;
	}

    fprintf(stderr,"FS: frmin == %6.3f, %6.3f, %6.3f\n", frmin[0], frmin[1], frmin[2]);
    fprintf(stderr,"FS: frmax == %6.3f, %6.3f, %6.3f\n", frmax[0], frmax[1], frmax[2]);
    
    for (k=0; k<3; k++) {
      Imin[k] = (int) (frmin[k] * (float) N[k]) - 1;
      Imax[k] = (int) (frmax[k] * (float) N[k]) + 1;
      fprintf(stderr,"FS: min[%d] = %6.3f   ;  max[%d] = %6.3f\n", 
	      k, frmin[k] * (float) N[k], k, frmax[k] * (float) N[k]);

      /* prevent  crash due to roundoff error */
      if (Imin[k] < -N[k]) Imin[k] = -N[k];
      if (Imax[k] >  N[k]) Imax[k] =  N[k];

      /* 
	 WARNING: the Imin/Imax must be -+N[k] as long as we don't have a non-periodic LCASI
	 interpolation. Setting Imin/Imax to -+N[k] assures the grid is periodic, otherwise with
	 the really minimal Imin/Imax it is not !!!
       */
      Imin[k] = -N[k];
      Imax[k] = +N[k];

      frImin[k] = (float) Imin[k] / N[k];
      frImax[k] = (float) Imax[k] / N[k];
    }

    fprintf(stderr,"FS: Imin == %d, %d, %d\n", Imin[0], Imin[1], Imin[2]);
    fprintf(stderr,"FS: Imax == %d, %d, %d\n", Imax[0], Imax[1], Imax[2]);

    fprintf(stderr,"FS: frImin == %f, %f, %f\n", frImin[0], frImin[1], frImin[2]);
    fprintf(stderr,"FS: frImax == %f, %f, %f\n", frImax[0], frImax[1], frImax[2]);

    mols->i = Imax[0]-Imin[0]+1;
    mols->j = Imax[1]-Imin[1]+1;
    mols->k = Imax[2]-Imin[2]+1;

    fprintf(stderr,"FS: mols(i,j,k) == %d %d %d\n", mols->i, mols->j, mols->k);

    /* malloc mols->molgrid space */
    if ( mols->molgrid ) xcFree_Tensor3f( mols->molgrid );
    mols->molgrid = xcMallocTensor3f( mols->i, mols->j, mols->k );
				
    for (k=0; k<3; k++) {
      /* this is wrong (too much content of the grid will be forced to be squeezed inside the BZ) */
      /*
      mols->lowcoor[k] = 
	frmin[0]*vec.recprim[0][k] +
	frmin[1]*vec.recprim[1][k] +
	frmin[2]*vec.recprim[2][k];
      */

      /* this is correct */
      mols->lowcoor[k] = 
	frImin[0]*vec.recprim[0][k] +
	frImin[1]*vec.recprim[1][k] +
	frImin[2]*vec.recprim[2][k];

      for (j=0; j<3; j++) {
	/* this is wrong (too much content of the grid will be forced to be squeezed inside the BZ) */
	/*
	mols->vec[k][j] = (frmax[k]-frmin[k]) * vec.recprim[k][j];
	*/

	/*  this is correct */
	mols->vec[k][j] = (frImax[k]-frImin[k]) * vec.recprim[k][j];
	mols->isoexpand.rep_vec[k][j] = vec.recprim[k][j];
      }
    }
    fprintf(stderr,"FS: mols->lowcoor == %f %f %f\n", 
	    mols->lowcoor[0], mols->lowcoor[1], mols->lowcoor[2]);
    /*fprintf(stderr,"FS: orig == %f %f %f\n", orig[0], orig[1], orig[2]);*/
	    
    for (i=Imin[0], i1=0; i<=Imax[0]; i++, i1++) {
      ii=i;
      if (i<0) ii=N[0]+i;
      for (j=Imin[1], j1=0; j<=Imax[1]; j++, j1++) {
	jj=j;
	if (j<0) jj=N[1]+j;
	for (k=Imin[2], k1=0; k<=Imax[2]; k++, k1++) {
	  kk=k;
	  if (k<0) kk=N[2]+k;	  
	  mols->molgrid[i1][j1][k1] = bagrid[ii][jj][kk];	  
	  /*fprintf(stderr,"FS: BZgrid(%d,%d,%d) <-- CELLgrid(%d,%d,%d)\n", i1, j1, k1,  ii, jj, kk);*/
	}
      }
    }
    xcFree_Tensor3f( bagrid );
  } else {
    /*******************/
    /* XCR_PARAPIPEDAL */
    /*******************/
    mols->i = grid->n[0];
    mols->j = grid->n[1];
    mols->k = grid->n[2];

    fprintf(stderr,"FS: mols(i,j,k) == %d %d %d\n", mols->i, mols->j, mols->k);

    if ( mols->molgrid ) xcFree_Tensor3f( mols->molgrid );    
    mols->molgrid = bagrid;

    for (k=0; k<3; k++) {
      mols->lowcoor[k] = 0.0;
      for (j=0; j<3; j++) {
	mols->vec[k][j] = (float) vec.recprim[k][j];
	mols->isoexpand.rep_vec[k][j] = vec.recprim[k][j];
      }
    }
  }
  return XC_OK;
}

void
CropBz( NEW_WIN_CONTEXT *wc, ISOSURFACE *iso ) 
{
  register int i, ip, j3, index;
  XYZ tri[3];

  for (i=0; i<iso->ntriangl; i++) {
    for (j3=0; j3<3; j3++) {
      index = iso->tri2verIN[i].ver_ind[j3];
      tri[j3] = iso->vertex[index];
    }      
    for (ip=0; ip < bz[BZ_PRIMCELL].npoly; ip++) {
      iso->triangl_status[i] = ClipSimple3( tri, 
					    bz[BZ_PRIMCELL].norm[ip][0],
					    bz[BZ_PRIMCELL].poly[ip][0]);
      if (!iso->triangl_status[i]) break;
    }
  }
}

void
CropBz_New( NEW_WIN_CONTEXT *wc, ISOSURFACE *iso ) 
{
  register int ip;
  XYZ          pos, nml;
  ISOSURFACE   *auxil1, *auxil2;

  /*  allocate out ISOSURFACE and initialize it !!! */
  auxil1 =  (ISOSURFACE *) calloc (1, sizeof(ISOSURFACE));
  
  auxil1->index           = iso->index;				   
  auxil1->revertnormals   = iso->revertnormals;
  auxil1->nvertex         = iso->nvertex;			   
  auxil1->ntriangl        = iso->ntriangl;			   
  auxil1->smooth_nstep    = iso->smooth_nstep;			   
  auxil1->smooth_weight   = iso->smooth_weight;			   
  auxil1->triangl_status  = iso->triangl_status;			   
  auxil1->vertex	  = iso->vertex;				   
  auxil1->vertex_orig     = iso->vertex_orig; 			   
  auxil1->normal	  = iso->normal;				   
  auxil1->color	          = iso->color;			   
  auxil1->tri2verIN	  = iso->tri2verIN;			   
  auxil1->nver2triIN	  = iso->nver2triIN;			   
  auxil1->ver2triIN	  = iso->ver2triIN;
  
  fprintf(stderr,"Number of IN triangles:  %d\n", iso->ntriangl);
  fprintf(stderr,"Number of IN vertices:   %d\n", iso->nvertex);  

  for (ip=0; ip < bz[BZ_PRIMCELL].npoly; ip++) {
    pos.x = bz[BZ_PRIMCELL].poly[ip][0][0];
    pos.y = bz[BZ_PRIMCELL].poly[ip][0][1];
    pos.z = bz[BZ_PRIMCELL].poly[ip][0][2];

    nml.x = bz[BZ_PRIMCELL].norm[ip][0][0];
    nml.y = bz[BZ_PRIMCELL].norm[ip][0][1];
    nml.z = bz[BZ_PRIMCELL].norm[ip][0][2];

    auxil2 = cryIsosurfClip(auxil1, nml, pos);

    xcFree ((void*) auxil1->triangl_status);
    xcFree ((void*) auxil1->vertex);
    xcFree ((void*) auxil1->vertex_orig); 
    xcFree ((void*) auxil1->normal);
    xcFree ((void*) auxil1->color);
    xcFree ((void*) auxil1->tri2verIN);
    xcFree ((void*) auxil1->nver2triIN);
    xcFree ((void*) auxil1->ver2triIN);
    xcFree ((void*) auxil1);
    auxil1 = auxil2;
  }
  
  /* copy isosurface and free unnecessary space */

  /* free .... */
  /*
  xcFree ((void*) iso->triangl_status);
  xcFree ((void*) iso->vertex);
  xcFree ((void*) iso->vertex_orig); 
  xcFree ((void*) iso->normal);
  xcFree ((void*) iso->color);
  xcFree ((void*) iso->tri2verIN);
  xcFree ((void*) iso->nver2triIN);
  xcFree ((void*) iso->ver2triIN);
  */

  /* copy */
  iso->index           = auxil1->index;				   
  iso->revertnormals   = auxil1->revertnormals;
  iso->nvertex         = auxil1->nvertex;			   
  iso->ntriangl        = auxil1->ntriangl;			   
  iso->smooth_nstep    = auxil1->smooth_nstep;			   
  iso->smooth_weight   = auxil1->smooth_weight;			   
  iso->triangl_status  = auxil1->triangl_status;			   
  iso->vertex	       = auxil1->vertex;				   
  iso->vertex_orig     = auxil1->vertex_orig; 			   
  iso->normal	       = auxil1->normal;				   
  iso->color	       = auxil1->color;			   
  iso->tri2verIN       = auxil1->tri2verIN;			   
  iso->nver2triIN      = auxil1->nver2triIN;			   
  iso->ver2triIN       = auxil1->ver2triIN;

  fprintf(stderr,"Number of OUT triangles: %d\n",   iso->ntriangl);
  fprintf(stderr,"Number of OUT vertices:  %d\n\n", iso->nvertex);
  fflush(stderr);
}

