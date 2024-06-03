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
* Source: $XCRYSDEN_TOPDIR/C/writeXSF.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <stdio.h>
#include <tcl.h>
#include "struct.h"
#include "molsurf.h"
#include "xcfunc.h"

/* writeXSF.c */
static void _writeMolecule(FILE *fp);
static void _writePeriodic(FILE *fp);
static void _printAtom(FILE *fp, int nat, double x, double y, double z, double *f, boolean printForce);

/*
 * usage: xc_writeXSF filename 
 */
int
XC_WriteXSFCmd(ClientData clientData, Tcl_Interp *interp,
	       int argc, const char *argv[])
{
  FILE *fp;

  if ( argc != 2 ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "wrong # of arguments: should be \"xc_writeXSF filename\"");
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  fp = fopen(argv[1], "w");
  if (!fp) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "couldn't write file %s", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( xcr.dim == 0 ) {
    /* MOLECULE */
    _writeMolecule(fp);
  } else {
    /* PERIODIC STRUCTURE */
    _writePeriodic(fp);
  }

  fclose(fp);
  return TCL_OK;
}



/* ------------------------------------------------------------------------- 
 *                 
 * xc_writebandXSF grid-identifier fermi-energy band-index BXSF-file       
 *                                                              
 * write BANDGRID3d as BXSF file                                
 * ------------------------------------------------------------------------- */

int 
XC_WriteBandXSFCmd(ClientData clientData, Tcl_Interp *interp,
		   int argc, const char *argv[])
{
  register int    i, j, k, ind;
  int             ident, band_ind;
  double          e_f;
  float           value;
  FILE            *fp;
  MOL_SURF        *mols;
  struct DATAGRID *grid;


  if ( argc != 5 ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "Usage: xc_writebandXSF grid-identifier fermi-energy band-index BXSF-file\n");
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }


  /* grid-index */

  if ( Tcl_GetInt(interp, argv[1], &ident) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\"", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  mols = FindMolSurf( ident );
  if (!mols) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"couldn't find MolSurface %d", ident);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  if (mols->isosurf_index < 0) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "couldn't find surface index");
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }


  /* fermi energy */

  if ( Tcl_GetDouble(interp, argv[2], &e_f) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted double, but got \"%s\"", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }



  /* band-index */

  if ( Tcl_GetInt(interp, argv[3], &band_ind) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\"", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }


  /* open file for writing */

  fp = fopen(argv[4],"w");
  if (fp == NULL) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "Can't write file \"%s\"", argv[4]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }


  /* "allocate" the grid */
  grid = FindDataGrid( mols->fs.gridindex );
  fseek(grid->fp, 
	grid->band_fpos[mols->fs.gridsubindex][mols->fs.bandindex], SEEK_SET);


  /* let's write a BXSF file */

  fprintf(fp, " BEGIN_INFO\n");
  fprintf(fp, "   Fermi Energy: %12.6f\n", e_f);
  fprintf(fp, " END_INFO\n");

  fprintf(fp, " BEGIN_BLOCK_BANDGRID_3D\n");
  fprintf(fp, " band_energies\n");
  fprintf(fp, " BANDGRID_3D_BANDS\n");
  /* number of bands */
  fprintf(fp, " 1\n"); 
  /* grid size */
  fprintf(fp, " %d %d %d\n", grid->n[0], grid->n[1], grid->n[2]);
  /* grid origin */
  fprintf(fp, " %15.8f %15.8f %15.8f\n", 
	  grid->orig[0], grid->orig[1], grid->orig[2]); 
  /* grid spanning vectors */
  for (i=0; i<3; i++)
    fprintf(fp, " %15.8f %15.8f %15.8f\n", 
	    grid->vec[i][0], grid->vec[i][1], grid->vec[i][2]);
  fprintf(fp, " BAND: %4d\n", band_ind);

  /* grid values */
  ind = 0;
  for (i=0; i<grid->n[0]; i++)
    for (j=0; j<grid->n[1]; j++)
      for (k=0; k<grid->n[2]; k++) {
	fread(&value, sizeof(float), 1, grid->fp);
	fprintf(fp, " %e", value);
	ind++;
	if (ind == 6) {
	  fprintf(fp,"\n");
	  ind = 0;
	}
      }
  if (ind) fprintf(fp, "\n");

  /* BXSF tail */
  fprintf(fp, " END_BANDGRID_3D\n");
  fprintf(fp, " END_BLOCK_BANDGRID_3D\n");

  /* close the file */
  fclose(fp);

  return TCL_OK;
}



/************************************************************************

static functions below

************************************************************************/

static void _writeMolecule(FILE *fp) {
  register int ia;  

  fprintf(fp, "ATOMS\n");
  for (ia=1; ia<=natoms; ia++) 
    _printAtom(fp, nat[ia], xat[ia]+mx, yat[ia]+my, zat[ia]+mz, fv[ia], xcr.lforce);
}

static void _writePeriodic(FILE *fp) {
  register int ia;  

  /* write periodic structure type */

  if      ( xcr.dim == 3 ) fprintf(fp, "CRYSTAL\n");
  else if ( xcr.dim == 2 ) fprintf(fp, "SLAB\n");
  else                     fprintf(fp, "POLYMER\n");

  /* write PRIMVEC */

  fprintf(fp, "PRIMVEC\n");
  for (ia=0; ia<3; ia++)
    fprintf(fp, "%15.10f  %15.10f  %15.10f\n",
	    vec.prim[ia][0], vec.prim[ia][1], vec.prim[ia][2]);

  /* write CONVVEC */

  if ( xcr.lconvvec ) 
    {
      fprintf(fp, "CONVVEC\n");
      for (ia=0; ia<3; ia++)
	fprintf(fp, "%15.10f  %15.10f  %15.10f\n",
		vec.conv[ia][0], vec.conv[ia][1], vec.conv[ia][2]);
    }          

  /* write PRIMCOORD */

  fprintf(fp, "PRIMCOORD\n");
  fprintf(fp, "%d 1\n", xcr.natr);
  for (ia = 0; ia < xcr.natr; ia++) 
    {
      _printAtom(fp, xcr.prim_nat[ia], xcr.prim_coor[ia][0], xcr.prim_coor[ia][1], xcr.prim_coor[ia][2], 
		 xcr.prim_forc[ia], xcr.prim_lforc);
    }
}

static void _printAtom(FILE *fp, 
		       int nat, double x, double y, double z, double *f, boolean printForce) {
  if ( ! printForce ) {
    fprintf(fp, "%3d    %15.10f  %15.10f  %15.10f\n", nat, x, y, z);
  } else {
    fprintf(fp, "%3d    %15.10f  %15.10f  %15.10f    %15.10f  %15.10f  %15.10f\n", nat, x, y, z, f[0], f[1], f[2]);
  }
}
