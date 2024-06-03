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
* Source: $XCRYSDEN_TOPDIR/C/xcIsoDataGrid.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

/*  #include <togl.h> */
#include <tk.h>
#include <GL/gl.h>
#include <stdlib.h>
#include <stdio.h>
/*  #include <string.h> */
/*  #include <math.h> */
/*  #include "xcGLparam.h" */
/*  #include "primitives.h" */
#include "struct.h"
#include "isosurf.h"
#include "xcfunc.h"
#include "system.h"

extern PLANEVERTEX ***plvertex;
extern ISOSTATE isostate;
/*****************************************************************************
 *                                                                           *
 * xc_isodatagrid <index> <subindex1 factor1> [subindex2 factor2 ...]        *
 *                                                                           *
 *****************************************************************************/
int
XC_IsoDataGridCmd(ClientData clientData, Tcl_Interp *interp,
		  int argc, const char *argv[])
{ 
  FILE   *gridFP;
  struct DATAGRID *grid;
  int    i, j, k, l, io, index, subindex, n_grid_blocks;
  double signfactor;
  float  value;
  char   **argvList = (char **) malloc (2 * sizeof(char *));
  char   *filename  = (char *) malloc (sizeof(char) * 256);

  if (argc % 2 != 0) {
    Tcl_SetResult(interp, "Usage: xc_isodatagrid <index> <subindex1 factor1> [subindex2 factor2 ...]\nor\nxc_isodatagrid info", TCL_STATIC);
    return XC_ERROR;
  }

  if ( strcmp(argv[1], "info" ) == 0 ) {
    /* return the list containing the info about grid structure */
    /* structure of list:

       No_of_grids {#D ident1 No_of_subgrids {subident1 ... subidentNo}
       .....................
       {#D identNo No_of_subgrids {subident1 ... subidentNo}
    */
    struct DATAGRID *grid;
    char *string = (char *) malloc(sizeof(char) * 512);
    Tcl_DString *dsPtr = (Tcl_DString *) Tcl_Alloc(sizeof(Tcl_DString));

    Tcl_DStringInit(dsPtr); /* check if dsPtr should be previously allocated */
    fflush(stdout);
    n_grid_blocks = GetNumberOfGridBlocks();
    sprintf(string,"%d", n_grid_blocks);
    Tcl_DStringAppendElement(dsPtr, string);
    for(i=0; i<n_grid_blocks; i++) {
      Tcl_DStringStartSublist(dsPtr);
      grid = FindDataGrid( i );
      sprintf(string,"2D");
      if ( grid->type == DATAGRID_3D ) sprintf(string,"3D");
      Tcl_DStringAppendElement(dsPtr, string);
      Tcl_DStringAppendElement(dsPtr, grid->ident);
      sprintf(string,"%d", grid->n_of_subgrids);
      Tcl_DStringAppendElement(dsPtr, string);
      Tcl_DStringStartSublist(dsPtr);
      for(j=0; j<grid->n_of_subgrids; j++)
	Tcl_DStringAppendElement(dsPtr, grid->subident[j]);
      Tcl_DStringEndSublist(dsPtr);
      Tcl_DStringEndSublist(dsPtr);
    }
    Tcl_DStringResult(interp, dsPtr);
    free((FREE_ARG) string);
    return TCL_OK;
  }

  if ( Tcl_GetInt(interp, argv[1], &index) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" in xc_isodatagrid command", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* find grid according to the index */
  if ( (grid = FindDataGrid( index )) == NULL ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "grid # %d not found", index);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return XC_ERROR;
  }

  /* initialize grid's signfactor & selected */
  for(i=0; i<DATAGRID_MAXSUBINDEX; i++) {
    grid->signfactor[i] = 0.0;
    grid->selected[i] = 0;
  }

  for(i=2; i<argc; i+=2) {
    if ( Tcl_GetInt(interp, argv[i], &subindex) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" in xc_isodatagrid command", argv[i]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( Tcl_GetDouble(interp, argv[i+1], &signfactor) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted double, but got \"%s\" in xc_isodatagrid command", argv[i+1]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    grid->selected[subindex] = 1;
    grid->signfactor[subindex] = signfactor; 
  }

  /* everything seems to be OK, proceed */  
  argvList[0] = "xc_iso";
  argvList[1] = "init";
  /*nargc = 2;*/
  if (XC_IsoCmd(clientData, interp, 2, (const char **)argvList) == TCL_ERROR) {
    Tcl_SetResult(interp, "an error occurred during initialization of xc_isodatagrid process", TCL_STATIC);
    free((FREE_ARG) argvList);
    return TCL_ERROR;
  }
  free((FREE_ARG) argvList);

  for (i=0; i<3; i++)
    isodata.points[ISOOBJ_BASE][0][i] = grid->orig[i];
  for (i=0; i<grid->type; i++) {
    isodata.vec[ISOOBJ_BASE][i].x = grid->vec[i][0];
    isodata.vec[ISOOBJ_BASE][i].y = grid->vec[i][1];
    isodata.vec[ISOOBJ_BASE][i].z = grid->vec[i][2];
  }  

  if (grid->type == DATAGRID_3D ) {
    float det = XYZ_det3x3(isodata.vec[ISOOBJ_BASE][0], 
			   isodata.vec[ISOOBJ_BASE][1], 
			   isodata.vec[ISOOBJ_BASE][2]);
    if (det > 0) {
      fprintf(stderr,"coordinate system orientation is right\n");
      isodata.cell_orientation = XC_RIGHT;
    } else {
      fprintf(stderr,"coordinate system orientation is left\n");
      isodata.cell_orientation = XC_LEFT;
    }	
  }  

  isostate.stateflag = ISO_INIT | ISO_STACK | ISO_SIGN | ISO_POINTS;
  /* must assing iso_data & iso_files yet */

  grd.nx = newgrd.nx = grid->n[0];
  grd.ny = newgrd.ny = grid->n[1];
  grd.nz = newgrd.nz = grid->n[2];

  /* now read the data and write to bin_vertex_fp */
  gridFP = MaybeOpenDataGridFile("r");

  if (grid->type == DATAGRID_2D) {
    isodata.dim[ISOOBJ_BASE]  =  2;
    plvertex[ISOOBJ_BASE]     =  xcMallocPLANEVERTEX( grd.nx, grd.ny );
    isostate.plvertex_malloc[ISOOBJ_BASE]  =  1;

    for (i=0; i<ISOLINE_MAXLEVEL; i++) {
      isostate.max_n_isoline2D[ISOOBJ_BASE][i] = grd.nx;
      isoline2D[ISOOBJ_BASE].segment[i] = 
	xcMallocLINE( isostate.max_n_isoline2D[ISOOBJ_BASE][i] ); 
    }
    isostate.isoline2D_malloc[ISOOBJ_BASE] = 1;

    for(j=0; j<grid->n[1]; j++)
      for(k=0; k<grid->n[0]; k++)
	plvertex[ISOOBJ_BASE][k][j].val = 0.0;
  }

  if (grid->type == DATAGRID_3D) {
    isodata.dim[ISOOBJ_BASE]   =  3;
    gridvertex                 =  xcMallocGRIDVERTEX( grd.nx, grd.ny, grd.nz );
    isostate.gridvertex_malloc =  1;

    for(j=0; j<grd.nz; j++)
      for(k=0; k<grd.ny; k++)
	for (l=0; l<grd.nx; l++)
	  gridvertex[l][k][j].val = 0.0;
  }

  for (i=0; i<DATAGRID_MAXSUBINDEX; i++)
    if (grid->selected[i]) {
      fseek(gridFP, grid->fpos[i], SEEK_SET);

      if (grid->type == DATAGRID_2D)
	for(j=0; j<grd.ny; j++)
	  for(k=0; k<grd.nx; k++) {
	    fread(&value, sizeof(float), 1, gridFP);
	    plvertex[ISOOBJ_BASE][k][j].val += grid->signfactor[i] * value;
	    plvertex[ISOOBJ_BASE][k][j].p[0] = 
	      ((float) k / (grd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].x +  
	      ((float) j / (grd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].x;
	    plvertex[ISOOBJ_BASE][k][j].p[1] =
	      ((float) k / (grd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].y +  
	      ((float) j / (grd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].y;
	    plvertex[ISOOBJ_BASE][k][j].p[2] =
	      ((float) k / (grd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].z +  
	      ((float) j / (grd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].z;

	    if (j == 0 && k == 0) {
	      isodata.min = plvertex[ISOOBJ_BASE][k][j].val;
	      isodata.max = plvertex[ISOOBJ_BASE][k][j].val;
	    } 
	    else if ( plvertex[ISOOBJ_BASE][k][j].val < isodata.min ) 
	      isodata.min = plvertex[ISOOBJ_BASE][k][j].val;
	    else if ( plvertex[ISOOBJ_BASE][k][j].val > isodata.max ) 
	      isodata.max = plvertex[ISOOBJ_BASE][k][j].val;	    
	  }

      if (grid->type == DATAGRID_3D)
	for (j=0; j<grd.nz; j++)
	  for (k=0; k<grd.ny; k++)
	    for (l=0; l<grd.nx; l++) {
	      fread(&value, sizeof(float), 1, gridFP);
	      gridvertex[l][k][j].val += grid->signfactor[i] * value;
	      gridvertex[l][k][j].p.x =
		((float) l / (grd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].x +  
		((float) k / (grd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].x +  
		((float) j / (grd.nz - 1)) * isodata.vec[ISOOBJ_BASE][2].x;  
	      gridvertex[l][k][j].p.y =   
		((float) l / (grd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].y + 
		((float) k / (grd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].y +
		((float) j / (grd.nz - 1)) * isodata.vec[ISOOBJ_BASE][2].y; 
	      gridvertex[l][k][j].p.z =   
		((float) l / (grd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].z +  
		((float) k / (grd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].z +  
		((float) j / (grd.nz - 1)) * isodata.vec[ISOOBJ_BASE][2].z;   

	      if ( l == 0 && k == 0 && j == 0 ) {
		isodata.min = gridvertex[l][k][j].val;
		isodata.max = gridvertex[l][k][j].val;
	      }
	      if ( gridvertex[l][k][j].val < isodata.min ) 
		isodata.min = gridvertex[l][k][j].val;
	      else if ( gridvertex[l][k][j].val > isodata.max ) 
		isodata.max = gridvertex[l][k][j].val;	    	    	      
	    }
    }
  CloseDataGridFile();

  /* if dim=3D assign also points, vectors & normals for PLANE1, PLANE2, 
     PLANE3 */
  if ( isodata.dim[ISOOBJ_BASE] == 3 ) {
    for(io=ISOOBJ_PLANE1; io<=ISOOBJ_PLANE3; io++)
      for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	  isodata.points[io][i][j] = isodata.points[ISOOBJ_BASE][i][j];

    isodata.vec[ISOOBJ_PLANE1][0] = isodata.vec[ISOOBJ_BASE][0];
    isodata.vec[ISOOBJ_PLANE1][1] = isodata.vec[ISOOBJ_BASE][1];
    isodata.vec[ISOOBJ_PLANE1][2] = isodata.vec[ISOOBJ_BASE][2];

    isodata.vec[ISOOBJ_PLANE2][0] = isodata.vec[ISOOBJ_BASE][2];
    isodata.vec[ISOOBJ_PLANE2][1] = isodata.vec[ISOOBJ_BASE][0];
    isodata.vec[ISOOBJ_PLANE2][2] = isodata.vec[ISOOBJ_BASE][1];

    isodata.vec[ISOOBJ_PLANE3][0] = isodata.vec[ISOOBJ_BASE][1];
    isodata.vec[ISOOBJ_PLANE3][1] = isodata.vec[ISOOBJ_BASE][2];
    isodata.vec[ISOOBJ_PLANE3][2] = isodata.vec[ISOOBJ_BASE][0];

    for(i=ISOOBJ_PLANE1; i<=ISOOBJ_PLANE3; i++) {
      isodata.dim[i] = 2;
      isodata.colnml[i] = 
	VertexNormal(isodata.vec[i][0], isodata.vec[i][1]);
    }
  } else {
    isodata.colnml[ISOOBJ_BASE] = 
      VertexNormal(isodata.vec[ISOOBJ_BASE][0], isodata.vec[ISOOBJ_BASE][1]);
  }

  /* now write data to bin_vertex_fp */
  if (isostate.bin_vertex_file_open) {
    fclose( isodata.bin_vertex_fp );
    isostate.bin_vertex_file_open = 0;
  }
  sprintf(filename, "%s%s%d", 
	  xc_system.scratch_dir, "xc_binVrt.", xc_system.pid);  
  if ( (isodata.bin_vertex_fp = fopen(filename,"w+")) == NULL ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"can't open file \"%s\" while trying to write bin_vertex_fp", filename);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    xcIsoError();
    return TCL_ERROR;
  }
  WriteBinVertexFP();

  /*
    this is a hack to allow just visualization of datagrids -->
    without atom
  */
  if ( natoms == 0 && nframes == 0 ) {
    /* 
       set (mx, my, mz) 
    */
    mx = 2.0 * grid->orig[0] + grid->vec[0][0] + grid->vec[1][0];
    my = 2.0 * grid->orig[1] + grid->vec[0][1] + grid->vec[1][1];
    mz = 2.0 * grid->orig[2] + grid->vec[0][2] + grid->vec[1][2];

    if ( grid->type == DATAGRID_3D ) {
      MVf.structsize = (double)
	DetermineParapipedSize(grid->vec[0], grid->vec[1], 
			       grid->vec[2], grid->orig);
      mx += grid->vec[2][0];
      my += grid->vec[2][1];
      mz += grid->vec[2][2];
    } else {
      MVf.structsize = (double)
	DetermineParalleSize(grid->vec[0], grid->vec[1], grid->orig);
    }

    mx *= 0.5;
    my *= 0.5;
    mz *= 0.5;
  }

  isostate.stateflag |= ISO_DATA | ISO_FILES;
  return TCL_OK;
}
