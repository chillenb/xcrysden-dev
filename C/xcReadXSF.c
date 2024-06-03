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
* Source: $XCRYSDEN_TOPDIR/C/xcReadXSF.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <tk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "struct.h"
#include "memory.h"
#include "xcfunc.h"

#define XSF_OPEN   0
#define XSF_UPDATE 1

extern struct DATAGRID *FindDataGrid(int index);
extern FILE *MaybeOpenDataGridFile(char *mode);
extern int GetNumberOfGridBlocks(void);
extern void NewGridList( int gridtype, FILE *gridFP );
extern int ReadBandGrid(FILE *fp, FILE *gridFP, int gridtype, char *ident);

int XC_ReadXSFCmd(ClientData clientData, Tcl_Interp *interp,
		  int argc, const char *argv[]);
int XC_ReadBandXSFCmd(ClientData clientData, Tcl_Interp *interp,
		      int argc, const char *argv[]);

/* ------------------------------------------------------------ *
 *                                                              *
 * xc_readXSF  <XSF-file>                                       *
 *                                                              *
 * ------------------------------------------------------------ */
int 
XC_ReadXSFCmd(ClientData clientData, Tcl_Interp *interp,
	      int argc, const char *argv[])
{
  FILE *fp;

  printf("argc=%d\n", argc);
  fflush(stdout);

  if ( argc != 2 ) {
    Tcl_SetResult(interp, "Usage: xc_readXSF <XSF-file>", TCL_STATIC);
    return TCL_ERROR;
  }

  /* Check if 'file' is OK */
  fp = fopen(argv[1],"r");
  if (fp == NULL) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "Can't open file \"%s\"",argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* command is OK 
   * LET'S READ A FILE
   */
  if (!ReadStructFile(fp, argv[1], FORMAT_XSF, XSF_OPEN)) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"Error reading file \"%s\"",argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    fclose(fp);
    return TCL_ERROR;
  }

  /* assing the sInfo global array */
  Set_sInfoArray( interp );  
  return TCL_OK;
}


/* ------------------------------------------------------------ *
 *                                                              *
 * xc_readbandXSF  <XSF-file>                                   *
 *                                                              *
 * read BANDGRID3d from (band)XSF file                          *
 * ------------------------------------------------------------ */
int 
XC_ReadBandXSFCmd(ClientData clientData, Tcl_Interp *interp,
		  int argc, const char *argv[])
{
  FILE *fp, *gridFP;
  char *line;
  int i, j, c;
  int first_grid_ind, last_grid_ind;
  struct DATAGRID *grid;
  char *string;  
  Tcl_DString *dsPtr = (Tcl_DString *) Tcl_Alloc(sizeof(Tcl_DString));

  if ( argc != 2 ) {
    Tcl_SetResult(interp, "Usage: xc_readbandXSF <(band)XSF-file>", TCL_STATIC);
    return TCL_ERROR;
  }

  /* Check if 'file' is OK */
  fp = fopen(argv[1],"r");
  if (fp == NULL) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "Can't open file \"%s\"",argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /********************* 
   *  command is OK    *
   * LET'S READ A FILE *
   *********************/

  /* maximum line width is 256 character */
  line   = (char *) xcMalloc( sizeof(char) * 256 );
  string = (char *) xcMalloc( sizeof(char) * 10000 );
  first_grid_ind = GetNumberOfGridBlocks();
  while((c = fscanf(fp,"%s\n",line)) != EOF) {

    /* maybe it is a comment line */
    if ( strncmp(line,"#",1) == 0) {
      /* read the rest of the line */
      fscanf(fp,"%[^\n]",line);
      continue;
    }

    /*************************************************************************/
    /* read 3D BANDGRID                                                      */
    /*************************************************************************/
    if ( strncmp(line, "BEGIN_BLOCK_BANDGRID3D", 22) == 0 ||
	 strncmp(line, "BEGIN_BLOCK_BANDGRID_3D", 23) == 0) {
      /* read comment line */
      fprintf(stderr,"BEGIN_BLOCK_BANDGRID_3D was read\n");
      if ( (c = fscanf(fp,"%s\n",line)) == EOF) return XC_ERROR;
      SetDataGridCommentLine( line );      
      gridFP = MaybeOpenDataGridFile("w+");
      NewGridList( DATAGRID_3D, gridFP );      
      /* read DATAGRID_3D_XXX line */
      if ( (c = fscanf(fp,"%s\n",line)) == EOF) return XC_ERROR;
      while (strncmp(line, "END_BLOCK_BANDGRID3D", 15) != 0 &&
	     strncmp(line, "END_BLOCK_BANDGRID_3D", 15) != 0 ) {
	if ( strncmp(line, "BANDGRID_3D", 11) == 0 
	     || strncmp(line, "BEGIN_BANDGRID_3D", 17) == 0 
	     || strncmp(line, "BANDGRID3D", 10) == 0 ) {
	  /*fprintf(stderr,"BANDGRID3D\n");*/
	  if ( ReadBandGrid( fp, gridFP, DATAGRID_3D, line) 
	       == XC_ERROR ) {
	    fprintf(stderr, "ERROR: Error reading BANDGRID_3D_ section, while reading ");
	    free(line);
	    return XC_ERROR;
	  }
	  /* read END_DATAGRID_3D line */
	  /*breakpoint("breakpoint");*/
	  if ( (c = fscanf(fp,"%s\n",line)) == EOF) return XC_ERROR;
	  if (strncmp(line, "END_BANDGRID_3D", 15) != 0 &&
	      strncmp(line, "END_BANDGRID3D", 14) != 0) return XC_ERROR;
	} else return XC_ERROR;
	/* read next line */ 
	if ( (c = fscanf(fp,"%s\n",line)) == EOF) return XC_ERROR;
      }
      xcr.lbandgrid3D = 1;
    }
    else if ( strncmp(line, "INFO_BEGIN", 6) == 0 ||
	      strncmp(line, "BEGIN_INFO", 8) == 0) {
      do {
	c = fscanf(fp,"%s",line);
	if (c == EOF) {
	  fprintf(stderr, "ERROR: unexpected end of file, while reading");
	  free(line);
	  return XC_ERROR;
	}	  
      }	while ( strncmp(line, "INFO_END", 8) != 0 &&
		strncmp(line, "END_INFO", 8) != 0);
    }
    /*************************************************************************/
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss), "expect to read BEGIN_BLOCK_BANDGRID_3D, but got %s, while executing %s %s\n", line, argv[0], argv[1]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }
  last_grid_ind = GetNumberOfGridBlocks();

  /*fprintf(stderr, "first_grid_ind = %d , last_grid_ind = %d\n", first_grid_ind, last_grid_ind);*/

  /*
    return INFO about grid structure; info structure looks like:
    -----------------------------------------------------------
    ngrid
    { grid_index 2D|3D grid_ident grid_nband grid_n_of_subgrids
    {subgrid1_ident subgrid2_ident ...}
    }
  */
  Tcl_DStringInit(dsPtr);
  sprintf(string,"%d ", last_grid_ind - first_grid_ind);
  Tcl_DStringAppendElement(dsPtr, string);
  for (i=first_grid_ind; i<last_grid_ind; i++) 
    {
      /* sublist-level #.1 */
      Tcl_DStringStartSublist(dsPtr);

      grid = FindDataGrid( i );

      sprintf(string,"%d ", grid->index);
      Tcl_DStringAppendElement(dsPtr, string);

      if ( grid->type == DATAGRID_3D ) 
	{
	  sprintf(string,"3D");
	} 
      else 
	{
	  sprintf(string,"2D");
	}

      Tcl_DStringAppendElement(dsPtr, string);
      Tcl_DStringAppendElement(dsPtr, grid->ident);

      sprintf(string,"%d ", grid->nband);
      Tcl_DStringAppendElement(dsPtr, string);    

      sprintf(string,"%d", grid->n_of_subgrids);
      Tcl_DStringAppendElement(dsPtr, string);

      /* sublist-level #.2 */
      Tcl_DStringStartSublist(dsPtr);
      for(j=0; j<grid->n_of_subgrids; j++) 
	{
	  /*fprintf(stderr, "grid #: %d , ident = %s\n", j, grid->subident[j]);*/
	  Tcl_DStringAppendElement(dsPtr, grid->subident[j]);
	}
      Tcl_DStringEndSublist(dsPtr); /* END: sublist-level #.2 */

      Tcl_DStringEndSublist(dsPtr); /* END: sublist-level #.1 */
    }
  Tcl_DStringResult(interp, dsPtr);
  free((FREE_ARG) string);

  return TCL_OK;
}


/* ------------------------------------------------------------ *
 *                                                              *
 * xc_gridvalue  min|max grid_index grid_subindex               *
 *                                                              *
 * returns the minimum|maximum value in the grid                *
 * ------------------------------------------------------------ */
int
XC_GridValueCmd(ClientData clientData, Tcl_Interp *interp,
		int argc, const char *argv[])
{
  int index, subindex;
  struct DATAGRID *grid;
  char *result = (char*) Tcl_Alloc( sizeof(char) * 48);
  if ( argc != 4 ) {
    Tcl_SetResult(interp, "Usage: xc_gridminvalue min|max grid_index griod_subindex", TCL_STATIC);      
    return TCL_ERROR;
  }
  if ( strcmp(argv[1], "min") != 0 && strcmp(argv[1], "max") != 0 ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "wrong value-type \"%s\", must be min or max",
	     argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( Tcl_GetInt(interp, argv[2], &index) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"expected integer, but got %s", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( Tcl_GetInt(interp, argv[3], &subindex) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"expected integer, but got %s", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  grid = FindDataGrid ( index );  
  if ( strcmp(argv[1], "min") == 0 ) {
    sprintf(result, "%f", grid->minvalue[subindex]);
  } else {
    sprintf(result, "%f", grid->maxvalue[subindex]);
  }

  Tcl_SetResult(interp, result, TCL_DYNAMIC);
  return TCL_OK;
}

