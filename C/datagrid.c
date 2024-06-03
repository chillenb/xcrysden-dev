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
 * Source: $XCRYSDEN_TOPDIR/C/datagrid.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "struct.h"
#include "isosurf.h"
#include "system.h"


extern PLANEVERTEX ***plvertex;

struct DATAGRID *grid, *gridPtr = NULL;

static int datagrid_index = 0;
static int datagrid_subindex = 0;
static char *comment_line;
static int comment_line_malloced = 0;

struct DATAGRID_FILE {
  int opened;
  FILE *fp;
  char *name;
  int name_malloced;
} dgf = {0, NULL, NULL, 0};

/* --- function prototypes --- */
static void AddToGridList(struct DATAGRID *g);
void CloseGridList(void);
struct DATAGRID *FindDataGrid(int index);
int GetNumberOfGridBlocks(void);
void NewGridList( int gridtype, FILE *gridFP );
FILE *MaybeOpenDataGridFile(char *mode);
void CloseDataGridFile(void);
void SetDataGridCommentLine(char *line);
int ReadDataGrid(FILE *fp, FILE *gridFP, int gridtype, char *ident);
int ReadBandGrid(FILE *fp, FILE *gridFP, int gridtype, char *ident);
int WriteDataGrid(FILE *fp, int gridtype, const char *ident, int obj);

static void AddToGridList(struct DATAGRID *g)
{
  g->ptr  = gridPtr;
  gridPtr = g;
}

void CloseGridList(void)
{
  struct DATAGRID *g = gridPtr, *f;
  
  while(g) {
    f = g->ptr;
    free((FREE_ARG) g);
    g = f;
  }
  datagrid_index = 0;
  datagrid_subindex = 0;
  comment_line_malloced = 0;
  CloseDataGridFile();
  dgf.opened = 0;
  dgf.fp     = NULL;
  dgf.name   = NULL;
  dgf.name_malloced = 0;

  gridPtr = NULL;
}

struct DATAGRID *FindDataGrid(int index)
{
   struct DATAGRID *g = gridPtr;
   while (g) {
     if (index == g->index)
       return g;
     g = g->ptr;
   }
   return NULL;
}


int
GetNumberOfGridBlocks(void)
{
  return datagrid_index;
}

void
NewGridList( int gridtype, FILE *gridFP ) 
{
  grid = (struct DATAGRID *) malloc ( sizeof(struct DATAGRID) );
  grid->fp     = gridFP;
  grid->type   = gridtype;
  grid->index  = datagrid_index++;
  grid->ident  = (char *) malloc( sizeof(char) * 512 );
  sprintf(grid->ident, "%s", comment_line);
  datagrid_subindex = 0;
  AddToGridList( grid );
}


FILE *MaybeOpenDataGridFile(char *mode) 
{  
  if (!dgf.name_malloced) {
    dgf.name = (char *) malloc ( sizeof(char) * 256 );
    sprintf(dgf.name, "%s%s%d", 
	    xc_system.scratch_dir, "xc_datagrid.", xc_system.pid);
    dgf.name_malloced = 1;
  }
  if (!dgf.opened) {
    dgf.fp = fopen(dgf.name, mode);
    dgf.opened = 1;
  }

  return dgf.fp;
}
    

void
CloseDataGridFile(void) {
  if (dgf.opened) {
    fclose(dgf.fp);
    dgf.opened = 0;
  }
}
 

void
SetDataGridCommentLine(char *line)
{
  if (!comment_line_malloced) {
    comment_line = (char *) malloc( sizeof(char) * 512 );
    comment_line_malloced = 1;
  }
  sprintf(comment_line, "%s", line);
}


int
ReadDataGrid(FILE *fp, FILE *gridFP, int gridtype, char *ident) 
{
  register int i, j, k;
  int len;
  float value;
  int n[3] = {0, 0, 0};
  float o[3] = { 0.0, 0.0, 0.0 };
  float v[3][3] = {
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0}};

  grid->fpos[datagrid_subindex] = ftell( gridFP );
  len = strlen(ident) + 10;
  grid->subident[datagrid_subindex] = (char *) calloc( len, sizeof(char) );
  sprintf(grid->subident[datagrid_subindex],"%s", ident);

  if ( gridtype == DATAGRID_2D ) {
    if ( fscanf(fp,"%d %d", n, n+1) != gridtype )
      return XC_ERROR;  
  }
  else if ( gridtype == DATAGRID_3D ) {
    if ( fscanf(fp,"%d %d %d", n, n+1, n+2) != gridtype )
      return XC_ERROR;
  }
    
  
  if ( fscanf(fp,"%f %f %f", o, o+1, o+2) != 3 ) return XC_ERROR;
  
  for (i=0; i<gridtype; i++)
    if ( fscanf(fp,"%f %f %f", v[i], v[i]+1, v[i]+2) != 3 ) return XC_ERROR;
  
  /* check if DATAGRID header match previous one */
  if (datagrid_subindex == 0) {
    for(i=0; i<3; i++) {
      grid->n[i] = n[i];
      grid->orig[i] = o[i];
      for(j=0; j<3; j++)
	grid->vec[i][j] = v[i][j];
    }
  } else {
    for(i=0; i<3; i++) {
      if ( n[i] != grid->n[i] || o[i] != grid->orig[i] ) return XC_ERROR;
      for(j=0; j<3; j++)
	if ( v[i][j] != grid->vec[i][j] ) return XC_ERROR;
    }
  }

  if ( gridtype == DATAGRID_2D) {
    for (i=0; i<grid->n[1]; i++)
      for (j=0; j<grid->n[0]; j++) {
	if ( fscanf(fp,"%f", &value) != 1 ) return XC_ERROR;
	fwrite(&value, sizeof(float), 1, gridFP);
      }
  } else if ( gridtype == DATAGRID_3D ) {
    len=0;
    for (i=0; i<grid->n[2]; i++)
      for (j=0; j<grid->n[1]; j++)
	for (k=0; k<grid->n[0]; k++) {
	  if ( fscanf(fp,"%f", &value) != 1 ) {
	    fprintf(stderr,"(error in reading datagrid values; only %d values read)\n", len);
	    return XC_ERROR;
	  }
	  fwrite(&value, sizeof(float), 1, gridFP);
	  len++;
	}
  }

  grid->n_of_subgrids = ++datagrid_subindex;

  return XC_OK;
}



int
ReadBandGrid(FILE *fp, FILE *gridFP, int gridtype, char *ident) 
{
  register int i, j, k, ib;
  int len, nband;
  float value;
  int n[3] = {0, 0, 0};
  float o[3] = { 0.0, 0.0, 0.0 };
  char band[5];
  float v[3][3] = {
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0}};

  grid->fpos[datagrid_subindex] = ftell( gridFP );
  len = strlen(ident) + 10;
  grid->subident[datagrid_subindex] = (char *) calloc( len, sizeof(char) );
  sprintf(grid->subident[datagrid_subindex],"%s", ident);

  fscanf(fp,"%d", &nband);
  grid->nband = nband;
  grid->band_index[datagrid_subindex] = 
    (int *) malloc( (grid->nband+1) * sizeof(int) );
  grid->band_fpos[datagrid_subindex] = 
    (long *) malloc( (grid->nband+1) * sizeof(long) );

  grid->lband = 1;

  if ( gridtype == DATAGRID_2D ) {
    if ( fscanf(fp,"%d %d", n, n+1) != gridtype )
      return XC_ERROR;  
  }
  else if ( gridtype == DATAGRID_3D ) {
    if ( fscanf(fp,"%d %d %d", n, n+1, n+2) != gridtype )
      return XC_ERROR;
  }

  /* read ORIGIN and VECTORS */
  if ( fscanf(fp,"%f %f %f", o, o+1, o+2) != 3 ) return XC_ERROR;
  for (i=0; i<gridtype; i++)
    if ( fscanf(fp,"%f %f %f", v[i], v[i]+1, v[i]+2) != 3 ) return XC_ERROR;

  /* check if DATAGRID header match previous one */
  if (datagrid_subindex == 0) {
    for(i=0; i<3; i++) {
      grid->n[i] = n[i];
      grid->orig[i] = o[i];
      for(j=0; j<3; j++)
	grid->vec[i][j] = v[i][j];
    }
  } else {
    for(i=0; i<3; i++) {
      if ( n[i] != grid->n[i] || o[i] != grid->orig[i] ) return XC_ERROR;
      for(j=0; j<3; j++)
	if ( v[i][j] != grid->vec[i][j] ) return XC_ERROR;
    }
  }
  
  for(ib=0; ib<grid->nband; ib++) {
    /* read BAND: */
    if ( fscanf(fp, "%5s %i", band, 
		&grid->band_index[datagrid_subindex][ib]) != 2 ) {
      fprintf(stderr, "ERROR: Error reading %d-th BAND !!!\n", ib);
      return XC_ERROR;
    }
    grid->band_fpos[datagrid_subindex][ib] = ftell( gridFP );

    grid->minvalue[ib] = +1.0e+10;
    grid->maxvalue[ib] = -1.0e+10;
   
    if ( gridtype == DATAGRID_2D) {
      for (i=0; i<grid->n[0]; i++)
	for (j=0; j<grid->n[1]; j++) {
	  if ( fscanf(fp,"%f", &value) != 1 ) return XC_ERROR;
	  fwrite(&value, sizeof(float), 1, gridFP);
	  if ( value < grid->minvalue[ib] ) 
	    grid->minvalue[ib] = value;
	  if ( value > grid->maxvalue[ib] ) 
	    grid->maxvalue[ib] = value;
	}
    } else if ( gridtype == DATAGRID_3D ) {
      for (i=0; i<grid->n[0]; i++)
	for (j=0; j<grid->n[1]; j++)
	  for (k=0; k<grid->n[2]; k++) {
	    if ( fscanf(fp,"%f", &value) != 1 ) return XC_ERROR;
	    fwrite(&value, sizeof(float), 1, gridFP);
	    if ( value < grid->minvalue[ib] ) 
	      grid->minvalue[ib] = value;
	    if ( value > grid->maxvalue[ib] ) 
	      grid->maxvalue[ib] = value;
	  }
    }
  }

  grid->n_of_subgrids = ++datagrid_subindex;

  return XC_OK;
}


int
WriteDataGrid(FILE *fp, int gridtype, const char *ident, int obj)
{
  char *keyword = (char *) malloc( sizeof(char) * 256 );
  char *key;
  register int i, j, k, ind;

  key = "DATAGRID_2D_";
  if (gridtype == DATAGRID_3D) key = "DATAGRID_3D_";
  sprintf(keyword, "%s%s", key, ident);

  fprintf(fp, "%s\n", keyword);
  if (gridtype == DATAGRID_2D) {
    int n, m;
    n = newgrd.nx;
    m = newgrd.ny;
    if ( obj == ISOOBJ_PLANE2 ) {
      n = newgrd.nz;
      m = newgrd.nx;
    }
    if ( obj == ISOOBJ_PLANE3 ) {
      n = newgrd.ny;
      m = newgrd.nz;
    }    
    fprintf(fp, "%d %d\n", n, m);
    fprintf(fp, "%e %e %e\n", isodata.points[obj][0][0],
	   isodata.points[obj][0][1], isodata.points[obj][0][2]);
    for (i=0; i<2; i++)
      fprintf(fp, "%e %e %e\n", isodata.vec[obj][i].x,
	     isodata.vec[obj][i].y, isodata.vec[obj][i].z);
    ind = 0;
    for (j=0; j<m; j++)
      for (k=0; k<n; k++) {
	fprintf(fp, " %e", plvertex[obj][k][j].val);
	ind++;
	if (ind == 6) {
	  fprintf(fp,"\n");
	  ind = 0;
	}
      }	
  } else if (gridtype == DATAGRID_3D) {
    int nx, ny, nz;
    nx = newgrd.nx;
    ny = newgrd.ny;
    nz = newgrd.nz;
    fprintf(fp, "%d %d %d\n", nx, ny, nz);
    fprintf(fp, "%e %e %e\n", isodata.points[obj][0][0],
	   isodata.points[obj][0][1], isodata.points[obj][0][2]);
    for (i=0; i<3; i++)
      fprintf(fp, "%e %e %e\n", isodata.vec[obj][i].x,
	     isodata.vec[obj][i].y, isodata.vec[obj][i].z);
    ind = 0;
    for (i=0; i<nz; i++)
      for (j=0; j<ny; j++)
	for (k=0; k<nx; k++) {
	  fprintf(fp, " %e", gridvertex[k][j][i].val);
	  ind++;
	  if (ind == 6) {
	    fprintf(fp,"\n");
	    ind = 0;
	  }
	}	
  }
  if (ind) fprintf(fp, "\n");
  if (gridtype == DATAGRID_2D) fprintf(fp, "END_DATAGRID_2D");
  if (gridtype == DATAGRID_3D) fprintf(fp, "END_DATAGRID_3D");

  return XC_OK;
}
