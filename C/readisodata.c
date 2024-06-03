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
 * Source: $XCRYSDEN_TOPDIR/C/readisodata.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <stdio.h>
#include <tcl.h>
#include <stdlib.h>
#include <ctype.h>
#include "struct.h"
#include "isosurf.h"

extern ISOSTATE isostate;
extern PLANEVERTEX ***plvertex;

/* readisodata.c */
int ReadIsoData(int framenum[2 * 100 ][4 ][2], int ntimes);
void WriteBinVertexFP(void);
int ReadBlock0(int type, int s0, int s1, int s2, int s3);

/* -----------------------------------------------------------------------
 * WARNING:: this routine is written up to  
 *           isodata.maxstack = 3. If more stacks will be 
 *           needed, you will have to modified this !!!!!!!!!!!!
 * ----------------------------------------------------------------------*/
int
ReadIsoData(int framenum[2 * ISODATA_MAXFRAME][ISODATA_MAXSTACK][2], 
	    int ntimes) 
{
  register int i, j, jj, k, s3, s2, s1, s0;
  int nframe = 0, nframe3[ISODATA_MAXSTACK];
  int n0[4], n1[4];
  int is_present;  
  int nblock = 0, nblock3, tmpblock;
  int last_read;
  int start_read = 0;  /* where to start reading */
  struct index {
    int block[ISODATA_MAXBLOCK];
    int s0[ISODATA_MAXBLOCK];
    int s1[ISODATA_MAXBLOCK];
    int s2[ISODATA_MAXBLOCK];
    int s3[ISODATA_MAXBLOCK];
  } read;

  /* -----------------------------------------------------------------------
   * framenum[][j][]
   * j=0    --> stack #3
   * j=1    --> stack #2
   * j=2    --> stack #1
   * j=3    --> stack #0
   * -----------------------------------------------------------------------*/

  /* first calculate the number of blocks;
   * block is sometning like frame in stack-order 0 
   */  
  for (i=0; i<isodata.nframe[0][3]-1; i++) {
    nframe3[i] = isodata.nframe[i][0] * isodata.nframe[i][1] *
      isodata.nframe[i][2];
    nframe    += nframe3[i];
  }

  /* now calculate which of blocks to read */
  for (i=0; i<ntimes; i++) {
    /* single frame in stack3 or range of frames ????? */
    for(j=0; j<=3; j++) {
      jj = 3 - j; /* framenum[][0][] stands for 3-th stack-order, 
                   * so this is needed
                   */
      n0[jj] = framenum[i][j][0];
      if ( n0[jj] < 0 ) {
	/* all n0[jj] must be defined (not to be -1 -> 
         * this is considered as an error )
         */
	sprintf(isodata.error, "some stacks were not properly initialized with \"xc_isodata\" command !!!");
	return XC_ERROR;
      }
      n1[jj] = framenum[i][j][1];
      if ( n1[jj] <= framenum[i][j][0] ) 
	n1[jj] = framenum[i][j][0];      
    }
    
    for (s3=n0[3]; s3 <= n1[3]; s3++)
      for (s2=n0[2]; s2 <= n1[2]; s2++)
	for (s1=n0[1]; s1 <= n1[1]; s1++)
	  for (s0=n0[0]; s0 <= n1[0]; s0++) {	
	    nblock3 = 0;
	    for (j=0; j<s3; j++)
	      nblock3 += nframe3[j];
	    tmpblock = nblock3;
	    if (s2 > 0) 
	      tmpblock += s2 * isodata.nframe[s3][1] * isodata.nframe[s3][0];
	    if (s1 > 0) 
	      tmpblock += s1 * isodata.nframe[s3][0];
	    tmpblock += s0;
	    if ( tmpblock < start_read ) {
	      sprintf(isodata.error, "xc_isodata command was not used properly, because frames are not specified in sequential fashion");
	      return XC_ERROR;
	    }
	    read.block[nblock] = tmpblock;
	    read.s0[nblock] = s0;
	    read.s1[nblock] = s1;
	    read.s2[nblock] = s2;
	    read.s3[nblock] = s3;	    
	    nblock++;
	  }

    last_read = 0;
    for (j=start_read; j<=read.block[nblock - 1]; j++) {
      is_present = RB0_DUMMY;
      for (k=last_read; k<nblock; k++) {
	if ( read.block[k] == j ) {
	  is_present = RB0_READ;
	  last_read = k;
	  break;
	}
      }
      if (!ReadBlock0(is_present, read.s0[k], read.s1[k], read.s2[k], read.s3[k])) 
	return XC_ERROR;
    }
    start_read = read.block[nblock - 1] + 1;
  }
  
  /* 
   * assign isodata.min_allowed[] | isodata.max_allowed
   */
  for (i=ISOOBJ_BASE; i< MAX_ISOOBJECTS; i++) {
    isodata.min_allowed[i] = isodata.min;
    isodata.max_allowed[i] = isodata.max;
  }

  /* 
     now I should write the bin_vertex_fp file
  */
  WriteBinVertexFP();

  return XC_OK;
}



void WriteBinVertexFP(void)
{
  register int i, j, k;
  if ( isodata.dim[ISOOBJ_BASE] == 2 ) {
    for (j=0; j<grd.ny; j++)
      for (k=0; k<grd.nx; k++)
	fwrite(&(plvertex[ISOOBJ_BASE][k][j].val), sizeof(float), 1, 
	       isodata.bin_vertex_fp);
  } else if ( isodata.dim[ISOOBJ_BASE] == 3 ) {
    for (i=0; i<grd.nz; i++)
      for (j=0; j<grd.ny; j++)
	for (k=0; k<grd.nx; k++)    
	  fwrite(&(gridvertex[k][j][i].val), sizeof(float), 1, 
		 isodata.bin_vertex_fp);
  }
}



int 
ReadBlock0(int type, int s0, int s1, int s2, int s3)
{
  /*  static FILE *fff;
  static int ifff = 1;*/

  /* type = -1; initialization
   * type = 0; dummy reading
   * type = 1; read data
   * type = 2; found out the number of frames in stack #0
   */
  register int i, j, k, n;  
  float value;
  float x0, y0, dx, dy, cosxy;
  float xa, ya, za, xb, yb, zb, xc, yc, zc;
  int   grdnx, grdny;
  static int first_time = 1;  /* used for grd structure */
  static int first_t = 1;     /* used for isodata struc */
  /* this is TEMP */
  /* FILE *file = fopen("bin.dat","a"); */
  
  /*  if ( ifff )  {
    fff = fopen("grid.dat","w"); 
    ifff = 0; 
  }*/

  if ( type == RB0_INIT ) {    
    /* if 2D then we will need "plvertex" structure, but if 3D then we will 
     * need "gridvertex" structure; reallocate memory for needed structure */
    if ( isodata.dim[ISOOBJ_BASE] == 2 ) {
      for (i=0; i<grd.nx; i++)
	for (j=0; j<grd.ny; j++) 	  
	  plvertex[ISOOBJ_BASE][i][j].val = 0.0;
    } 
    else if ( isodata.dim[ISOOBJ_BASE] == 3 ) {
      for (i=0; i<grd.nx; i++)
	for (j=0; j<grd.ny; j++) 
	  for(k=0; k<grd.nz; k++) 
	    gridvertex[i][j][k].val = 0.0;
    }
    first_t    = 1;
    first_time = 1;

    /* if binary file is opened, rewind it */
    if ( isostate.bin_file_open ) 
      rewind( isodata.bin_fp );
  }
  /* type == 2 ==> findout the number of frames in stack 0 
   * -----------------------------------------------------
   * first two lines are of form:
   * CRYSTAL95:
   *  -%- MAPN    0  0  0  0  0  0            **********               
   *  9  13  0.000000  0.000000  0.472432  0.472432  0.000000
   * ...
   * CRYSTAL98:
   * -%-2MAPN  10  10 0.60279E+00 0.60279E+00 0.50000E+00
   * 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 3.83615E+00 3.83615E+00
   * 3.83615E+00 0.00000E+00 3.83615E+00
   * 7.17593E+04 3.69350E+00 7.58363E-01 1.12803E-01 3.65208E-02 3.65208E-02
   * 1.12803E-01 7.58363E-01 3.69350E+00 7.17593E+04 3.69350E+00 3.69350E+00
   * ...
   *
   * CRYSTAL03:
   * -%-2MAPN  10  10 0.60279E+00 0.60279E+00 0.50000E+00
   * 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 3.83615E+00 3.83615E+00
   * 3.83615E+00 0.00000E+00 3.83615E+00 n_atoms dimension
   * 7.17593E+04 3.69350E+00 7.58363E-01 1.12803E-01 3.65208E-02 3.65208E-02
   * 1.12803E-01 7.58363E-01 3.69350E+00 7.17593E+04 3.69350E+00 3.69350E+00
   * ...
   * foreach atom: nat symb x y z
   * lattice_vec1
   * lattice_vec2
   * lattice_vec3
   *
   * CRYSTAL06:
   * -%-0MAPN   28   28 0.20611E+00 0.20611E+00 0.00000E+00
   * 0.00000E+00 0.00000E+00 0.00000E+00 5.56500E+00 0.00000E+00 0.00000E+00
   * 5.56500E+00 5.56500E+00 0.00000E+00      16   3
   */  
  else if ( type == RB0_FIND ) {
    int c = 0, nn = 0;
    int natoms, dim, nat, ia;
    float x_, y_, z_, lat_;
    static int oldnx, oldny;
    char word1[8], word2[8]; 
    while ( c != EOF ) {
      /* read first line */
      /* sprintf(word1,"----",NULL); */
      /* sprintf(word2,"----",NULL); */
      if ( fscanf(isodata.fp,"%4c %4c",word1,word2) == EOF ) break;
      /* check if it is CRYSTAL98|03 ??? */
      if ( isdigit(word1[3]) ) {

	fprintf(stderr,"crystal_version == %d\n", crystal_version);

	if ( (crystal_version == 95) || (crystal_version == 98) || (crystal_version == 3) ) 
	  {
	    if ( (n = fscanf(isodata.fp,"%4d%4d %f %f %f",
			     &grdnx,&grdny,&dx,&dy,&cosxy)) != 5 ) {
	      /* maybe we are at the end of file */
	      while( !((c=fgetc(isodata.fp)) == '\n' || c == EOF) ) {
		;
	      }
	      /* maybe it's end of file */
	      if ( c == EOF ) break;
	      else {
		sprintf(isodata.error, "ERROR: Error reading unit 25");
		return XC_ERROR;
	      }
	    }
	  } 
	else {
	  /*
	    CRYSTAL06 or later
	  */
	  if ( (n = fscanf(isodata.fp,"%5d%5d %f %f %f",
			   &grdnx,&grdny,&dx,&dy,&cosxy)) != 5 ) {
	    /* maybe we are at the end of file */
	    while( !((c=fgetc(isodata.fp)) == '\n' || c == EOF) ) {
	      ;
	    }
	    /* maybe it's end of file */
	    if ( c == EOF ) break;
	    else {
	      sprintf(isodata.error, "ERROR: Error reading unit 25");
	      return XC_ERROR;
	    }
	  }
	}
      }
      c = fgetc(isodata.fp);
      while( !(c== '\n' || c == EOF) ) {
	c = fgetc(isodata.fp);
      }
      /* maybe it's end of file */
      if ( c == EOF ) break;
      /* read second line */
      if ( crystal_version == 95 ) {
	if ( (n = fscanf(isodata.fp,"%4d%4d%f %f %f %f %f",
			 &(grd.nx),&(grd.ny),&x0,&y0,&dx,&dy,&cosxy)) != 7 ) {
	  sprintf(isodata.error, "ERROR: Error reading unit 25");
	  return XC_ERROR;
	}      
      }
      else if ( crystal_version == 98 ) {
	grd.nx = grdnx;
	grd.ny = grdny;
	if ( (n = fscanf(isodata.fp,"%f %f %f %f %f %f %f %f %f",
			 &xa, &ya, &za, &xb, &yb, &zb, &xc, &yc, &zc)) != 9 ) {
	  sprintf(isodata.error, "ERROR: Error reading unit 25");
	  return XC_ERROR;
	}
      } else if ( crystal_version == 3 || crystal_version == 6 || crystal_version == 9 || crystal_version == 14) {
	grd.nx = grdnx;
	grd.ny = grdny;
	if ( (n = fscanf(isodata.fp,"%f %f %f %f %f %f %f %f %f %d %d",
			 &xa, &ya, &za, &xb, &yb, &zb, &xc, &yc, &zc,
			 &natoms, &dim)) != 11 ) {
	  sprintf(isodata.error, "ERROR: Error reading unit 25");
	  return XC_ERROR;
	}
	fprintf(stderr,"MAPN: %f %f %f %f %f %f %f %f %f %d %d\n",
		xa, ya, za, xb, yb, zb, xc, yc, zc,
		natoms, dim);
      }
            
      newgrd.nx = grd.nx;
      newgrd.ny = grd.ny;
      
      /* check if our reading-file is OK */
      if ( !(first_time) && 
	   (oldnx != grd.nx || oldny != grd.ny) ) {
	sprintf(isodata.error, "ERROR: Error reading unit 25");	
	return XC_ERROR;
      } 
      oldnx = grd.nx;
      oldny = grd.ny;
      first_time = 0;
      
      /* read data from isodata.fp and write to binary file */
      for (j=0; j<grd.ny; j++)
	for (k=0; k<grd.nx; k++) {	
	  fscanf(isodata.fp, "%e", &value);
	  fwrite(&value, sizeof(float), 1, isodata.bin_fp);
	  /* fprintf(file,"%e  ", value); */
	}
      
      if ( crystal_version == 3 || crystal_version == 6 || crystal_version == 9 || crystal_version == 14 ) {
	/* read atoms and lattice */
	for (ia=0; ia<natoms; ia++)
	  {
	    if ( fscanf(isodata.fp, "%d %s %f %f %f", &nat, word1, &x_, &y_, &z_) != 5 ) {
	      sprintf(isodata.error, "ERROR: Error reading unit 25");
	      return XC_ERROR;
	    }
	  }
	for (ia=0; ia<9; ia++)
	  fscanf(isodata.fp, "%f", &lat_);
      }

      /* now read till the end of line */
      while( !((c=fgetc(isodata.fp)) == '\n' || c == EOF )) {
	;
      }
      nn++;
    }
  
    if ( nn == isodata.nframe[s3][1] * isodata.nframe[s3][2] ) 
      isodata.nframe[s3][0] = 1;
    else if ( nn == 2 * isodata.nframe[s3][1] * isodata.nframe[s3][2] )
      isodata.nframe[s3][0] = 2;
    else {
      sprintf(isodata.error, "invalid number of blocks in unit 25");
      return TCL_ERROR;
    }
    /* fprintf(file,"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",NULL);
       fclose(file);*/
  }
  else if ( type == RB0_DUMMY || type == RB0_READ ) {    

    /* read data; 
     * don't change the order of loops, because order is important */
    /* CHANGES:: j in k zamenjam v for() */

    for (j=0; j<grd.ny; j++) {
      for (k=0; k<grd.nx; k++) {
	/* fseek(isodata.bin_fp, sizeof(float), SEEK_CUR); */
        fread(&value, sizeof(float), 1, isodata.bin_fp);  
	/* update vertex value */
	if ( type == RB0_READ ) {
	  /* ISO LINE/PLANE */
	  if ( isodata.dim[ISOOBJ_BASE] == 2 ) {
	    plvertex[ISOOBJ_BASE][k][j].val += value * 
	      isodata.framesign[3][s3] * 
	      isodata.framesign[2][s2] *
	      isodata.framesign[1][s1] * 
	      isodata.framesign[0][s0];
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

	    else if ( plvertex[ISOOBJ_BASE][k][j].val < isodata.min )  isodata.min = plvertex[ISOOBJ_BASE][k][j].val;
	    else if ( plvertex[ISOOBJ_BASE][k][j].val > isodata.max )  isodata.max = plvertex[ISOOBJ_BASE][k][j].val;
	  }
	  /* ISO SURFACE */	  
	  else if ( isodata.dim[ISOOBJ_BASE] == 3 ) {	    
	    gridvertex[k][j][s1].val += value * 
	      isodata.framesign[3][s3] * 
	      isodata.framesign[2][s2] *
	      isodata.framesign[1][s1] * 
	      isodata.framesign[0][s0];
	    gridvertex[k][j][s1].p.x =
	      ((float) k  / (grd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].x +  
	      ((float) j  / (grd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].x +  
	      ((float) s1 / (grd.nz - 1)) * isodata.vec[ISOOBJ_BASE][2].x;  
	    gridvertex[k][j][s1].p.y =   
	      ((float) k  / (grd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].y + 
	      ((float) j  / (grd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].y +
	      ((float) s1 / (grd.nz - 1)) * isodata.vec[ISOOBJ_BASE][2].y; 
	    gridvertex[k][j][s1].p.z =   
	      ((float) k  / (grd.nx - 1)) * isodata.vec[ISOOBJ_BASE][0].z +  
	      ((float) j  / (grd.ny - 1)) * isodata.vec[ISOOBJ_BASE][1].z +  
	      ((float) s1 / (grd.nz - 1)) * isodata.vec[ISOOBJ_BASE][2].z;   

	    /*fprintf(fff,"%d %d %d\n%f;       %f %f %f\n",k,j,s1,
		    gridvertex[k][j][s1].val,
		    gridvertex[k][j][s1].p.x,
		    gridvertex[k][j][s1].p.y,
		    gridvertex[k][j][s1].p.z);*/
	    
	    if ( first_t && k == 0 && j == 0 ) {
	      isodata.min = gridvertex[j][k][s1].val;
	      isodata.max = gridvertex[j][k][s1].val;
	    }
	    if ( gridvertex[k][j][s1].val < isodata.min ) 
	      isodata.min = gridvertex[k][j][s1].val;
	    else if ( gridvertex[k][j][s1].val > isodata.max ) 
	      isodata.max = gridvertex[k][j][s1].val;	    	    
	  }
	} /* if */
      } /* for k */
    } /* for j */
  } /* if */

  if ( type == RB0_READ ) first_t = 0;
    
  return XC_OK;  
}
	       
