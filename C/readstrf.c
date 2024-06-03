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
 * Source: $XCRYSDEN_TOPDIR/C/readstrf.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <math.h>
#include <ctype.h>
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tk.h>

#include "struct.h"
#include "3D.h"
#include "bz.h"
#include "vector.h"
#include "memory.h"
#include "xcfunc.h"
#include "getline.h"

#define XSF_OPEN   0
#define XSF_UPDATE 1
#define LINE_FMT   "%[^\n\r]"

StructSize         ss;
extern AtomicLabel *atomLabel, globalAtomLabel;
extern short       *do_not_display_atomlabel;
extern double      rcov[MAXNAT + 1];
extern char        *element[];
extern ForceVector FV;

/* FUNCTION PROTOTYPES */
int ReadStructFile(FILE *fp, const char *file, int format, int mode);
void MakeBonds(void);
void MallocCoor(void);
void FreeCoor(void);
int ReadXSF(FILE *fp);
int ReadVec(FILE *fp, double vec[][4]);
static int ReadCoor( FILE *fp, int natr, int ncell, const int celltype );
static int ReadAtoms( FILE *fp , long int atompos );
static int ReadFrames( FILE *fp , long int framepos );
static int ReadVoronoi(FILE *fp, float poly[][WIGNERSEITZ_MAXVERTEX][3], 
		       float norm[][WIGNERSEITZ_MAXVERTEX][3], int nvert[],
		       float *Max, boolean *xcrl, int *npoly);
void FindMaxRad(void);
static void CorrectVectors( double vec[][4], int dim );
static void CheckAtoms(int *flag1, int *nat1, 
		       double *xat1, double *yat1, double *zat1,
		       double (*fv1)[3]);
static int ReadXYZ(FILE *fp);
static int ReadPDB(FILE *fp);
/*static void StringTrailLeftWhiteSpaces(char *str, size_t size);*/
void ParseAtomName(char *atomName);
static void parsePDBAtomRecord(char *line, int *nat1, 
			       double *x, double *y, double *z);
static void parsePDBHetatmRecord(char *line, int *n, 
				 double *x, double *y, double *z);
static void xcReallocAtomLabels(void);

static void ToEndLine(FILE *fp);

/* extern functions prototypes */
/* --- auxilproc.c ---*/
extern double dist6(double x1, double x2, 
		    double y1, double y2, double z1, double z2);
extern double dist3(double x, double y, double z);

/* --- detnsplit.c --- */
extern int DetNBonds(int *nsplit);

/* --- 3D.c ---*/
extern void MakeCylinderCoor(void);

/* --- datagrid.c --- */
extern void NewGridList( int gridtype, FILE *gridFP );
extern FILE *MaybeOpenDataGridFile(char *mode);
extern int  ReadDataGrid(FILE *fp, FILE *gridFP, int gridtype, char *line);
extern void CloseDataGridFile(void);
extern void SetDataGridCommentLine(char *line);

/* --- xcColorScheme.c --- */
extern void LoadAtmCol(const int flag);

/* -- forces.c -- */
extern void MallocForceVectors(void);

/* -- xcForce.c -- */
extern void BuildForceVectors(ForceVector *fvPtr);

extern void Set_mx_my_mz(double sumx, double sumy, double sumz);

/* xcdebug.c */
extern void xcdebug(const char *text);
extern void xcErrDebug(const char *text);
extern void breakpoint(const char *text);

/* ------------------------------------------------------------------------ */



int 
ReadStructFile(FILE *fp, const char *file, int format, int mode)
{
  register int i, j; 
  double sumx = 0.0, sumy = 0.0, sumz = 0.0, vsize1 = 0.0, vsize2 = 0.0;
  double v1[3], v2[3], v3[3];
  GLboolean latom;

  vec.prim[0][0] = vec.prim[1][1] = vec.prim[2][2] = vec.prim[3][3] =
    vec.conv[0][0] = vec.conv[1][1] = vec.conv[2][2] = vec.conv[3][3] = 1.0;

  vec.prim[0][1] = vec.prim[0][2] = vec.prim[0][3] =
    vec.prim[1][0] = vec.prim[1][2] = vec.prim[1][3] =
    vec.prim[2][0] = vec.prim[2][1] = vec.prim[2][3] =
    vec.prim[3][0] = vec.prim[3][1] = vec.prim[3][2] =
    vec.conv[0][1] = vec.conv[0][2] = vec.conv[0][3] =
    vec.conv[1][0] = vec.conv[1][2] = vec.conv[1][3] =
    vec.conv[2][0] = vec.conv[2][1] = vec.conv[2][3] =
    vec.conv[3][0] = vec.conv[3][1] = vec.conv[3][2] = 0.0;
    
  xcr.celltype = XCR_NOCELL;

  /* reset all xcrys-file related parameters */
  xcr.ldimgroup   = 0;
  xcr.lconvvec    = 0;
  xcr.lprimvec    = 0;
  xcr.lrecconvvec = 0;
  xcr.lrecprimvec = 0;
  xcr.lprimwigner = 0;
  xcr.lconvwigner = 0;
  xcr.lprimbz     = 0;
  xcr.lconvbz     = 0;
  xcr.lconvcoor   = 0;
  xcr.lprimcoor   = 0;
  xcr.lforce      = 0;
  xcr.prim_lforc = 0;
  /*xcr.lHbond      = 0;*/
  if ( mode == XSF_OPEN ) {
    xcr.ldatagrid2D = 0;
    xcr.ldatagrid3D = 0;
    xcr.lbandgrid2D = 0;
    xcr.lbandgrid3D = 0;
  }
  xcr.dim = 0;

  xcr.nunit[0] = 0;
  xcr.nunit[1] = 0;
  xcr.nunit[2] = 0;

  MVf.structsize = 0.0;

  if ( format == FORMAT_XSF )
    if ( !ReadXSF(fp) ) {
      /* just append file name to error massage from ReadXSF(fp) */
      fprintf(stderr,"\"%s\"\n",file);
      return XC_ERROR;
    }
  if ( format == FORMAT_XYZ )
    if ( !ReadXYZ(fp) ) {
      /* just append file name to error massage from ReadXSF(fp) */
      fprintf(stderr,"\"%s\"\n",file);
      return XC_ERROR;
    }
  if ( format == FORMAT_PDB )
    if ( !ReadPDB(fp) ) {
      /* just append file name to error massage from ReadXSF(fp) */
      fprintf(stderr,"\"%s\"\n",file);
      return XC_ERROR;
    }

  /* MOVING TO GEOMETRICAL-MASS-CENTER */
  for(i = 1; i <= natoms; i++)
    {
      sumx += xat[i];
      sumy += yat[i];
      sumz += zat[i];
    }
  
  /*--> THIS ..*/
  Set_mx_my_mz (sumx, sumy, sumz);
  mx = mx_my_mz.mx;
  my = mx_my_mz.my;
  mz = mx_my_mz.mz;
  /* this ... */

  /* initialize max.* */
  max.x = 0.0;
  max.y = 0.0;
  max.z = 0.0;
  max.r = 0.0;

  ss.minX = ss.maxX = mx;
  ss.minY = ss.maxY = my;
  ss.minZ = ss.maxZ = mz;
  atm.natomkind = 0;
  
  for(i = 1; i <= natoms; i++) {

    if ( xat[i] < ss.minX ) ss.minX = xat[i];
    if ( yat[i] < ss.minY ) ss.minY = yat[i];
    if ( zat[i] < ss.minZ ) ss.minZ = zat[i];
    if ( xat[i] > ss.maxX ) ss.maxX = xat[i];
    if ( yat[i] > ss.maxY ) ss.maxY = yat[i];
    if ( zat[i] > ss.maxZ ) ss.maxZ = zat[i];
    
    latom = 0; /* if it's 0 than new element (atm.atomkind[]) */
    xat[i] = xat[i] - mx;
    yat[i] = yat[i] - my;
    zat[i] = zat[i] - mz;  
    
    if (ABS (xat[i]) > max.x) max.x = ABS (xat[i]);
    if (ABS (yat[i]) > max.y) max.y = ABS (yat[i]);
    if (ABS (zat[i]) > max.z) max.z = ABS (zat[i]);

    /* determine which atomic elements are in structure */
    
    for(j=0; j<atm.natomkind; j++) 
      if ( nat[i] == atm.atomkind[j] ) latom = 1;
    if (!latom) {
      atm.atomkind[atm.natomkind] = nat[i];
      atm.natomkind++;
    }
  }
  if ( xcr.shape == XCR_CELL ) {
    if ( xcr.dim == 1 ) {
      max.x   = ss.maxX = 2.0 * mx;
      ss.minX = 0.0;      
    } else if ( xcr.dim == 2 ) {
      max.x   = ss.maxX = 2.0 * mx;
      max.y   = ss.maxY = 2.0 * my;
      ss.minX = ss.minY = 0.0;
    } else if ( xcr.dim == 3 ) {
      max.x   = ss.maxX = 2.0 * mx;
      max.y   = ss.maxY = 2.0 * my;
      max.z   = ss.maxZ = 2.0 * mz;
      ss.minX = ss.minY = ss.minZ = 0.0;      
    }
  }
  FindMaxRad();

  /* --- FRAMES --- FRAMES --- FRAMES --- */
  for(i = 1; i <= nframes; i++)
    {
      *(xframe + i)  -= mx;
      *(yframe + i)  -= my;
      *(zframe + i)  -= mz; 
      *(xframe2 + i) -= mx;
      *(yframe2 + i) -= my;
      *(zframe2 + i) -= mz; 
      if (ABS (*(xframe	 + i)) > max.x) max.x = ABS (*(xframe  + i));
      if (ABS (*(yframe	 + i)) > max.y) max.y = ABS (*(yframe  + i));
      if (ABS (*(zframe	 + i)) > max.z) max.z = ABS (*(zframe  + i));
      if (ABS (*(xframe2 + i)) > max.x) max.x = ABS (*(xframe2 + i));
      if (ABS (*(yframe2 + i)) > max.y) max.y = ABS (*(yframe2 + i));
      if (ABS (*(zframe2 + i)) > max.z) max.z = ABS (*(zframe2 + i));      
    }
  
  /* structsize is longest distance from origin */
  /* structure is sqrt( max.x^2 + max.y^2 + max.z^2 ) 
   */
  MVf.structsize = sqrt( max.x * max.x + max.y * max.y + max.z * max.z );
  
  /* ---- MAKE BONDS ---- */
  MakeBonds();  

  /* ---- MAKE Hydrogen Bonds ---- */
  make_H_Bonds();
  
  /* VERY IMPORTANT !!!!! */
  /* ---- NOW ASSIGN "nobjects" variable ---- */
  nobjects = natoms + nbonds + nframes;
  
  /* malloc "coor" & "coor3D" and stuff that goes nearby */  
  MallocCoor();

  atm.col = (GLfloat (*)[3]) realloc(atm.col, sizeof(float [3]) * (natoms + 1 + N_SELECT_OBJ));
  
  LoadAtmCol(1);
  
  /* --- now calculate coordinates for 3Dbonds & 3Dframes --- */
  MakeCylinderCoor();

  /* --- if we come so far everything is OK --- */  

  /***********************************************************/

  /* correct the structusize also according to "crystal cages"; take
     the largest of PRIMVEC and CONVEC, so that when shifting between
     the two CELL displays the structure will not "size-breath"
   */

  vsize1 = 0.0;
  vsize2 = 0.0;
  
  if ( xcr.lprimvec ) {
    
    CorrectVectors( vec.prim, xcr.dim );
    
    for (i=0; i<3; i++) {
      v1[i] = (xcr.nunit[0] <= 1 ? 1.0 : (double)xcr.nunit[0]) * vec.prim[0][i];
      v2[i] = (xcr.nunit[1] <= 1 ? 1.0 : (double)xcr.nunit[1]) * vec.prim[1][i];
      v3[i] = (xcr.nunit[2] <= 1 ? 1.0 : (double)xcr.nunit[2]) * vec.prim[2][i];
    }
    vsize1 = papipedMaxDiagSize ( v1, v2, v3 );
  }
  if ( xcr.lconvvec ) {
    
    CorrectVectors( vec.conv, xcr.dim );

    for (i=0; i<3; i++) {
      v1[i] = (xcr.nunit[0] <= 1 ? 1.0 : (double)xcr.nunit[0]) * vec.conv[0][i];
      v2[i] = (xcr.nunit[1] <= 1 ? 1.0 : (double)xcr.nunit[1]) * vec.conv[1][i];
      v3[i] = (xcr.nunit[2] <= 1 ? 1.0 : (double)xcr.nunit[2]) * vec.conv[2][i];
    }      
    vsize2 = papipedMaxDiagSize ( v1, v2, v3 );
  }
  
  vsize1 = ( vsize1 > vsize2 ? vsize1 : vsize2 );
  MVf.structsize = ( vsize1 > MVf.structsize ? vsize1 : MVf.structsize );
  
  return XC_OK;
}


void
MakeBonds(void)
{
  int i, j, k, max_ibn, *ibn;
  int numb = 0; /* numb == seq. number of bond; 
		   0 because we increase ++numb, before usage */
  /* double fabs(); */
  double dis, len, fdis;
  int nsplit; /* criteria for additional bond spliting */

  /*C     Find Bonds 
    C     numb.......counter of bonds
    C     nbonds.....number of bonds
    C     dis.....distance between i & j atoms
    C     len.....lenght of sum of covalent radii of i & j atom
    C     fdis....lenght of bond of i atom in i-j bond
    */

  ibn = (int*) malloc (sizeof(int) * (natoms+1));

  /* determine "nbonds" */
  DetNBonds(&nsplit); /* this routine returns also an overestimated "nbonds" */

  max_ibn = nsplit * MAX_BONDS;

  {
    int nb1 = nbonds+1;
    natbond = (int *) calloc(nb1, sizeof(int));
    sqnbond = (int *) calloc(nb1, sizeof(int));
    bondend = (int *) calloc(nb1, sizeof(int));
    
    /* first end of a bond */
    xbond = (double *) calloc(nb1, sizeof(double));
    ybond = (double *) calloc(nb1, sizeof(double));
    zbond = (double *) calloc(nb1, sizeof(double));
    
    /* second end of a bond */
    xbond2 = (double *) calloc(nb1, sizeof(double));
    ybond2 = (double *) calloc(nb1, sizeof(double));
    zbond2 = (double *) calloc(nb1, sizeof(double));
  }

  for(i = 1; i <= natoms; i++) ibn[i]=0;
  
  /* make bonds & also determine if additional spliting is needed */
  for(i = 1; i <= natoms - 1 ; i++) 
    {
      /* 
	 don't make a bond for dummy X-atom 
      */
      if (nat[i] == 0) continue;

      for(j = i + 1; j <= natoms; j++)
	{		
	  /* 
	     don't make a bond for dummy X-atom 
	  */
	  if (nat[j] == 0) continue;
	  
	  /* calculate distance between two atoms */
	  dis = dist6(xat[i], xat[j],
		      yat[i], yat[j],
		      zat[i], zat[j]);
	  /*printf("ReadStructFile> dis = %f !!! \n", dis);
	    fflush(stdout); */   
	  
	  if ( dis < MINDIS ) {
	    if ( dis < MINTOL ) 
	      fprintf(stderr,"WARNING: Atom %d and atom %d overlap !!!\n",i,j);
	    else 
	      fprintf(stderr,
		      "WARNING: Atom %d and atom %d are very close !!!\n",i,j);
	    continue;
	  }
	  len = rcov[nat[i]] + rcov[*(nat + j)];
	  /*printf("ReadStructFile> len = %f !!! \n", len);
	    fflush(stdout);*/
	  if ( dis < len )
	    { 
	      for(k = 0; k < nsplit; k++) {	
		/* only BONDs to ATOM i */

		++ibn[i];

		if ( ibn[i] < max_ibn ) {

		  ++numb;

		  /* FIRST END OF A BOND */
		  *(bondend + numb) = BOND_ATOM_TO_MIDBOND;
		  fdis = rcov[nat[i]] / (len * (double) nsplit);

		  xbond[numb] = xat[i] + (double) k * fdis * (xat[j] - xat[i]);
		  ybond[numb] = yat[i] + (double) k * fdis * (yat[j] - yat[i]);
		  zbond[numb] = zat[i] + (double) k * fdis * (zat[j] - zat[i]);

		  /* this is needed for -O3 compilation option (don't ask me why !!!) */
		  if ( isnan(xbond[numb]) ) {
		    xbond[numb]=xat[i] + (double) k * fdis * (xat[j] - xat[i]);
		  }
		  if ( isnan(ybond[numb]) ) {
		    ybond[numb]=yat[i] + (double) k * fdis * (yat[j] - yat[i]);
		  }
		  if ( isnan(zbond[numb]) ) {
		    zbond[numb]=zat[i] + (double) k * fdis * (zat[j] - zat[i]);
		  }

		  /* SECOND END OF A BOND */
		  xbond2[numb] = xat[i] + (double) (k + 1) * fdis * (xat[j] - xat[i]);
		  ybond2[numb] = yat[i] + (double) (k + 1) * fdis * (yat[j] - yat[i]);
		  zbond2[numb] = zat[i] + (double) (k + 1) * fdis * (zat[j] - zat[i]);

		  /* this is needed for -O3 compilation option (don't ask me why !!!) */
		  if ( isnan(xbond2[numb]) ) {
		    xbond2[numb]=xat[i] + (double) (k+1) * fdis * (xat[j] - xat[i]);
		  }
		  if ( isnan(ybond2[numb]) ) {
		    ybond2[numb]=yat[i] + (double) (k+1) * fdis * (yat[j] - yat[i]);
		  }
		  if ( isnan(zbond2[numb]) ) {
		    zbond2[numb]=zat[i] + (double) (k+1) * fdis * (zat[j] - zat[i]);
		    
		  }

		  /* instead of color a 'bnat' will be assigned !!! */
		  *(natbond + numb) = nat[i];
		  *(sqnbond + numb) = i;

		  /* zbmid == middle Z - comp of bond */
		  /**(zbmid + numb) = 0.5 * (zbond[numb] + zbond2[numb]);*/
		}

		
		/* only BONDs to ATOM j */

		++ibn[j];

		if ( ibn[j] < max_ibn ) {
		  ++numb;

		  /* FIRST END OF A BOND */
		  *(bondend + numb) = BOND_MIDBOND_TO_ATOM;
		  fdis = rcov[*(nat + j)] / (len * nsplit);
		  
		  xbond2[numb] = xat[j] + (double) k * fdis * (xat[i] - xat[j]);
		  ybond2[numb] = yat[j] + (double) k * fdis * (yat[i] - yat[j]);
		  zbond2[numb] = zat[j] + (double) k * fdis * (zat[i] - zat[j]);

		  /* this is needed for -O3 compilation option (don't ask me why !!!) */
		  if ( isnan(xbond2[numb]) ) {
		    xbond2[numb]=xat[j] + (double) k * fdis * (xat[i] - xat[j]);
		  }
		  if ( isnan(ybond2[numb]) ) {
		    ybond2[numb]=yat[j] + (double) k * fdis * (yat[i] - yat[j]);
		  }
		  if ( isnan(zbond2[numb]) ) {
		    zbond2[numb]=zat[j] + (double) k * fdis * (zat[i] - zat[j]);
		  }
		  
		  /* SECOND END OF A BOND */
		  xbond[numb] = xat[j] + (double) (k + 1) * fdis * (xat[i] - xat[j]);
		  ybond[numb] = yat[j] + (double) (k + 1) * fdis * (yat[i] - yat[j]);
		  zbond[numb] = zat[j] + (double) (k + 1) * fdis * (zat[i] - zat[j]);

		  /* this is needed for -O3 compilation option (don't ask me why !!!) */
		  if ( isnan(xbond[numb]) ) {
		    xbond[numb]=xat[j] + (double) (k + 1) * fdis * (xat[i] - xat[j]);
		  }
		  if ( isnan(ybond[numb]) ) {
		    ybond[numb]=yat[j] + (double) (k + 1) * fdis * (yat[i] - yat[j]);
		  }
		  if ( isnan(zbond[numb]) ) {
		    zbond[numb]=zat[j] + (double) (k + 1) * fdis * (zat[i] - zat[j]);
		  }
		  
		  /* instead of color a 'bnat' will be assigned !!!!!!!!!*/
		  *(natbond + numb) = *(nat + j);
		  *(sqnbond + numb) = j;
		  
		  /* zbmid == middle Z - comp of bond */
		  /**(zbmid + numb) = 0.5 * (zbond[numb] + zbond2[numb]);*/
		}
	      } /* for(k=...) */
	    } /* if (dis < len) */
	} /* for(j=...) */
    }
  
  /* ---- ASSIGN NBONDS ---- */
  nbonds = numb;

  free(ibn);
}

/* ---- malloc coor & things that go nearby ---- */
void
MallocCoor(void)
{
  int maxobj;

  maxobj = nobjects + N_SELECT_OBJ + 1;

  /* allocate memory for *coor */
  coor = (AtomBond *) malloc( sizeof(AtomBond) * maxobj);  
  zorient = (double *) malloc( sizeof(double) * maxobj);
  iwksp = (int *) malloc( sizeof(int) * maxobj);

  /* allocate also memory for BondFrame3D coor3D */
  coor3D = (BondFrame3D *) malloc( sizeof(BondFrame3D) * 
				   (nbonds + nframes + 1) );
}


/* ========================================================================= */
/*                          READ XCRYSDEN STRUCTURE FILE                     */
/* ========================================================================= */
int
ReadXSF(FILE *fp)
{
  FILE *gridFP;
  int c;
  int latom = 0, lframe = 0, lend = 0;
  int n, idummy;
  long int atompos = 0;  /* position where atom coor starts in the file */
  long int framepos = 0; /* position where frame coor starts in the file */
  size_t linesize=256;
  static char *line = 0;

  if (line==0) {
    line = malloc( sizeof(char) * linesize );
  }

  while((c = fscanf(fp,"%s\n",line)) != EOF) {
    lend = 0;

    /* if # (comment character -- ignore line) */
    if ( strncmp(line,"#",1) == 0) {
      /* read the rest of the line */
      fscanf(fp,"%[^\n]",line);
      continue;
    }

    /* if INFO */
    if ( strncmp(line,"INFO",4) == 0 || 
	 strncmp(line,"BEGIN_INFO", 8) == 0) {
      char word[80];
      latom  = 0; 
      lframe = 0;
      lend   = 1;
      do {
	fscanf(fp,"%s\n",line);
	
	if ( strncmp(line,"title",5) == 0 ) {
	  /* this is dummy, but allowed in order to allow user specified
	     some title/comment for him-self */
	  ;
	}
	else if ( strncmp(line,"nunit",5) == 0 ) {
	  n = fscanf(fp,"%d %d %d", 
		     &xcr.nunit[0], &xcr.nunit[1], &xcr.nunit[2]);
	  if ( n != 3 ) {
	    fprintf(stderr, "error parsing INFO.nunit ");
	    return XC_ERROR;
	  }
	}
	else if ( strncmp(line,"shape",5) == 0 ) {
	  if ( (n = fscanf(fp,"%s", word)) != 1 ) {
	    fprintf(stderr, "error parsing INFO.shape ");	 
	    return XC_ERROR;
	  }
	  if ( strncmp(word,"para",4) == 0 ) 
	    xcr.shape = XCR_PARAPIPEDAL;	  
	  else if ( strncmp(word,"hexa",4) == 0 ) 
	    xcr.shape = XCR_HEXAGONAL;
	  else {
	    fprintf(stderr, "error parsing INFO.shape value ");	  
	    return XC_ERROR;
	  }
	}
	else if ( strncmp(line,"unit",4) == 0 ) {
	  if ( (n = fscanf(fp,"%s", word)) != 1 ) {
	    fprintf(stderr, "error parsing INFO.unit ");	  
	    return XC_ERROR;
	  }
	  if ( strncmp(word,"cell",4) == 0 ) 
	    xcr.shape = XCR_CELL;
	  else if ( strncmp(word,"tr_a",4) == 0 ) 
	    xcr.shape = XCR_TR_ASYM;
	  else {
	    fprintf(stderr, "error parsing INFO.unit section ");	  
	    return XC_ERROR;
	  }
	}
	else if ( strncmp(line,"celltype",8) == 0 ) {
	  if ( (n = fscanf(fp,"%s", word)) != 1 ) {
	    fprintf(stderr, "error parsing INFO.celltype ");	  
	    return XC_ERROR;
	  }
	  if ( strncmp(word,"prim",4) == 0 ) /* primitiv cell; "primcell" */
	    xcr.celltype = XCR_PRIMCELL;
	  else if ( strncmp(word,"conv",4) == 0 ) /* "convcell" */
	    xcr.celltype = XCR_CONVCELL;
	  else {
	    fprintf(stderr, "error parsing INFO.celltype section ");	  
	    return XC_ERROR;
	  }
	}	
	else if ( strncmp(line,"#",1) == 0 ) {
	  /* comment inside INFO record; read the rest of the line */
	  fscanf(fp,"%[^\n]",line);
	}
	else if ( strncmp(line, "Fermi", 5) == 0 ) {
	  /* this is used in BXSF for reading the fermi energy; read
	     the rest of the line */
	  fscanf(fp,"%[^\n]",line);
	}
	else if ( strncmp(line,"END_INFO",8) == 0 ||
		  strncmp(line,"INFO_END",8) == 0) {
	  ;
	}
	else {
	  fprintf(stderr, "error parsing INFO section ");	  
	  return XC_ERROR;
	}
      } while( strncmp(line,"END_INFO",8) != 0 &&
	       strncmp(line,"INFO_END",8) != 0 );
    }
    
    /* if DIM-GROUP */
    else if ( strncmp(line,"DIM",3) == 0 ) {
      latom = 0; 
      lframe = 0;
      /* read dimensionality and group of structure */
      n = fscanf(fp,"%d %d", &xcr.dim, &xcr.groupn);
      if ( n != 2 ) {
	fprintf(stderr,
		"ERROR: Error reading DIM-GROUP section, while reading ");
	/*xcFree((void*)line);*/
	return XC_ERROR;
      }
      if ( xcr.dim < 0 || xcr.dim > 3 ) {
	fprintf(stderr,	"ERROR: DIM is out of range, while reading ");
	/*xcFree((void*)line);*/
	return XC_ERROR;
      }
      xcr.ldimgroup = 1;
    }    	
    
    /* MOLECULE */
    else if ( strncmp(line,"MOLECULE",8) == 0 ) {
      /* do nothing */
      latom = 0; 
      lframe = 0;
    }
	
    /* POLYMER */
    else if ( strncmp(line,"POLYMER",7) == 0 ) {
      latom = 0; 
      lframe = 0;
      xcr.dim    = 1;
      xcr.groupn = 1;
      xcr.ldimgroup = 1;
    }

    /* SLAB */
    else if ( strncmp(line,"SLAB",4) == 0 ) {
      latom = 0; 
      lframe = 0;
      xcr.dim    = 2;
      xcr.groupn = 1;
      xcr.ldimgroup = 1;
    }
    
    /* CRYSTAL */
    else if ( strncmp(line,"CRYSTAL",7) == 0 ) {
      latom = 0; 
      lframe = 0;
      xcr.dim    = 3;
      xcr.groupn = 1;
      xcr.ldimgroup = 1;
    }

    /* if PRIMVEC */
    else if ( strncmp(line,"PRIMV",5) == 0 ) {
      latom = 0; 
      lframe = 0;
      if ( !ReadVec( fp, vec.prim ) ) {
	fprintf(stderr, "ERROR: Error reading PRIMVEC section, while reading ");
	/*xcFree((void*)line);*/
	return XC_ERROR;
      }
      /* we have read the primitive vectors, so xcr.lprimvec must be set to 1
       */
      xcr.lprimvec = 1;
    }

    /* if CONVVEC */
    else if ( strncmp(line,"CONVV",5) == 0 ) {
      latom = 0; 
      lframe = 0;
      if ( !ReadVec( fp, vec.conv ) ) {	
	fprintf(stderr, "ERROR: Error reading CONVVEC section, while reading ");
	/*xcFree((void*)line);*/
	return XC_ERROR;
      }
      xcr.lconvvec = 1;
    }

    /* if RECIP-PRIMVEC */
    else if ( strncmp(line,"RECIP-PRIMV",11) == 0 ) {
      latom = 0; 
      lframe = 0;
      if ( !ReadVec( fp, vec.recprim ) ) {
	fprintf(stderr, "ERROR: Error reading RECIP-PRIMVEC section, while reading ");
	/*xcFree((void*)line);*/
	return XC_ERROR;
      }
      /* we have read the primitive vectors, so xcr.lprimvec must be set to 1
       */
      xcr.lrecprimvec = 1;
    }

    /* if RECIP-CONVVEC */
    else if ( strncmp(line,"RECIP-CONVV",11) == 0 ) {
      latom = 0; 
      lframe = 0;
      if ( !ReadVec( fp, vec.recconv ) ) {	
	fprintf(stderr,
		"ERROR: Error reading RECIP-CONVVEC section, while reading ");
	/*xcFree((void*)line);*/
	return XC_ERROR;
      }
      xcr.lrecconvvec = 1;
    }

    /* if WIGNER-SEITZ-PRIMCELL */
    else if ( strncmp(line,"WIGNER-SEITZ-PRIMCELL",21) == 0 ) {
      latom = 0; 
      lframe = 0;
      if ( ReadVoronoi(fp, wsp.poly, wsp.norm, wsp.nvert, 
		       &(wsp.max), &(xcr.lprimwigner), &wsp.npoly) == XC_ERROR ) {
	fprintf(stderr, "error parsing WIGNER-SEITZ-PRIMCELL ");	  
	return XC_ERROR;
      }
    }

    /* if WIGNER-SEITZ-CONVCELL */
    else if ( strncmp(line,"WIGNER-SEITZ-CONVCELL",21) == 0 ) {
      latom = 0; 
      lframe = 0;
      if ( ReadVoronoi(fp, wsc.poly, wsc.norm, wsc.nvert, 
		       &(wsc.max), &(xcr.lconvwigner), &wsc.npoly) == XC_ERROR ) {
	fprintf(stderr, "error parsing WIGNER-SEITZ-CONVCELL ");
	return XC_ERROR;
      }
    }

    /* if BRILLOUIN-ZONE-PRIMCELL */
    else if ( strncmp(line,"BRILLOUIN-ZONE-PRIMCELL",23) == 0 ) {
      latom = 0; 
      lframe = 0;
      if ( ReadVoronoi(fp, bz[BZ_PRIMCELL].poly, bz[BZ_PRIMCELL].norm, 
		       bz[BZ_PRIMCELL].nvert, &(bz[BZ_PRIMCELL].max), 
		       &(xcr.lprimbz), &(bz[BZ_PRIMCELL].npoly)) == XC_ERROR ) {
	fprintf(stderr, "error parsing BRILLOUIN-ZONE-PRIMCELL ");	  
	return XC_ERROR;
      }
    }

    /* if BRILLOUIN-ZONE-CONVCELL */
    else if ( strncmp(line,"BRILLOUIN-ZONE-CONVCELL",23) == 0 ) {
      latom = 0; 
      lframe = 0;
      if ( ReadVoronoi(fp, bz[BZ_CONVCELL].poly, bz[BZ_CONVCELL].norm, 
		       bz[BZ_CONVCELL].nvert, &(bz[BZ_CONVCELL].max), 
		       &(xcr.lconvbz), &(bz[BZ_CONVCELL].npoly)) == XC_ERROR ) {
	fprintf(stderr, "error parsing BRILLOUIN-ZONE-CONVCELL ");	  
	return XC_ERROR;
      }
    }

    /* if PRIMCOORD */
    else if ( strncmp(line,"PRIMC",5) == 0 ) {
      latom = 0; 
      lframe = 0;
      n = fscanf(fp,"%d %d", &xcr.natr, &idummy);
      if ( n != 2 ) {
	fprintf(stderr,
		"ERROR: Error reading PRIMCOORD section, while reading ");
	/*xcFree((void*)line);*/
	return XC_ERROR;
      }
      ToEndLine(fp);      
      fprintf(stderr,"reading PRIMCOORD\n");
      if ( ReadCoor( fp, xcr.natr, 1, XCR_PRIMCELL ) == XC_ERROR ) {
	/*xcFree((void*)line);*/
	fprintf(stderr, "error parsing PRIMCOORD ");
	return XC_ERROR;
      }
      lend = 1;
      xcr.lprimcoor = 1;
    }

    /* if CONVCOORD */
    else if ( strncmp(line,"CONVC",5) == 0 ) {
      latom = 0; 
      lframe = 0;
      n = fscanf(fp,"%d %d", &xcr.natr, &xcr.ncell);
      if ( n != 2 ) {
	fprintf(stderr,
		"ERROR: Error reading CONVCOORD section, while reading ");
	/*xcFree((void*)line);*/
	return XC_ERROR;
      }
      ToEndLine(fp);      
      fprintf(stderr,"reading CONVCOORD\n");
      if ( ReadCoor( fp, xcr.natr, xcr.ncell, XCR_CONVCELL ) == XC_ERROR ) {
	/*xcFree((void*)line);*/
	fprintf(stderr, "error parsing CONVCOORD ");
	return XC_ERROR;
      }
      lend = 1;
      xcr.lconvcoor = 1;
    }    

    /* if ATOM */
    else if ( strncmp(line,"ATOM",4) == 0 ) {
      /* remember position */
      atompos = ftell(fp);
      latom = 1;
      lframe = 0;
    }

    /* if FRAME */
    else if ( strncmp(line,"FRAM",4) == 0 ) {
      /* remember position */
      framepos = ftell(fp);
      latom = 0;
      lframe = 1;
    }

    /* VECTORS3D */
/*     else if ( strncmp(line,"VECTORS3D", 14) == 0 ) { */
/*       latom  = 0; */
/*       lframe = 0; */
/*       if ( (n = fscanf(fp, "%d", &nvectors)) != 1 ) { */
/* 	fprintf(stderr, */
/* 		"ERROR: Error reading VECTORS3D section, while reading "); */
/* 	xcFree((void*)line);       */
/* 	return XC_ERROR; */
/*       } */
/*       ReadVectors3D (fp, nvectors); */
/*     } */

    /* read 2D DATAGRID */
    else if ( strncmp(line, "BEGIN_BLOCK_DATAGRID2D", 22) == 0 ||
	      strncmp(line, "BEGIN_BLOCK_DATAGRID_2D", 23) == 0 ) {
      latom = 0;  /* prevent counting of atoms and frames */
      lframe = 0;
      lend = 1;
      /* read comment line */
      if ( (c = fscanf(fp, LINE_FMT, line)) == EOF) {
	fprintf(stderr, "error parsing BLOCK_DATAGRID2D section (#1) ");	  
	return XC_ERROR;
      }
      SetDataGridCommentLine( line );
      gridFP = MaybeOpenDataGridFile("w");
      NewGridList( DATAGRID_2D, gridFP );
      /* read DATAGRID_2D_XXX line */
      if ( (c = fscanf(fp,"%s\n",line)) == EOF) {
	fprintf(stderr, "error parsing BLOCK_DATAGRID_2D section (#2) ");	  
	return XC_ERROR;
      }
      while (strncmp(line, "END_BLOCK_DATAGRID2D", 20) != 0 &&
	     strncmp(line, "END_BLOCK_DATAGRID_2D", 21) != 0 ) {
	if ( strncmp(line, "DATAGRID_2D_", 12) == 0 ||
	     strncmp(line, "DATAGRID2D_", 11) == 0 ||
	     strncmp(line, "BEGIN_DATAGRID_2D", 17) == 0 ) {      
	  if ( ReadDataGrid( fp, gridFP, DATAGRID_2D, line) == XC_ERROR ) {
	    fprintf(stderr, "ERROR: Error reading DATAGRID_2D_ section, while reading ");
	    /*xcFree((void*)line);*/
	    return XC_ERROR;
	  }
	  /* read END_DATAGRID_2D line */
	  breakpoint("breakpoint");
	  if ( (c = fscanf(fp,"%s\n",line)) == EOF) {
	    fprintf(stderr, "error parsing BLOCK_DATAGRID_2D section (#3) ");	  
	    return XC_ERROR;
	  }
	  if ( strncmp(line, "END_DATAGRID_2D", 15) != 0 &&
	       strncmp(line, "END_DATAGRID2D", 14) != 0 ) {
	    fprintf(stderr, "expecting END_DATAGRID_2D keyword, but it wasn't found ");
	    return XC_ERROR;
	  }
	} else {
	  fprintf(stderr, "error parsing BLOCK_DATAGRID_2D section (#4) ");
	  return XC_ERROR;
	}
	/* read next line */ 
	if ( (c = fscanf(fp,"%s\n",line)) == EOF) {
	  fprintf(stderr, "premature end of file in BLOCK_DATAGRID_2D section ");
	  return XC_ERROR;
	}
      }
      xcr.ldatagrid2D = 1;
    }

    /* read 3D DATAGRID */
    else if ( strncmp(line, "BEGIN_BLOCK_DATAGRID3D", 22) == 0 ||
	      strncmp(line, "BEGIN_BLOCK_DATAGRID_3D", 23) == 0) {
      latom = 0;  /* prevent counting of atoms and frames */
      lframe = 0;
      lend = 1;
      /* read comment line */
      if ( (c = fscanf(fp, LINE_FMT, line)) == EOF) {
	fprintf(stderr, "premature end of file in BLOCK_DATAGRID_3D section ");
	return XC_ERROR;
      }
      SetDataGridCommentLine( line );
      gridFP = MaybeOpenDataGridFile("w");
      NewGridList( DATAGRID_3D, gridFP );
      /* read DATAGRID_3D_XXX line */
      if ( (c = fscanf(fp,"%s\n",line)) == EOF) {
	fprintf(stderr, "premature end of file in BLOCK_DATAGRID_3D section ");
	return XC_ERROR;
      }
      while (strncmp(line, "END_BLOCK_DATAGRID3D", 20) != 0 &&
	     strncmp(line, "END_BLOCK_DATAGRID_3D", 21) != 0) {
	if ( strncmp(line, "DATAGRID_3D_", 12) == 0 ||
	     strncmp(line, "DATAGRID3D_", 11) == 0 ||
	     strncmp(line, "BEGIN_DATAGRID_3D", 17) == 0 ) {      
	  if ( ReadDataGrid( fp, gridFP, DATAGRID_3D, line) == XC_ERROR ) {
	    fprintf(stderr, "ERROR: Error reading DATAGRID_3D_ section, while reading ");
	    /*xcFree((void*)line);*/
	    return XC_ERROR;
	  }
	  /* read END_DATAGRID_3D line */
	  breakpoint("breakpoint");
	  if ( (c = fscanf(fp,"%s\n",line)) == EOF) {
	    fprintf(stderr, "premature end of file in BLOCK_DATAGRID_3D section ");	    
	    return XC_ERROR;
	  }
	  if ( strncmp(line, "END_DATAGRID_3D", 15) != 0 &&
	       strncmp(line, "END_DATAGRID3D", 14) != 0 ) {
	    fprintf(stderr, "expecting END_DATAGRID_3D keyword, but it wasn't found ");
	    return XC_ERROR;
	  }
	} else {
	  fprintf(stderr, "error parsing BLOCK_DATAGRID_3D section ");
	  return XC_ERROR;
	}
	/* read next line */ 
	if ( (c = fscanf(fp,"%s\n",line)) == EOF) {
	  fprintf(stderr, "premature end of file in BLOCK_DATAGRID_3D section ");
	  return XC_ERROR;
	}
      }
      xcr.ldatagrid3D = 1;
    }

    /*************************************************************************/
    /* read 3D BANDGRID                                                      */
    /*************************************************************************/
    else if ( strncmp(line, "BEGIN_BLOCK_BANDGRID3D", 22) == 0 ||
	      strncmp(line, "BEGIN_BLOCK_BANDGRID_3D", 22) == 0) {
      latom  = 0;  /* prevent counting of atoms and frames */
      lframe = 0;
      lend   = 1;
      /* read comment line */
      if ( (c = fscanf(fp,"%s\n",line)) == EOF) {
	fprintf(stderr, "premature end of file in BLOCK_BANDGRID_3D");
	return XC_ERROR;
      }
      SetDataGridCommentLine( line );      
      gridFP = MaybeOpenDataGridFile("w");
      NewGridList( DATAGRID_3D, gridFP );      
      /* read DATAGRID_3D_XXX line */
      if ( (c = fscanf(fp,"%s\n",line)) == EOF) {
	fprintf(stderr, "premature end of file in BLOCK_BANDGRID_3D");
	return XC_ERROR;
      }
      while (strncmp(line, "END_BLOCK_BANDGRID3D", 15) != 0 &&
	     strncmp(line, "END_BLOCK_BANDGRID_3D", 16) != 0) {
	if ( strncmp(line, "BANDGRID_3D", 11) == 0 
	     || strncmp(line, "BEGIN_BANDGRID_3D", 17) == 0 
	     || strncmp(line, "BANDGRID3D", 10) == 0 ) {      
	  if ( ReadBandGrid( fp, gridFP, DATAGRID_3D, line) 
	       == XC_ERROR ) {
	    fprintf(stderr, "ERROR: Error reading BANDGRID_3D_ section, while reading ");
	    /*xcFree((void*)line);*/
	    return XC_ERROR;
	  }
	  /* read END_DATAGRID_3D line */
	  /*breakpoint("breakpoint");*/
	  if ( (c = fscanf(fp,"%s\n",line)) == EOF) {
	    fprintf(stderr, "premature end of file in BLOCK_BANDGRID_3D");
	    return XC_ERROR;
	  }
	  if (strncmp(line, "END_BANDGRID_3D", 15) != 0 &&
	      strncmp(line, "END_BANDGRID3D", 14) != 0) {
	    fprintf(stderr, "expecting END_BANDGRID_3D keyword, but it wasn't found ");
	    return XC_ERROR;
	  }
	} else {
	  fprintf(stderr, "error parsing END_BANDGRID_3D section ");
	  return XC_ERROR;
	}
	/* read next line */ 
	if ( (c = fscanf(fp,"%s\n",line)) == EOF) {
	  fprintf(stderr, "premature end of file in BLOCK_BANDGRID_3D");
	  return XC_ERROR;
	}
      }
      xcr.lbandgrid3D = 1;
    }
    /*************************************************************************/

    if ( latom ) ++natoms;
    if ( lframe ) ++nframes;

    if ( !lend )
      while((c = getc(fp)) != '\n' && c != EOF)
	;
  }
  /* if datagrid file was opended, close it */
  CloseDataGridFile();

  /* generate prim_fcoor and conv_fcoor (i.e. fractional coordinates) */
  if ( xcr.lrecprimvec && xcr.lprimcoor ) {
    double a, b, c;
    register int i;

    xcr.prim_fcoor = (double (*)[3]) 
      realloc(xcr.prim_fcoor,(size_t)xcr.natr * sizeof(double [3]));

    for (i=0; i<xcr.natr; i++) {
      /* use transposed reciprocal vectors */
      a = 
	xcr.prim_coor[i][0]*vec.recprim[0][0] + 
	xcr.prim_coor[i][1]*vec.recprim[0][1] + 
	xcr.prim_coor[i][2]*vec.recprim[0][2];
      b = 
	xcr.prim_coor[i][0]*vec.recprim[1][0] + 
	xcr.prim_coor[i][1]*vec.recprim[1][1] + 
	xcr.prim_coor[i][2]*vec.recprim[1][2];
      c = 
	xcr.prim_coor[i][0]*vec.recprim[2][0] + 
	xcr.prim_coor[i][1]*vec.recprim[2][1] + 
	xcr.prim_coor[i][2]*vec.recprim[2][2];
      /* fractional coordinates in range [0,1] */
      if (a < 0.0) a+=1.0;
      if (b < 0.0) b+=1.0;
      if (c < 0.0) c+=1.0;
      xcr.prim_fcoor[i][0] = a;
      xcr.prim_fcoor[i][1] = b;
      xcr.prim_fcoor[i][2] = c;
    }
  }
  if ( xcr.lrecconvvec && xcr.lconvcoor ) {
    double a, b, c;
    register int i, size;

    size = xcr.natr * xcr.ncell;
    /* xcr.conv_fcoor = (double (*)[3]) 
       realloc(xcr.conv_fcoor,(size_t)size * sizeof(double [3])); */
    
    for (i=0; i<size; i++) {
      /* use transposed reciprocal vectors */
      a = 
	xcr.conv_fcoor[i][0]*vec.recconv[0][0] + 
	xcr.conv_fcoor[i][1]*vec.recconv[0][1] + 
	xcr.conv_fcoor[i][2]*vec.recconv[0][2];
      b = 				      
	xcr.conv_fcoor[i][0]*vec.recconv[1][0] + 
	xcr.conv_fcoor[i][1]*vec.recconv[1][1] + 
	xcr.conv_fcoor[i][2]*vec.recconv[1][2];
      c = 				      
	xcr.conv_fcoor[i][0]*vec.recconv[2][0] + 
	xcr.conv_fcoor[i][1]*vec.recconv[2][1] + 
	xcr.conv_fcoor[i][2]*vec.recconv[2][2];
      /* fractional coordinates in range [0.1] */
      if (a < 0.0) a+=1.0;
      if (b < 0.0) b+=1.0;
      if (c < 0.0) c+=1.0;
      xcr.conv_fcoor[i][0] = a;
      xcr.conv_fcoor[i][1] = b;
      xcr.conv_fcoor[i][2] = c;
    }
  }
  /* put somewhere the numbers of atoms */
  fprintf(stderr,"Number of Atoms:  %d\n", natoms);
  fprintf(stderr,"Number of Frames: %d\n", nframes);
  
  /* now read atoms */
  if ( natoms > 0 )
    if (!ReadAtoms( fp, atompos )) {
      fprintf(stderr, "error reading Atoms ");	  
      return XC_ERROR;
    }

  /* if nframes > 0 also read frames */
  if ( nframes > 0 ) 
    if (!ReadFrames( fp, framepos )) {
      fprintf(stderr, "error reading Frames ");	  
      return XC_ERROR;
    }
        
  /*xcFree((void*)line);*/

  current_file_format = FORMAT_XSF;

  return XC_OK;
}


int
ReadVec(FILE *fp, double vec[][4])
{
  int n0, n1, n2;

  /* LoadIdentity( vec ); */
  /* lets read vectors in ROW-MAJOR mode */
  n0 = fscanf(fp,"%lf %lf %lf",&vec[0][0], &vec[0][1], &vec[0][2]);
  n1 = fscanf(fp,"%lf %lf %lf",&vec[1][0], &vec[1][1], &vec[1][2]);
  n2 = fscanf(fp,"%lf %lf %lf",&vec[2][0], &vec[2][1], &vec[2][2]);
    
  if ( n0 == 3 && n1 == 3 && n2 == 3 ) return XC_OK;
  
  else {
    fprintf(stderr, "error reading lattice vectors ");
    return XC_ERROR;
  }
}


/* read PRIMCOORD/CONVCOORD */
static int
ReadCoor( FILE *fp, int natr, int ncell, const int celltype )
{
  int i, n, size;
  size_t linesize=256;
  double (*coor)[3], (*forc)[3];
  int  *atn;
  static char *line = 0;

  if (line==0) {
    line = malloc(sizeof(char)*linesize);
  }
  
  /* natr  -- number of atoms in prim. cell */
  /* ncell -- number of particle per cell (FC=4, RC=3) */

  size = natr * ncell;
  if ( celltype == XCR_PRIMCELL ) {
    coor = (double (*)[3]) 
      realloc(xcr.prim_coor,(size_t)size * sizeof(double [3]));
    forc = (double (*)[3]) 
      realloc(xcr.prim_forc,(size_t)size * sizeof(double [3]));
    atn  = (int *)
      realloc(xcr.prim_nat,(size_t)size * sizeof(int));
  } else if ( celltype == XCR_CONVCELL ) {
    /* 
       BEWARE: 
       
       for the time being we never use the conv_coor, hence we will
       read conv_coor into conv_fcoor and transform them later to
       fractional coordinates !!!
    */
    coor = (double (*)[3]) 
      realloc(xcr.conv_fcoor,(size_t)size * sizeof(double [3]));
    atn  = (int *)
      realloc(xcr.conv_nat,(size_t)size * sizeof(int));
  }    

  for(i=0; i<size; i++) 
    {
      getline(&line, &linesize, fp);
      n = sscanf(line,"%d %lf %lf %lf %lf %lf %lf", 
		 &atn[i], &coor[i][0], &coor[i][1], &coor[i][2],
		 &forc[i][0], &forc[i][1], &forc[i][2]);
            
      if ( n < 4 ) {
	fprintf(stderr,"ERROR: Error reading XSF file (n=%d)",n);
	return XC_ERROR;
      }      
      if (n==7 && celltype == XCR_PRIMCELL) xcr.prim_lforc=1;
    }

  if ( celltype == XCR_PRIMCELL ) {
    xcr.prim_coor = coor;
    xcr.prim_forc = forc;
    xcr.prim_nat  = atn;
  } else if ( celltype == XCR_CONVCELL ) {
    xcr.conv_fcoor = coor;
    xcr.conv_nat  = atn;
  }    

  return XC_OK;
}


static int
ReadAtoms( FILE *fp , long int atompos )
{
  int sizeIA, sizeDA; /* for size of malloc */
  int i, n, j, k;
  int ndel;  
  double dis;
  double *xat1, *yat1, *zat1, (*fv1)[3];
  int *nat1, *flag1;
  size_t linesize = 256;
  static char *line = 0;

  if (line == 0) {
    line = (char *) malloc( sizeof(char) * linesize );
  }

  /* go to atoms position !!!!! */
  fseek(fp, atompos, 0);

  if ( natoms == 0 )
    {
      fprintf(stderr,"WARNING: No atoms find !!!   File: ");
      return XC_ERROR;
    }
  /* read NAT,X,Y,Z from file to *xat1, *yat1, *zat1, *nat1*/
  /* allocate memory for NAT,X,Y,Z */ 
  sizeIA = sizeof(int) * (natoms + 1);
  sizeDA = sizeof(double) * (natoms + 1);
    
  /*sqn = (int *) malloc(sizeIA);*/
  nat1  = (int *)    malloc(sizeIA);
  flag1 = (int *)    malloc(sizeIA);
  xat1  = (double *) malloc(sizeDA);
  yat1  = (double *) malloc(sizeDA);
  zat1  = (double *) malloc(sizeDA);
  fv1   = (double (*)[3]) xcCalloc( (size_t) natoms + 1, sizeof(double [3]) );

  /* first read coordinates and nat; use *xat1, ..., *nat1 variables,
   * than determine if some of atoms overlap, and delete overlaping atoms,
   * and than write to *xat, ..., *nat
   */
  for(i = 1; i <= natoms; i++) 
    {      
      *(flag1 + i) = 1;
      getline(&line,&linesize,fp);
      n = sscanf(line,"%d %lf %lf %lf %lf %lf %lf",
		 nat1 + i, xat1 + i, yat1 + i, zat1 + i,
		 &fv1[i][0], &fv1[i][1], &fv1[i][2]);
      if ( n < 4 ) {
	fprintf(stderr,"ERROR: Error reading file ");
	return XC_ERROR;
      }      

      if (n==7) xcr.lforce=1;

      /* In crystal95 atomic number can be greater than 100, for example:
       * 208 means oxygen with pseudo-potential
       * Thatwhy:
       */
      nat1[i] =  nat1[i] % 100;
 
      /*if ( nat1[i] > MAXNAT ) {
       *fprintf(stderr,"WARNING: Atomic number greater then %d !!!   File: ", 
       *MAXNAT);
       *return XC_ERROR;
       *}
       */
    }

  /* if two atoms overlap, than we will have troubles when making bonds;
   * if two atoms are to close just disregard one of them
   */
  ndel=0;
  for(i = 1; i<= natoms - 1; i++)
    {
      if ( nat1[i] == 0 ) continue; /* dummy atom */
      if (*(flag1 + i)) {
	for(j = i + 1; j <= natoms; j++) {
	  if ( flag1[j] == 0 ) continue; /* atom already deleted */
	  if ( nat1[j]  == 0 ) continue; /* dummy atom */
	  dis = dist6(xat1[i], *(xat1 + j),
		      yat1[i], *(yat1 + j),
		      zat1[i], *(zat1 + j));
	  if ( dis < MINDIS ) {
	    if ( dis < MINTOL ) 
	      fprintf(stderr,"WARNING: Atom %d and atom %d overlap !!!",i,j);
	    else 
	      fprintf(stderr,
		      "WARNING: Atom %d and atom %d are very close !!!",i,j);
	    fprintf(stderr,"   Atom %d deleted !!!\n",j);
	    ++ndel;
	    *(flag1 + j) = 0;
	  }
	}
      }
    }
  
  fprintf(stderr,"natoms: %d , ndel: %d\n", natoms, ndel);
  fprintf(stderr,"Filtered number of atoms: %d\n", natoms-ndel);

  /* allocate memory for NAT,X,Y,Z */ 
  sizeIA = sizeof(int) * (natoms - ndel + 1);
  sizeDA = sizeof(double) * (natoms - ndel + 1);
    
  nat = (int *)    malloc(sizeIA);
  xat = (double *) malloc(sizeDA);
  yat = (double *) malloc(sizeDA);
  zat = (double *) malloc(sizeDA);
  /*fv  = (double (*)[3]) malloc((size_t) (sizeDA * 3));*/
  fv  = (double (*)[3]) xcCalloc((size_t)  (natoms - ndel + 1), sizeof(double [3]) );

  /*
    TODO: this should be reallocated, but then we need to feed in the zeroes,
          make a routine that will do that
  */
  xcReallocAtomLabels();  

  /*  
      for(i=0; i<N_LIST_PRIMITIVES; i++)
      primitive_Lists[i] = (GLuint *) 
      malloc((size_t) (sizeof(GLuint) * (natoms - ndel + 1)) );
  */

  /* write nat,xat,yat,zat from nat1,xat1,..., if and only if
   * flag1 is 1
   */
  j=0;
  for(i = 1; i <= natoms; i++)
    if ( *(flag1 + i) ) {
      ++j;
      *(nat + j) = nat1[i];
      xat[j] = xat1[i];
      yat[j] = yat1[i];
      zat[j] = zat1[i];
      for(k=0; k<3; k++) fv[j][k] = fv1[i][k];
    }
  /* now, natoms is natoms - ndel */            
  natoms -= ndel;

  /* this is for forces -- TESTING */
  MallocForceVectors();
  BuildForceVectors(&FV);
  /* END */

  xcFree( (void*) xat1 );
  xcFree( (void*) yat1 );
  xcFree( (void*) zat1 );
  xcFree( (void*) nat1 );
  xcFree( (void*) flag1 );
  xcFree( (void*) fv1 );
  return XC_OK;
}



static int
ReadFrames( FILE *fp , long int framepos )
{
  int sizeIA, sizeDA; /* for size of malloc */
  int i, n;
  
  /* go to atoms position !!!!! */
  fseek(fp, framepos, 0);

  /* read NAT,X,Y,Z from file */
  /* allocate memory for NAT,X,Y,Z */ 
  sizeIA = sizeof(int) * (nframes + 1);
  sizeDA = sizeof(double) * (nframes + 1);
    
  /*sqn = (int *) malloc(sizeIA);*/
  frametype = (int *) malloc(sizeIA);
  xframe = (double *) malloc(sizeDA);
  yframe = (double *) malloc(sizeDA);
  zframe = (double *) malloc(sizeDA);
  xframe2 = (double *) malloc(sizeDA);
  yframe2 = (double *) malloc(sizeDA);
  zframe2 = (double *) malloc(sizeDA);

  /* read coordinates and nat */
  for(i = 1; i <= nframes; i++) 
    {
      n = fscanf(fp,"%d %lf %lf %lf %lf %lf %lf",frametype + i, 
		 xframe + i, yframe + i, zframe + i,
		 xframe2 + i, yframe2 + i, zframe2 + i);     
      if ( n < 7 ) {
	fprintf(stderr,"ERROR: Error reading file ");
	return XC_ERROR;
      }      
    }
  return XC_OK;
}


static int
ReadVoronoi(FILE *fp, float poly[][WIGNERSEITZ_MAXVERTEX][3], 
	    float norm[][WIGNERSEITZ_MAXVERTEX][3], int nvert[],
	    float *Max, boolean *xcrl, int *npoly)
{
  int eos, i, j, nv;
  float max = 0;
  /* double fabs(); */
  char line1[256];

  /* VORONOI section ends with XXXXXXXXXXXXX_END,
     so this will do */
  *npoly = 0;
  eos = 0;
  while( !eos ) { 	
    fscanf( fp, "%s", line1);
    if ( (nv = atoi(line1)) != 0 ) {
      /* t.k */
      for(i = 0; i < nv; i++) {
	if ( fscanf( fp, "%f %f %f %f %f %f",
		     &poly[*npoly][i][0],
		     &poly[*npoly][i][1],
		     &poly[*npoly][i][2],
		     &norm[*npoly][i][0],
		     &norm[*npoly][i][1],
		     &norm[*npoly][i][2] ) == 6 ) {
	  for(j = 0; j < 3; j++) {
	    if ( max < (float) ABS (poly[*npoly][i][j]) ) 
	      max = (float) ABS ( poly[*npoly][i][j]);
	  }
	  normalizepvf( &norm[*npoly][i][0],
			&norm[*npoly][i][1],
			&norm[*npoly][i][2] );
	} else {
	  fprintf(stderr, "ERROR: Error reading VORONOI-POLYHEDRA section, while reading ");
	  return XC_ERROR;
	}
      }
      nvert[*npoly] = nv;
      (*npoly)++;
    } else {    
      *Max = max;
      *xcrl = 1;
      eos = 1;
    }
  }

  return XC_OK;
}


void 
FindMaxRad(void)
{
  int i;
  /* which atom have maximum radius -- this information is needed for
   * determination of glOrtho
   */

  /* if (0) { */
  /*   /\* this is too complicated, because the viewport changes when radii change *\/ */
  /*   max.r = 0.0; */
  /*   for(i=0; i<atm.natomkind; i++) */
  /*     if (atrad[atm.atomkind[i]] > max.r) max.r = atrad[atm.atomkind[i]];  */
  
  /*   /\* this is too simple *\/ */
  /*   max.r = 1.5; */
  /* } */

  max.r = 0.0;
  for(i=0; i<atm.natomkind; i++) {
    if (rvdw[atm.atomkind[i]] > max.r) max.r = rvdw[atm.atomkind[i]];   
    if (rcov[atm.atomkind[i]] > max.r) max.r = rcov[atm.atomkind[i]];
  }
}


/* maybe we do not have crystal & some vectors that were just 
 * read are meaningless, coorect it by the size of the strucure
 */
static void
CorrectVectors( double vec[][4], int dim ) 
{
  if ( dim <= 2 ) {
    vec[2][0] = 0.0;
    vec[2][1] = 0.0;
    vec[2][2] = max.z;
    if (abs(max.z) < MINTOL) vec[2][2] = 1.0;
  }
  if ( dim <= 1 ) {
    vec[1][0] = 0.0;
    vec[1][1] = max.y;
    vec[1][2] = 0.0;
    if (abs(max.y) < MINTOL) vec[1][1] = 1.0;
  }
  if ( dim == 0 ) {
    vec[0][0] = max.x;
    vec[0][1] = 0.0;
    vec[0][2] = 0.0;
    if (abs(max.x) < MINTOL) vec[0][0] = 1.0;
  }
}


static int 
ReadPDB(FILE *fp) 
{
  int i, sizeIA, sizeDA;
  int *nat1, *flag1;
  double *xat1, *yat1, *zat1, (*fv1)[3];
  char line[256];

  natoms = 0;
  /* the file will have to be readed twice, first time just count # of atoms */
  while ( fgets(line, sizeof(line), fp) != NULL ) {
    if ( strncmp(line, "ATOM", 4) == 0 || strncmp(line, "HETATM", 6) == 0) {
      /* ATOM section */
      natoms++;
    }
  }

  rewind(fp);
  /* malloc here */
  sizeIA = sizeof(int) * (natoms + 1);
  sizeDA = sizeof(double) * (natoms + 1);
  
  nat1  = (int *)    malloc(sizeIA);
  flag1 = (int *)    malloc(sizeIA);
  xat1  = (double *) malloc(sizeDA);
  yat1  = (double *) malloc(sizeDA);
  zat1  = (double *) malloc(sizeDA);
  /* fv1 dummy so far !!! */
  fv1   = (double (*)[3]) malloc( sizeof(double [3]) );

  i=0;
  while ( fgets(line, sizeof(line), fp) != NULL ) {
    if ( strncmp(line, "ATOM", 4) == 0 ) {
      i++;
      *(flag1 + i) = 1;    
      parsePDBAtomRecord(line, nat1 + i, xat1 + i, yat1 + i, zat1 + i);
    }
    else if ( strncmp(line, "HETATM", 6) == 0 ) {
      i++;
      *(flag1 + i) = 1;    
      parsePDBHetatmRecord(line, nat1 + i, xat1 + i, yat1 + i, zat1 + i);
    }
  }

  CheckAtoms(flag1, nat1, xat1, yat1, zat1, fv1 );

  current_file_format = FORMAT_PDB;

  return XC_OK;
}


static int
ReadXYZ(FILE *fp) 
{
  int n, na, i, j, sizeIA, sizeDA;
  int *nat1, *flag1;
  double *xat1, *yat1, *zat1, qc, (*fv1)[3];
  char atmn[80];
  char line[256];

  /* read first line; number of atoms */
  if ( fgets(line, sizeof(line), fp) == NULL ) {
    fprintf(stderr, 
	    "ERROR unexpected end of file; while reading file ");
    return XC_ERROR;
  }
  sscanf(line, "%d", &na);
  if ( na < 1 ) {
    fprintf(stderr, "number of atoms lower then one !!!; ERROR reding file ");
    return XC_ERROR;
  } else {
    natoms = na;
  }

  /* read comment */
  if ( fgets(line, sizeof(line), fp) == NULL ) {
    fprintf(stderr, "ERROR while reading file; unexpected end of file ");
    return XC_ERROR;
  }

  /* malloc here */
  sizeIA = sizeof(int) * (natoms + 1);
  sizeDA = sizeof(double) * (natoms + 1);
    
  nat1  = (int *)    malloc(sizeIA);
  flag1 = (int *)    malloc(sizeIA);
  xat1  = (double *) malloc(sizeDA);
  yat1  = (double *) malloc(sizeDA);
  zat1  = (double *) malloc(sizeDA);
  fv1  = (double (*)[3]) malloc( sizeof(double [3]) * (natoms + 1) );

  for(i=1; i<=natoms; i++) {
    *(flag1 + i) = 1;    
    if ( fgets(line, sizeof(line), fp) == NULL ) {
      fprintf(stderr, 
	      "ERROR while reading file; unexpected end of file ");
      return XC_ERROR;
    }

    if ( (n = sscanf(line,"%s %lf %lf %lf %lf %lf %lf %lf\n",
		     atmn, &xat1[i], &yat1[i], &zat1[i], 
		     &qc, &fv1[i][0], &fv1[i][1], &fv1[i][2])) < 4 ) {
      fprintf(stderr, "ERROR while XYZ reading file !!!");
      xcFree((void*) nat1);
      xcFree((void*) xat1);
      xcFree((void*) yat1);
      xcFree((void*) zat1);
      xcFree((void*) flag1);
      return XC_ERROR;
    }
    if (n==7) xcr.lforce=1;
    
    /* parse atomic name */
    ParseAtomName(atmn);      

    /* find atomic number */
    nat1[i]=100; /* if atomic number will not be found, 100 will be used */
    for(j=0; j<=100; j++) {
      if ( strcmp(element[j], atmn) == 0 ) nat1[i] = j;
    }
  }
  CheckAtoms(flag1, nat1, xat1, yat1, zat1, fv1 );

  current_file_format = FORMAT_XYZ;

  return XC_OK;
}


/*
static void
StringTrailLeftWhiteSpaces(char *str, size_t size) {
  register int i, ii, left=1;
  char *name = (char *) malloc( sizeof(char) * size );

  strncpy(name, str, size);

  for (i=0, ii=0; i<size; i++) {
    if ( left && *(name + i) == ' ' ) {
      ;
    } else {
      left = 0;
      *(str + ii++) = *(name + i);
    }
  }
  xcFree((void*) name);
}
*/


void
ParseAtomName(char *atomName) {
  atomName[0] = toupper(atomName[0]);
  atomName[1] = tolower(atomName[1]);
  if ( isalpha( atomName[1] ) ) {
    atomName[2] = '\0';
  } else {
    atomName[1] = '\0';
  }
}


static void
CheckAtoms(int *flag1, int *nat1, double *xat1, double *yat1, double *zat1,
	   double (*fv1)[3]) 
{
  int i, j, ndel, sizeIA, sizeDA;
  double dis;

  /* if two atoms overlap, then we will have troubles when making bonds;
   * if two atoms are to close just disregard one of them
   */
  ndel=0;
  for(i = 1; i<= natoms - 1; i++)
    {
      if ( nat1[i] == 0 ) continue; /* dummy atom */
      if (*(flag1 + i)) {
	for(j = i + 1; j <= natoms; j++) {
	  if ( nat1[j] == 0 ) continue; /* dummy atom */
	  dis = dist6(xat1[i], *(xat1 + j),
		      yat1[i], *(yat1 + j),
		      zat1[i], *(zat1 + j));
	  if ( dis < MINDIS ) {
	    if ( dis < MINTOL ) 
	      fprintf(stderr,"WARNING: Atom %d and atom %d overlap !!!\n",i,j);
	    else 
	      fprintf(stderr,
		      "WARNING: Atom %d and atom %d are very close !!!\n",i,j);
	    fprintf(stderr,"         Atom %d deleted !!!\n",j);
	    ++ndel;
	    *(flag1 + j) = 0;
	  }
	}
      }
    }

  /* allocate memory for NAT,X,Y,Z */ 
  sizeIA = sizeof(int) * (natoms - ndel + 1);
  sizeDA = sizeof(double) * (natoms - ndel + 1);
  
  nat = (int *)    malloc(sizeIA);
  xat = (double *) malloc(sizeDA);
  yat = (double *) malloc(sizeDA);
  zat = (double *) malloc(sizeDA);
  fv  = (double (*)[3]) malloc( sizeof(double [3]) * (natoms + 1) );
  xcReallocAtomLabels();  
  
  /* write nat,xat,yat,zat from nat1,xat1,..., if and only if
   * flag1 is 1
   */
  j=0;
  for(i = 1; i <= natoms; i++)
    if ( *(flag1 + i) ) {
      ++j;
      *(nat + j) = nat1[i];
      xat[j] = xat1[i];
      yat[j] = yat1[i];
      zat[j] = zat1[i];
      if (xcr.lforce) {
	fv[j][0] = fv1[i][0];
	fv[j][1] = fv1[i][1];
	fv[j][2] = fv1[i][2];
      }
    }
  /* now, natoms is natoms - ndel */            
  natoms -= ndel;

  xcFree((void*) xat1 );
  xcFree((void*) yat1 );
  xcFree((void*) zat1 );
  xcFree((void*) nat1 );
  xcFree((void*) flag1 );
  xcFree((void*) fv1 );
}

static void
parsePDBAtomRecord(char *line, int *n, double *x, double *y, double *z) 
{
  /* ATOMFORMAT "(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,1x,i3)" */
  /*                       12                      31,39,47            */
  /* "%*6s%5i%*1c%4s%*15c%8.3f%8.3f%8.3f" */
  register int j;
  char coor[10];
  char atmn[3]; 
  /* atom name */
  /*strncmp(atmn, line + 11, 4);*/
  /*strncmp(atmn, line + 11, 2);*/

  sscanf(line, "%*6s%*5i%*1c%1c%1c", &atmn[0], &atmn[1]);
  atmn[2] = '\0';

  if ( atmn[0] == ' ' ) {
    atmn[0] = atmn[1];
    atmn[1] = '\0';
  }
  if (atmn[0] == 'A') 
    {
      if (    toupper(atmn[1]) != 'L'    /* Al */
	   && toupper(atmn[1]) != 'R'    /* Ar */
	   && toupper(atmn[1]) != 'S'    /* As */
	   && toupper(atmn[1]) != 'G'    /* Ag */
	   && toupper(atmn[1]) != 'U'    /* Au */
	   && toupper(atmn[1]) != 'T'    /* At */
	   && toupper(atmn[1]) != 'C'    /* Ac */
	   && toupper(atmn[1]) != 'M'    /* Am */
	   ) {
	atmn[0] = 'X';
	atmn[1] = '\0';
      }
    } 
  else if (atmn[0] == 'D') {
    atmn[0] = 'H';
    atmn[1] = '\0';
  }
  /*fprintf(stderr,"PDB>  ATOM = %s\n", atmn);*/

  /*StringTrailLeftWhiteSpaces(atmn,2);*/
  /*ParseAtomName(atmn);*/

  /* find atomic number */
  *n=100; /* if atomic number will not be found, 100 will be used */
  for(j=0; j<=100; j++) {
    if ( strcmp(element[j], atmn) == 0 ) {
      *n = j;
      break;
    }
  }

  /* X coor */
  strncpy(coor, line + 30, 8);
  sscanf(coor, "%lf", x);

  /* X coor */
  strncpy(coor, line + 38, 8);
  sscanf(coor, "%lf", y);

  /* X coor */
  strncpy(coor, line + 46, 8);
  sscanf(coor, "%lf", z);
  
  /*fprintf(stderr,"PDB>       Coords = %f %f %f\n", *x, *y, *z);*/
}



/* ------------------------------------------------------------------------ */

static void
parsePDBHetatmRecord(char *line, int *n, double *x, double *y, double *z) 
{
  register int j;
  char coor[10];
  char atmn[3], resn[4];
  char _atmn[3], *__atmn;

  __atmn = _atmn;

  sscanf(line, "%*6s%*5i%*1c%1c%1c%*3s%3s", &__atmn[0], &__atmn[1], resn);
  __atmn[2] = '\0';
  if ( strncmp(resn, "   ", 3) == 0 ) strncpy(resn, "UNK\0", 4);

  /*
    
  This is the BABEL recipe:

  * 1. if there is no residue name then we assume the atom label is
  * the name of the element.

  * 2. if the residue name and the atom label are the same then we
  * will assume that the atom label is the name of the element.


  * 3. generally, if there is a letter in the first column then the
  * atom label represents an element with a 2 letter element name.
  * The exceptions so far to this is when the residue name is
  * ADR, COA, FAD, GPG, NAD, NAL or NDP.  In these cases there is a
  * letter in the first column and the second letter is the element
  * name.
  
  */  

  if ( !strncmp ( resn, "ADR", 3 ) || !strncmp ( resn, "COA", 3 ) ||
       !strncmp ( resn, "FAD", 3 ) || !strncmp ( resn, "GPG", 3 ) ||
       !strncmp ( resn, "NAD", 3 ) || !strncmp ( resn, "NAL", 3 ) ||
       !strncmp ( resn, "NDP", 3 ) ) {
    __atmn++;
    __atmn[1] = '\0';
  } else if (isdigit(__atmn[0])) {
    __atmn++;
    __atmn[1] = '\0';   
  } else if (__atmn[0] != ' ') {
    __atmn[2] = '\0';
    if (isalpha(__atmn[1]) && isupper(__atmn[1])) __atmn[1] = tolower(__atmn[1]);
  } else {
    if ( __atmn[0] == ' ' ) __atmn[0] = __atmn[1];
    __atmn[1] = '\0';
  }
  /*fprintf(stderr,"PDB>  HETATM = %s\n", __atmn);*/

  strcpy(atmn,__atmn);

  /* find atomic number */
  *n=100; /* if atomic number will not be found, 100 will be used */
  for(j=0; j<=100; j++) {
    if ( strcmp(element[j], atmn) == 0 ) {
      *n = j;
      break;
    }
  }

  /* X coor */
  strncpy(coor, line + 30, 8);
  sscanf(coor, "%lf", x);

  /* X coor */
  strncpy(coor, line + 38, 8);
  sscanf(coor, "%lf", y);

  /* X coor */
  strncpy(coor, line + 46, 8);
  sscanf(coor, "%lf", z);

  /*fprintf(stderr,"PDB>       Coords = %f %f %f\n", *x, *y, *z);*/
}


static void xcReallocAtomLabels(void) {
  static int natoms_old = 0;
  
  if (natoms_old == 0) {
    natoms_old = natoms;
  }
  if ( !atomLabel) {
    atomLabel = (AtomicLabel*) xcCalloc((natoms + 1), sizeof(AtomicLabel));
    do_not_display_atomlabel = (short*) xcCalloc((natoms + 1), sizeof(short));
  } 

  else if (natoms_old == 0) {
    /* something wrong */
    xcError("     Something wrong in xcReallocAtomLabels:\n     atomLabel already aloccatedm but natoms_old == 0");
  }
  
  else if ( natoms > natoms_old ){
    /* NOTE: if natoms <= natoms_old --> do nothing, there is enought
       of space, therefore here we only consider natoms > natoms_old !!! */
    AtomicLabel *_labels = (AtomicLabel*) xcCalloc((natoms + 1), sizeof(AtomicLabel));
    short       *_do_not = (short*)       xcCalloc((natoms + 1), sizeof(short));
    int         i;

    for (i=1; i<=natoms_old; i++) {
      _labels[i].base        = atomLabel[i].base;               
      _labels[i].width	     = atomLabel[i].width;	      
      _labels[i].height      = atomLabel[i].height;	      
      _labels[i].do_display  = atomLabel[i].do_display;	      
      _labels[i].label	     = atomLabel[i].label;	      
      _labels[i].tkfont      = atomLabel[i].tkfont;           
      
      COPY_V(3,_labels[i].bright_color, atomLabel[i].bright_color);  
      COPY_V(3,_labels[i].dark_color,   atomLabel[i].dark_color);    
      
      _do_not[i] = do_not_display_atomlabel[i];
    }

    xcFree((void*)atomLabel);
    xcFree((void*)do_not_display_atomlabel);
    
    atomLabel = _labels;
    do_not_display_atomlabel = _do_not;
    
    /***************/
    natoms_old = natoms;
    /***************/
  }
}


static void ToEndLine(FILE *fp)
{
  int c;
  while ((c = getc(fp)) != '\n' && c != EOF) { ; }
}


