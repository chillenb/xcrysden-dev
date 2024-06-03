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
 * Source: $XCRYSDEN_TOPDIR/C/xsf2xsf.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MIN(x,y) ( (x)<(y) ? (x) : (y) )


#define LINE_FMT   "%[^\n\r]"


char *element[] = {
  "X",
  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O", 
  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S", 
  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", 
  "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", 
  "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", 
  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", 
  "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", 
  "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", 
  "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", 
  "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", 
  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", 
  "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", 
  "Bk", "Cf", "Es", "Fm"};


/* xsf2xsf.c */
int main(int argc, char *argv[]);
int xsf_readFile(char *inpfile, char *outfile, int reduceTo);
void xsf_dummyReadCoor(FILE *inpFP, int natoms);
void xsf_readPrintCoor(char *line, FILE *outFP);
void PrintVec(FILE *outFP, double vec[][4]);
int ReadVec(FILE *fp, double vec[][4]);
void ParseAtomName(char *atomName);
void EndLine(FILE *inpFP, FILE *outFP);
void ToEndLine(FILE *inpFP);

/*
  Usage: xsf2xsf inputFile outputFile [reduceTo]
*/
int
main(int argc, char *argv[])  
{
  int reduceTo = 3;

  if ( argc != 3 && argc != 4 ) {
    fprintf(stderr,"\nUsage: xsf2xsf inputFile outputFile [reduceTo]\n\n");
    exit(1);
  }

  if (argc == 4) {
    if ( ! isdigit ((int)argv[3][0]) ) {
      fprintf(stderr,"\nUsage: xsf2xsf inputFile outputFile [reduceTo]\n\n");
      exit(1);
    }
    reduceTo = (int) strtol(argv[3], NULL, 10);
  }
    
  xsf_readFile(argv[1], argv[2], reduceTo);
  /* xsf_readFile(stdin, stdout, reduceTo); */
  exit(0);
}

/* ========================================================================= */
/*                          READ XCRYSDEN STRUCTURE FILE                     */
/* ========================================================================= */
int
xsf_readFile(char *inpfile, char *outfile, int reduceTo)
{
  FILE *inpFP, *outFP;
  char *line = (char *) malloc( sizeof(char) * 20000 ); /* max. line length */
  char *_line = (char *) malloc( sizeof(char) * 20000 ); /* max. line length */
  char *word = (char *) malloc( sizeof(char) * 256 ); /* max. word length */
  char *unit = (char *) malloc( sizeof(char) * 256 ); /* max. word length */
  int  i, n, c, latom = 0;

  int dim = 0, groupn = 1, natr, ncell = 1, idummy;
  double prim[4][4], conv[4][4];

#ifdef WIN32
  inpFP = fopen(inpfile, "r");
  outFP = fopen(outfile, "w");
#else
  inpFP = fopen(inpfile, "r");
  outFP = fopen(outfile, "w");
#endif
  if ( inpFP == NULL || outFP == NULL) {
    fprintf(stderr,"error opening files: %s and %s\n", inpfile, outfile);
    exit(1);
  }
  
  while((c = fscanf(inpFP, LINE_FMT, line)) != EOF) {

    /*fprintf(stderr,"|%s| len=%d\n", line, strlen(line));*/

    word[0] = '\0';
    if (strlen(line) != 0) {
      sscanf(line, "%s %s", word, unit);
    }
    
    /* check for comment-line "#" */
    if ( strncmp(word, "#", 1) == 0 ) {
      /* ignore this line */
      ToEndLine(inpFP);
      line[0] = '\0';
      continue;
    }
    
    /* if INFO */
    if ( strncmp(word,"ANIMS",5) == 0 
	 || strncmp(word,"INFO",4) == 0 
	 || strncmp(word, "BEGIN_", 6) == 0 
	 || strncmp(word, "END_",4) == 0 
	 || strncmp(word, "RECIP-", 6) == 0 
	 || strncmp(word, "WIGNER-", 7) == 0 
	 || strncmp(word, "BRILLOUIN-",10) == 0
	 || strncmp(word, "FRAM",4) == 0
	 || strncmp(word, "DATAGRID", 8) == 0 
	 || strncmp(word, "BANDGRID", 8) == 0 ) {
      latom = 0;
      fprintf(outFP,"%s\n", line);
    } 

    /* if DIM-GROUP */
    else if ( strncmp(word,"DIM",3) == 0 ) {
      latom = 0; 
      /* read dimensionality and group of structure */
      ToEndLine(inpFP);
      fscanf(inpFP, LINE_FMT, _line);
      n = sscanf(_line,"%d %d", &dim, &groupn);
      if ( n < 1 ) {
	fprintf(stderr,
		"ERROR: Error reading DIM-GROUP section ");
	exit(1);
      }
      if ( dim < 0 || dim > 3 ) {
	fprintf(stderr,
		"ERROR: DIM is out of range ");
	exit(1);
      }
      dim = MIN(reduceTo, dim);
      if ( reduceTo > 0 ) {
	if ( reduceTo != 3 ) groupn = 1;
	fprintf(outFP, "DIM-GROUP\n%d %d\n", dim, groupn);
      }
    }  	
    
    /* MOLECULE */
    else if ( strncmp(word,"MOLECULE",8) == 0 ) {
      /* do nothing */
      reduceTo = 0;
      latom = 0; 
    }
	
    /* POLYMER */
    else if ( strncmp(word,"POLYMER",7) == 0 ) {
      latom = 0; 
      if (reduceTo >= 1) fprintf(outFP, "POLYMER\n");
    }

    /* SLAB */
    else if ( strncmp(word,"SLAB",4) == 0 ) {
      latom = 0; 
      if (reduceTo == 1 ) fprintf(outFP, "POLYMER\n");
      else if (reduceTo >= 2) fprintf(outFP, "SLAB\n");
    }
    
    /* CRYSTAL */
    else if ( strncmp(word,"CRYSTAL",7) == 0 ) {
      latom = 0; 
      if (reduceTo == 1 ) fprintf(outFP, "POLYMER\n");
      else if (reduceTo == 2) fprintf(outFP, "SLAB\n");
      else if (reduceTo != 0) fprintf(outFP, "CRYSTAL\n");
    }

    /* if PRIMVEC */
    else if ( strncmp(word,"PRIMV",5) == 0 ) {
      latom = 0;
      ToEndLine(inpFP); /* t.k. */
      if ( !ReadVec( inpFP, prim ) ) {
	fprintf(stderr,
		"ERROR: Error reading PRIMVEC section ");
	exit(1);
      }
      if (reduceTo >= 1) {
	if ( strncmp(unit,"angs",4) == 0 )
	  /* strip "angstrom" unit (if we have xsf2) */
	  fprintf(outFP, "%s\n", word);
	else
	  fprintf(outFP, "%s\n", line);
	PrintVec (outFP, prim);
      }
    }

    /* if CONVVEC */
    else if ( strncmp(word,"CONVV",5) == 0 ) {
      latom = 0; 
      ToEndLine(inpFP); /* t.k. */
       if ( !ReadVec( inpFP, conv ) ) {	
	fprintf(stderr,
		"ERROR: Error reading CONVVEC section ");
	exit(1);
      }
      if (reduceTo > 2) {
	/* Maybe it is safer to print the CONVVEC only for crystals */
	if ( strncmp(unit,"angs",4) == 0)
	  /* strip "angstrom" unit (if we have xsf2) */
	  fprintf(outFP, "%s\n", word);
	else
	  fprintf(outFP, "%s\n", line);
	PrintVec (outFP, conv);
      }
    }

    /* if PRIMCOORD */
    else if ( strncmp(word,"PRIMC",5) == 0 ) {
      latom = 0; 
      ToEndLine(inpFP); 
      fscanf(inpFP, LINE_FMT, _line);
      n = sscanf(_line,"%d %d", &natr, &idummy);
      if ( n < 1 ) {
	fprintf(stderr,
		"ERROR: Error reading PRIMCOORD section ");
	exit(1);
      }
      if (reduceTo == 0) {
	int ind;
	n = sscanf(line,"%s %d", word, &ind);
	if ( n==1 ) fprintf(outFP,"ATOMS\n");
	else if ( n==2 ) fprintf(outFP,"ATOMS %d\n", ind);
	else {
	  fprintf(stderr,
		  "ERROR: Error reading PRIMCOORD section ");
	  exit(1);
	}
      } else {
	if ( strncmp(unit,"angs",4) == 0)
	  /* strip "angstrom" unit (if we have xsf2) */
	  fprintf(outFP, "%s\n%d 1\n", word, natr);
	else
	  fprintf(outFP, "%s\n%d 1\n", line, natr);
      }
      for (i=0; i<natr; i++) {
	EndLine(inpFP, outFP);
	fscanf(inpFP, LINE_FMT, line);
	xsf_readPrintCoor (line, outFP);
      }
    }
    
    /* if CONVCOORD */
    else if ( strncmp(word,"CONVC",5) == 0 ) {
      latom = 0;
      EndLine(inpFP, outFP); 
      fscanf(inpFP, LINE_FMT, _line);
      n = sscanf(_line,"%d %d", &natr, &ncell);
      if ( n < 1 ) {
	fprintf(stderr,
		"ERROR: Error reading CONVCOORD section ");
	exit(1);
      }
      /* ignore this entries; do not print ... */
      xsf_dummyReadCoor (inpFP, natr*ncell);
    }    
      
    /* if ATOMS */
    else if ( strncmp(word,"ATOM",4) == 0 ) {
      latom = 1;
      fprintf(outFP, "%s\n", line);
    }
    
    else {
      if (latom == 0) {
	if (line[0] != '\r') fprintf(outFP, "%s\n", line);
      } else {
	/* check if not emppty line */
	if ( strlen(word) > 0 ) {
	  xsf_readPrintCoor (line, outFP);
	}
      }
    }
    EndLine(inpFP, outFP);
    line[0] = '\0';
  }

  free(unit);
  free(word);
  free(_line);
  free(line);

  fclose(inpFP);
  fclose(outFP);
  exit(0);
}


void xsf_dummyReadCoor(FILE *inpFP, int natoms) {
  int i;
  char *atm  = (char *) malloc( sizeof(char) * 10 );
  char *line = (char *) malloc( sizeof(char) * 256 ); /* 256 is max. width */
  double x, y, z, fx, fy, fz;
  char c;
  
  for (i=0; i<natoms; i++) {
    while((c = getc(inpFP)) != '\n' && c != EOF) {;}
    fscanf(inpFP, LINE_FMT, line);
    sscanf(line, "%s %lf %lf %lf %lf %lf %lf", atm, &x, &y, &z, &fx, &fy, &fz);  
  }
  free(line);
  free(atm);
}

void xsf_readPrintCoor (char *line, FILE *outFP) {
  int j, n, _nat;
  char *atm = (char *) malloc( sizeof(char) * 10 );
  double x, y, z, fx, fy, fz;
  
  n = sscanf(line, "%s %lf %lf %lf  %lf %lf %lf", 
	     atm, &x, &y, &z, &fx, &fy, &fz);

  if ( isdigit(atm[0]) ) {
    if ( n==4 )  
      fprintf(outFP, "%s %14.9f %14.9f %14.9f\n", atm, x, y, z);
    else if ( n==7 ) {
      fprintf(outFP, "%s %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f\n",
	      atm, x, y, z, fx, fy, fz);
    } else {
      fprintf(stderr, "ERROR: wrong number of fields (%d) in atomic coordinates, should be 4 or 7", n);
      exit(1);
    }
  } else {
    /* convert atomic-symbol to atomic-number */
    /* find atomic number */
    ParseAtomName(atm);
    
    _nat = 0; /* if atomic number will not be found, 0 will be used */
    for(j=0; j<100; j++) {
      if ( strncmp(element[j], atm, 2) == 0 ) {
	_nat = j;
	break;
      }
      if (strlen(element[j]) == 1) {
	if (strncmp(element[j], atm, 1) == 0 ) _nat = j;
      }
    }
    if ( n==4 )  
      fprintf(outFP, "%3d  %14.9f %14.9f %14.9f\n", _nat, x, y, z);
    else
      fprintf(outFP, "%3d  %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f\n",
	      _nat, x, y, z, fx, fy, fz);      
  }
  free(atm);
}
			
		
void 
PrintVec(FILE *outFP, double vec[][4]) {
  fprintf(outFP, "%14.9f %14.9f %14.9f\n", 
	  vec[0][0], vec[0][1], vec[0][2]);
  fprintf(outFP, "%14.9f %14.9f %14.9f\n", 
	  vec[1][0], vec[1][1], vec[1][2]);
  fprintf(outFP, "%14.9f %14.9f %14.9f\n", 
	  vec[2][0], vec[2][1], vec[2][2]);
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
      
  if ( n0 == 3 && n1 == 3 && n2 == 3 ) return 1;
  
  else return 0;
}

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

/*--->*/
void EndLine(FILE *inpFP, FILE *outFP) {
  int c;
  while((c = getc(inpFP)) != '\n' && c != EOF) {
    if (c != '\r' ) fprintf(outFP,"%c", c);
  }
  /*fprintf(outFP, "\n");*/
}
void ToEndLine(FILE *inpFP)
{
  int c;
  while ((c = getc(inpFP)) != '\n' && c != EOF) { ; }
}
