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
 * Source: $XCRYSDEN_TOPDIR/C/detnsplit.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "struct.h"
#include "memory.h"

extern double rcov[MAXNAT + 1];
/* functions prototypes */
extern double dist6(double x1, double x2, 
		    double y1, double y2, double z1, double z2);
extern double dist3(double x, double y, double z);
extern int normalizepv(double *x, double *y, double *z);

int
DetNBonds(int *nsplit)
{
  int ibn, i, j, ii, jj; /* nn is increased before usage */
  int numb = 0; /* numb is ++numb before usage */
  double (*xb)[MAX_BONDS], (*yb)[MAX_BONDS],  (*zb)[MAX_BONDS];
  double (*xb2)[MAX_BONDS], (*yb2)[MAX_BONDS],  (*zb2)[MAX_BONDS];
  double angl, minangl = PI; /* parameters to detemine min angle */
  double dis, len, fdis, 
    minbdis = 10.0, maxfdis = 0.0; /* bond lenght parameters */
  double dx1, dx2, dy1, dy2, dz1, dz2; 

  xb = (double (*)[MAX_BONDS]) malloc( sizeof(double [MAX_BONDS]) * (natoms + 1) );
  yb = (double (*)[MAX_BONDS]) malloc( sizeof(double [MAX_BONDS]) * (natoms + 1) );
  zb = (double (*)[MAX_BONDS]) malloc( sizeof(double [MAX_BONDS]) * (natoms + 1) );
  xb2 = (double (*)[MAX_BONDS]) malloc( sizeof(double [MAX_BONDS]) * (natoms + 1) );
  yb2 = (double (*)[MAX_BONDS]) malloc( sizeof(double [MAX_BONDS]) * (natoms + 1) );
  zb2 = (double (*)[MAX_BONDS]) malloc( sizeof(double [MAX_BONDS]) * (natoms + 1) );
  
  /* make bonds & also determine if additional spliting is needed */
  for(i = 1; i <= natoms; i++) {
    /* 
       don't make a bond for dummy X-atom 
    */
    if (nat[i] == 0) continue;
    
    ibn = 0; /* if ibn > 0 -> measure bond angle */
    for(j = 1; j <= natoms; j++)
      {		
	/* 
	   don't make a bond for dummy X-atom 
	*/
	if (nat[j] == 0) continue;
 
	if ( j == i ) continue; /* do not compare atom with it's self */

	/* calculate distance between two atoms */
	dis = dist6(*(xat + i), *(xat + j),
		    *(yat + i), *(yat + j),
		    *(zat + i), *(zat + j));
	if ( dis < MINDIS ) continue;	
	/* len == sum of covalent radius from atom i & j */
	len = rcov[*(nat + i)] + rcov[*(nat + j)];
	if ( dis < len )
	  {
	    ++ibn; /* criteria if bond angle is measured */	    
	    if ( ibn == (int) MAX_BONDS/2.0 ) {
	      fprintf(stderr,
		      "WARNING: atom \"%d\" has huge number of bonds !!!\n",
		      i);
	      /* continue; */
	    }
	    /* each bond belong to both atoms, so we'll split/devide 
	       a bond to both atoms, criteria for spliting is fdis 
	       parameter */

	    /* only BONDs to ATOM i */
	    ++numb;

	    /* FIRST END OF A BOND */
	    fdis = dis * rcov[*(nat + i)] / len;
	    if ( dis < minbdis ) minbdis = dis; /* needed for additional  */
	    if ( fdis > maxfdis ) maxfdis = fdis; /* bond spliting criteria */
	    xb[i][ibn] = *(xat + i);
	    yb[i][ibn] = *(yat + i);
	    zb[i][ibn] = *(zat + i);
	    /* SECOND END OF A BOND */
	    xb2[i][ibn] = *(xat + i) + fdis * (*(xat + j) - *(xat + i)) / dis;
	    yb2[i][ibn] = *(yat + i) + fdis * (*(yat + j) - *(yat + i)) / dis;
	    zb2[i][ibn] = *(zat + i) + fdis * (*(zat + j) - *(zat + i)) / dis;
	  } /* if (dis < len) */

	if ( ibn >= MAX_BONDS - 1 ) {
	  fprintf(stderr, "WARNING: atom \"%d\" has more than %d bonds. Further bonds will be neglected !!!\n", i, MAX_BONDS);
	  /* atom i has too many bonds, stop at this point */
	  break;
	}
      } /* for(j=... */

    /* now measure bond angles between bonds of atom i */
    for(ii = 1; ii <= ibn - 1; ii++)
      for(jj = ii + 1; jj <= ibn; jj++)
	{
	  dx1 = xb2[i][ii] - xb[i][ii];
	  dy1 = yb2[i][ii] - yb[i][ii]; 
	  dz1 = zb2[i][ii] - zb[i][ii];
	  dx2 = xb2[i][jj] - xb[i][jj];
	  dy2 = yb2[i][jj] - yb[i][jj];
	  dz2 = zb2[i][jj] - zb[i][jj];

	  if (!normalizepv(&dx1, &dy1, &dz1)) continue;
	  if (!normalizepv(&dx2, &dy2, &dz2)) continue;

	  angl = acos ( dx1 * dx2 + dy1 * dy2 + dz1 * dz2 );

	  if ( angl < minangl )
	    {
	      minangl = angl;
	      if ( angl < 0.087) /* this is appr. 5 degrees */
		{
		  minangl = 0.53; /* this is appr. 30 degrees */
		  /* fprintf(stderr,"WARNING: Bond angle (%f degrees) between atom %d and atom %d is very small !!!\n", angl, i, j); */
		}		  
	    }
	} 
    
  } /* for(i=... */

  /* DETERMINE IF ADDITIONAL SPLITING OF BONDS is NECESSARY */
  if ( minangl < 0.53 ) {
    /* set minangl to cca. 30 degrees */
    minangl = 0.53;
  }
  *nsplit = (int) (maxfdis * cos(minangl)) / (2.0 * minbdis) + 1;
  if ( *nsplit <= 0 ) *nsplit = 1;

  /* now I now how many bonds is there */
  nbonds = numb * *nsplit;
  fprintf(stderr,"Estimated number of bonds = %d\n", nbonds);
  fflush(stderr);

  xcFree((void*)zb2);
  xcFree((void*)yb2);
  xcFree((void*)xb2);
  xcFree((void*)zb);
  xcFree((void*)yb);
  xcFree((void*)xb);

  return XC_OK;
}
