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
 * Source: $XCRYSDEN_TOPDIR/C/hbonds.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "xcfunc.h"
#include "memory.h"

/* 
   variable
*/

H_Bond hbonds = {
  0,                        /* number of H-bonds                             */
  10,                       /* max_n                                         */
  (int*)NULL,               /* H_like_list: H-only                           */
  (int*)NULL,               /* O_like_list: N, O, F, Cl                      */
  {1.0, 1.0, 1.0, 1.0},     /* color of the H-bonds                          */
  1.3,                      /* minimum H-distance                            */
  2.5,                      /* maximum H-distance                            */
  120.0,                    /* minimum angle                                 */
  0xAAAA,                   /* line pattern                                  */
  2,                        /* pattern size (i.e. how long pattern appears)  */
  4.0,                      /* line width                                    */
  (VEC3d*) NULL,            /* start_coor                                    */
  (VEC3d*) NULL             /* end_coor                                      */
};
extern GLuint tempDisable3Dlist, tempEnable3Dlist;

/*
  function prototypes
*/
static short int _atomFromList(int nat, int *list);


/*
  global functions
*/

void
make_H_Bonds(void)
{
  register int ih, _ih, jo;
  int n;
  double dis, len, angle, chemBondVec[3], HbondVec[3];

  if (hbonds.start_coor == (VEC3d*) NULL) 
    {
      hbonds.start_coor = (VEC3d*) xcCalloc(hbonds.max_n, sizeof(VEC3d));
      hbonds.end_coor   = (VEC3d*) xcCalloc(hbonds.max_n, sizeof(VEC3d));
    }

  /* 
     NOTE: 1st loop(ih) is over H-like atoms, 
     2nd loop(jo) is over O-like atoms 
  */ 

  n = 0;
  for (ih = 1; ih <= natoms; ih++) 
    {
      /* check if ith-atom is Hydrogen like */

      if ( _atomFromList(nat[ih], hbonds.H_like_list ) == XC_FALSE) continue;
      
      /* NOTE: 
	 -----
         H-bond is formed only when X-H....X angle is close to 180
	 degrees. For this reason we have to find the chemical bonds
	 of the ih-atom and see if a potential H-bond forms the
	 angle close to 180 degree-->only then we assign it as an H-bond
      */
      for (_ih = 1; _ih <= natoms; _ih++) 
	{
	  if (ih == _ih) continue;
	  dis = dist6(xat[ih], xat[_ih],
		      yat[ih], yat[_ih],
		      zat[ih], zat[_ih]);
	  len = rcov[nat[ih]] + rcov[nat[_ih]];
	  
	  if ( dis < len ) 
	    {
	      /* chemical bond */
	      chemBondVec[0] = xat[_ih] - xat[ih];
	      chemBondVec[1] = yat[_ih] - yat[ih];
	      chemBondVec[2] = zat[_ih] - zat[ih];
	      /* hydrogen forms only one chemicla bond */
	      break;
	    }
	}


      for (jo = 1; jo <= natoms; jo++) 
	{
	  if (ih == jo) continue;
	  
	  /* check if jth-atom is O-like */
	  
	  if ( _atomFromList(nat[jo], hbonds.O_like_list ) == XC_FALSE) continue;
	  
	  /* check if H-bond distance criteria is meet */
	  
	  dis = dist6(xat[ih],xat[jo],  yat[ih],yat[jo],  zat[ih],zat[jo]);

	  if ( dis < hbonds.length_max && dis > hbonds.length_min ) {	    
	    /* 
	       The distance criteria for hydrogen bond is met !!!
	       Check the angle criteria
	    */
	    HbondVec[0] = xat[jo] - xat[ih];
	    HbondVec[1] = yat[jo] - yat[ih];
	    HbondVec[2] = zat[jo] - zat[ih];
	    angle = RAD2DEG * acos(DOT_V3(HbondVec, chemBondVec) / 
				   (NORM_V3(HbondVec) * NORM_V3(chemBondVec)));
	    
	    if (ABS(angle) > hbonds.angle_min) 
	      {
		/* we have a hydrogen bond */
		
		hbonds.start_coor[n][0] = xat[ih];
		hbonds.start_coor[n][1] = yat[ih];
		hbonds.start_coor[n][2] = zat[ih];
		
		hbonds.end_coor[n][0]   = xat[jo];
		hbonds.end_coor[n][1]   = yat[jo];
		hbonds.end_coor[n][2]   = zat[jo];	    
		n++;
	      }

	    if (n >= hbonds.max_n) 
	      {
		/* reallocate */
		hbonds.max_n *= 2;
		hbonds.start_coor = (VEC3d*)
		  xcRealloc(hbonds.start_coor, hbonds.max_n * sizeof(VEC3d));
		hbonds.end_coor   = (VEC3d*)
		  xcRealloc(hbonds.end_coor, hbonds.max_n * sizeof(VEC3d));	      
		/*
		  hbonds.chemBond   = (double (*)[3]) 
		  xcRealloc(hbonds.chemBond, 
		  (size_t) hbonds.max_n * sizeof(double [3]));
		*/
		/*
		  hbonds.start_coor = xcReallocVEC3d(hbonds.start_coor, hbonds.max_n, hbonds.max_n*2);
		  hbonds.end_coor   = xcReallocVEC3d(hbonds.end_coor,   hbonds.max_n, hbonds.max_n*2);
		*/
	      }
	  }
	}
    }

  hbonds.n = n;
}


void
xcRenderHbonds3D(void)
{
  register int ib;
  
  glCallList    ( tempDisable3Dlist );
  glEnable      ( GL_LINE_STIPPLE   );
  glLineStipple ( hbonds.line_patternsize, hbonds.line_pattern );
  glColor3fv    ( hbonds.color      );
  glLineWidth   ( hbonds.line_width   );

  glBegin(GL_LINES);
    for (ib=0; ib < hbonds.n; ib++) {
      glVertex3dv( hbonds.start_coor[ib] );
      glVertex3dv( hbonds.end_coor[ib]   );
    }
  glEnd();

  glDisable  ( GL_LINE_STIPPLE  );   
  glCallList ( tempEnable3Dlist );
}



/*
  static functions
*/

static short int _atomFromList(int nat, int *list) {
  while (*list > 0) {
    if (nat == *list) {
      return XC_TRUE;
    }
    list++;
  }
  return XC_FALSE;
}


/*
VEC3d *xcReallocVEC3d(VEC3d *v, int old_nmemb, int new_nmemb) {
  register int i;
  VEC3d *ptr = (VEC3d*) xcCalloc(new_nmemb, sizeof(VEC3d));
  fprintf(stderr,"xcReallocVEC3d vec = %d\n", v);
  for (i=0; i<old_nmemb; i++) 
    {
      ptr[i][0] = v[i][0];
      ptr[i][1] = v[i][1];
      ptr[i][2] = v[i][2];
    }
  fprintf(stderr,"xcReallocVEC3d vec = %d\n", v);
  free(v);
  return ptr;
}
*/
