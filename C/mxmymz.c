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
 * Source: $XCRYSDEN_TOPDIR/C/mxmymz.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <math.h>
#include <stdio.h>
#include "struct.h"
#include "xcfunc.h"

void Set_mx_my_mz(double sumx, double sumy, double sumz);

void Set_mx_my_mz(double sumx, double sumy, double sumz) {
  int i, j;
  double _vec[3][3];

  /* 
     ------------------------------------------------------------------------
     if we have a molecule -- the center is the mass center
     if we have a crystal  -- the center is the cell center
     ------------------------------------------------------------------------
  */
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      _vec[i][j] = 0.0;
      if ( xcr.celltype == XCR_CONVCELL ) {
	_vec[i][j] = vec.conv[i][j];
      } else {
	_vec[i][j] = vec.prim[i][j];
      }
    }
  }
  
  if (xcr.dim == 0) {

    mx_my_mz.mx = sumx / ((double) natoms);
    mx_my_mz.my = sumy / ((double) natoms);
    mx_my_mz.mz = sumz / ((double) natoms);
    
  } else if (xcr.dim == 1) {
    
    mx_my_mz.mx = 0.5 * (double)xcr.nunit[0] * _vec[0][0];
    mx_my_mz.my = sumy / ((double) natoms);
    mx_my_mz.mz = sumz / ((double) natoms);

  } else if (xcr.dim == 2) {
    
    /*
      here we suppose that non-periodic dimension goes along Z axis !!!
      (if this is not fulfilled we are in trouble !!!)
    */
    /*
    fprintf(stderr,"---------------:\n");
    fprintf(stderr,"Lattice vectors:\n");
    fprintf(stderr,"---------------:\n");
    fprintf(stderr,"  %12.8f %12.8f %12.8f\n", 
	    _vec[0][0], _vec[0][1], _vec[0][2]);
    fprintf(stderr,"  %12.8f %12.8f %12.8f\n", 
	    _vec[1][0], _vec[1][1], _vec[1][2]);
    fprintf(stderr,"  %12.8f %12.8f %12.8f\n", 
	    _vec[2][0], _vec[2][1], _vec[2][1]);
    fprintf(stderr,"  %3d %3d %3d\n", 
	    xcr.nunit[0], xcr.nunit[1], xcr.nunit[2]);
    */
    
    mx_my_mz.mx = 0.5 * ( 
			 ((double) xcr.nunit[0]) * _vec[0][0] +
			 ((double) xcr.nunit[1]) * _vec[1][0] );
    
    mx_my_mz.my = 0.5 * (
			 ((double) xcr.nunit[0]) * _vec[0][1] + 
			 ((double) xcr.nunit[1]) * _vec[1][1] );

    mx_my_mz.mz = sumz / ((double) natoms);

    /*      
    fprintf(stderr,"(++,my,mz) = (%f,%f,%f)\n", mx_my_mz.mx, mx_my_mz.my, mx_my_mz.mz);
    */

  } else {

    mx_my_mz.mx = 0.5 * ((double)xcr.nunit[0] * _vec[0][0] + 
		 (double)xcr.nunit[1] * _vec[1][0] +
		 (double)xcr.nunit[2] * _vec[2][0]);
    
    mx_my_mz.my = 0.5 * ((double)xcr.nunit[0] * _vec[0][1] + 
		 (double)xcr.nunit[1] * _vec[1][1] +
		 (double)xcr.nunit[2] * _vec[2][1]);
    
    mx_my_mz.mz = 0.5 * ((double)xcr.nunit[0] * _vec[0][2] + 
		 (double)xcr.nunit[1] * _vec[1][2] +
		 (double)xcr.nunit[2] * _vec[2][2]);
  }
  /*
  fprintf(stderr,"(mx,my,mz) = (%f,%f,%f)\n", mx_my_mz.mx, mx_my_mz.my, mx_my_mz.mz);
  fprintf(stderr,"NUNITs     = (%d,%d,%d)\n", xcr.nunit[0], xcr.nunit[1], xcr.nunit[2]);
  */

  /* work-around for BUG connected with gcc -O3 compilation */
  if ( isnan(mx_my_mz.mx) )  mx_my_mz.mx = sumx / ((double) natoms);
  if ( isnan(mx_my_mz.my) )  mx_my_mz.my = sumy / ((double) natoms);
  if ( isnan(mx_my_mz.mz) )  mx_my_mz.mz = sumz / ((double) natoms);
}
