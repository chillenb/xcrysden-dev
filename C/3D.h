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
 * Source: $XCRYSDEN_TOPDIR/C/3D.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

/* definition that are used for 3D structures/displayModes */

GLdouble NewMat[4][4];
GLdouble OldMat[4][4];
GLdouble MajorMat[4][4];
GLdouble xcGLvec[16];	      		      

    
/*****************************************************************************/
typedef struct {
  int flag;    /* flag for BOND or FRAME */
  int bondFlag; 
  int nat;     /* corresponding atomic number 
		* from NAT all other atribbutes are assigned */
  int sqn;     /* sequential number of atom; from sqn other attribbutes are
		  assigned */
  GLuint list; /* identifier for BallList to use when making/rendering balls */
  double x1;   
  double y1;   /* for ATOMS only X1 Y1 Z1 are needed */
  double z1;
  double xrvb;  /* XRotationVectorBond; */
  double yrvb;  /* (xrvb, yrvb, zrvb) is vector around which I rotatte */
  double zrvb;  /* auxCylinder for fibond; lenght of bond is bondl */
  double fibond;
  double bondl;
} BondFrame3D;
BondFrame3D *coor3D;

