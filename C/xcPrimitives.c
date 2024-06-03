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
 * Source: $XCRYSDEN_TOPDIR/C/xcPrimitives.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <GL/gl.h>
#include "vector.h"
#include "primitives.h"

/* xcPrimitives.c */
void xcSolidCage(CellCage cage);
void xcWireCage(CellCage cage);
void xcSolidVector(RenderVectors vec);
void xcParallelogram(int type, float p[4][3], float nml[3]);

extern void xcSolidCylinder(GLdouble radius, GLdouble height);
extern void xcSolidVector( RenderVectors vec );
extern void xcSolidCone(GLdouble baseradius, GLdouble topradius, GLdouble height);

void
xcSolidCage( CellCage cage )
{
  register int i, j;

  for (i=0; i<6; i++) {
    glBegin( GL_QUADS );
      for (j=0; j<4; j++) {
	glNormal3fv( cage.normal[i][j] );
	glVertex3fv( cage.vertex[i][j] );
      }
    glEnd();
  }
}
 
  
void
xcWireCage( CellCage cage )
{
  register int i, j;

  glDisable( GL_LIGHTING );
  for (i=0; i<6; i++) {
    glBegin( GL_LINE_LOOP );
      for (j=0; j<4; j++) {
	glVertex3fv( cage.vertex[i][j] );
      }
    glEnd();
  }
  glEnable( GL_LIGHTING );
}
 
  
void
xcSolidVector( RenderVectors vec )
{
  glPushMatrix();
    glRotated( vec.arrfi, vec.arrx, vec.arry, vec.arrz );
    xcSolidCylinder( vec.vecthick, vec.vecl );
  glPopMatrix();
  /* now arrow of vector */
  glPushMatrix();
    glTranslated( vec.coor[0][1], vec.coor[1][1], vec.coor[2][1] );
    glRotated( vec.arrfi, vec.arrx, vec.arry, vec.arrz );
    xcSolidCone( vec.arrthick, 0.0, vec.arrl );
  glPopMatrix();
}


void
xcParallelogram( int type, float p[4][3], float nml[3] )
{
  register int i;

  if ( type == XCPRIM_WIRE ) {
    glBegin( GL_LINE_LOOP );
      glNormal3fv( nml );
      for(i=0; i<4; i++) {
	glVertex3fv( p[i] );
      }
    glEnd();      
  }
  else if ( type == XCPRIM_SOLID || type == XCPRIM_SOLID_AND_BORDER ) {
    glBegin( GL_QUADS );
      glNormal3fv( nml );
      for(i=0; i<4; i++) {
	glVertex3fv( p[i] );
      }
    glEnd();
  }
  if ( type == XCPRIM_SOLID_AND_BORDER ) {
    float color[4] = { 1.0, 1.0, 1.0, 1.0 };
    /* disable transparency */
    glDisable(GL_BLEND); 
    glDisable( GL_LIGHTING );
    glDepthMask(GL_TRUE);  
    /* here border is defined; it is white and the width = 2.0 */
    glLineWidth( 2.0 );
    glColor4fv( color );
    glBegin( GL_LINE_LOOP );
      /* front side */
      glNormal3fv( nml );
      for(i=0; i<4; i++) {
	glVertex3fv( p[i] );
      }
    glEnd();
    glBegin( GL_LINE_LOOP );
      /* back side */
      glNormal3f( -nml[0], -nml[1], -nml[2] );
      for(i=3; i>=0; i--) {
	glVertex3fv( p[i] );
      }
    glEnd();      
  }
  glEnable( GL_LIGHTING );
}
