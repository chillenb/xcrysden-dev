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
 * Source: $XCRYSDEN_TOPDIR/C/xcLabels.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <stdio.h>
#include <string.h>
#include "struct.h"
#include "xcfunc.h"

/* --- function prototypes --- */
void makeRasterFont(void);
void makeAtomLabels(void);
void makeXYZLabels(void);
static void arrow(void);
void makeCrdList(void);
void makeTemp3D2DList(void);

extern GLubyte rasters[][13];
extern char *element[];

/* extern GLuint xcGenLists( GLsizei i ); */

GLuint fontOffset;
RasterFontSize rf = {8, 13, 0.0, 0.0, 0.0, 0.0};


void 
makeRasterFont(void)
{
    GLuint i;
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); /* 1 for byte-alignment */

    fontOffset = xcGenLists (128);
    for (i = 32; i < 127; i++) {
	glNewList(i+fontOffset, GL_COMPILE);
	  glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, rasters[i-32]);
	glEndList();
    }
}


/*****************************************************************************/
/* make lists for atom labels (all elements)
 */
void
makeAtomLabels(void)
{ 
  int i;

  atomlabelOffset = xcGenLists(MAXNAT + 1);
  for(i=0; i < MAXNAT + 1; i++) {
    glNewList(atomlabelOffset + i, GL_COMPILE);
      glPushAttrib (GL_LIST_BIT);
      glListBase(fontOffset);
      glCallLists(strlen(element[i]), GL_UNSIGNED_BYTE, (GLubyte *) element[i]);
      glPopAttrib ();
    glEndList();
  }
}


/* also make lists for x,y,z label that go with coordinate sistem, that will be
 * displayed
 */
void
makeXYZLabels(void)
{ 
  int i;
  char *xyz[3] = { "x", "y", "z" };

  xyzlabelOffset = xcGenLists(3);
  for(i=0; i < 3; i++) {    
    glNewList(xyzlabelOffset + i, GL_COMPILE);
      glPushAttrib (GL_LIST_BIT);
      glListBase(fontOffset);
      glCallLists(1, GL_UNSIGNED_BYTE, (GLubyte *) xyz[i]);
      glPopAttrib ();
    glEndList();
  }
}


#define T 0.8
#define S 0.0625 

GLuint arrowList = -1;
GLuint crdList = -1;
GLuint xyplaneList = -1;

static void
arrow(void)
{
  int i;
  double arrow[18][3] = {
    /* arrow for X axis */
    { T, S, 0},
    { T,-S, 0},
    {1.0, 0.0, 0.0},
    { T, 0,-S},
    { T, 0, S},
    {1.0, 0.0, 0.0},
    /* arrow for Y axis */
    {-S, T, 0},
    { S, T, 0},
    {0.0, 1.0, 0.0},
    { 0, T, S},
    { 0, T, -S},
    {0.0, 1.0, 0.0},
    /* arrow for Z axis */
    { S, 0, T},
    {-S, 0, T},
    {0.0, 0.0, 1.0},
    { 0,-S, T},
    { 0, S, T},
    {0.0, 0.0, 1.0}};

  arrowList = xcGenLists(1);

  glNewList( arrowList, GL_COMPILE );       
  /* X axis */
    glBegin(GL_LINES);
      glVertex3d( 0.0, 0.0, 0.0);
      glVertex3d( 1.0, 0.0, 0.0);
    glEnd();
    glBegin(GL_TRIANGLES);
      for(i=0; i<6; ) {    
        glVertex3dv( arrow[i++] );
	glVertex3dv( arrow[i++] );
	glVertex3dv( arrow[i++] );
      }
    glEnd();
    /* Y axis */
    glBegin(GL_LINES);
      glVertex3d( 0.0, 0.0, 0.0);
      glVertex3d( 0.0, 1.0, 0.0);
    glEnd();
    glBegin(GL_TRIANGLES);
      for(; i<12; ) {    
        glVertex3dv( arrow[i++] );
	glVertex3dv( arrow[i++] );
	glVertex3dv( arrow[i++] );
      }
    glEnd();
    /* Z axis */
    glBegin(GL_LINES);
      glVertex3d( 0.0, 0.0, 0.0);
      glVertex3d( 0.0, 0.0, 1.0);
    glEnd();
    glBegin(GL_TRIANGLES);
      for(; i<18; ) {    
        glVertex3dv( arrow[i++] );
	glVertex3dv( arrow[i++] );
	glVertex3dv( arrow[i++] );
      }    
    glEnd();
  glEndList();
}


void
makeCrdList(void)
{
  /* load arrowList */
  arrow();
    
  crdList = xcGenLists(1);
  printf("makeCrdList> crdList=%d\n", crdList);
  fflush(stdout);
  /* first make XYZ axes */
  glNewList( crdList, GL_COMPILE );
    glCallList(arrowList);
    /* now X label */
    glRasterPos3f( 1.1,-0.1, 0.0);
    glCallList(xyzlabelOffset);  
    /* now Y label */
    glRasterPos3f(-0.1, 1.1, 0.0);
    glCallList(xyzlabelOffset + 1);
    /* now Z label */
    glRasterPos3f(-0.1,-0.1, 1.1);
    glCallList(xyzlabelOffset + 2);
  glEndList();

  xyplaneList = xcGenLists(1);
  glNewList( xyplaneList, GL_COMPILE );
    glBegin(GL_POLYGON);
      glVertex3f(-.2, .2, 0.2);
      glVertex3f( .2, .2, 0.2);
      glVertex3f( .2,-.2, 0.2);
      glVertex3f(-.2,-.2, 0.2);
    glEnd();
  glEndList();
}
    

GLuint tempDisable3Dlist, tempEnable3Dlist;
/* when we are displaying coor-sist, we are switching between modes; it's 
 * more efficient to "pack" this switching in displayLists
 */
void
makeTemp3D2DList(void)
{
  /* when we are displaying coord-sist in 3D mode, 
   * we must Disable some things
   */
  tempDisable3Dlist = xcGenLists(1);
  glNewList( tempDisable3Dlist, GL_COMPILE );
    glDisable( GL_LIGHTING );
    glDisable( GL_LIGHT0 );
    glDisable( GL_DITHER );
  glEndList();

  /* than we must go back */
  tempEnable3Dlist = xcGenLists(1);
  glNewList( tempEnable3Dlist, GL_COMPILE );
    glEnable( GL_LIGHTING );
    glEnable( GL_LIGHT0 );
    glEnable( GL_DITHER );
  glEndList();
}
  

    
