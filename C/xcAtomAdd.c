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
* Source: $XCRYSDEN_TOPDIR/C/xcAtomAdd.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <togl.h> 
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tk.h>
#include <GL/glu.h>
#include "struct.h"
#include "displayfunc.h"
#include "vector.h"
#include "xcfunc.h"
#include "memory.h"

#define ADDLINE_PAT 0xAAAA  /* addline stipple pattern */
#define ADDLINE_WIDTH 2.0   /* width of "addline" */


/*GLuint BasicVectorsList, AtomAddList;*/
GLdouble points[8][3];
GLdouble addatompos[3];
extern GLuint tempDisable3Dlist;
extern GLuint tempEnable3Dlist;

static RenderVectors renvec[3];
static double mxx, myy, mzz;

/* --- FUNCTIONS PROTOTYPES --- */
int XC_AtomAddCmd(ClientData clientData, Tcl_Interp *interp,
		  int argc, const char *argv[]);
static void SetBasicVectors( double vec[][4], int n );
static void xcAddAtomUpdate(double af, double bf, double cf, char *addatomdata);
void xcDisplayAddAtomBox(void);
static void RenderBasicVectors(void) ;

/*static void xcMakeAtomAddList(void);*/


/* --- xcDisplayFunc.c --- */
extern void xcSolidCylinder (GLdouble radius, GLdouble height);
extern void xcSolidCone (GLdouble baseradius, GLdouble topradius, 
			 GLdouble height);
/* extern void (*xcDisplay)(void); */

/* --- xcSelect.c --- */
extern void GetSelCylinderPar(double x21, double y21, double z21, 
			      double *xrvb, double *yrvb, double *zrvb, 
			      double *fibond, double *bondl);



/* XC_AtomAddCmd --> inplementation of 'xc_atomadd' custom Tcl command,
   which is used for Cell-Adding type of ATOMINSE
   crystal command
   * ----------------------------------------------------------------------------
   * Usage: xc_atomadd <toglName> begin|update|clean ?AF? ?BF? ?CF? 
   *		
   *              xc_atomadd begin -- before we begin atomadd procedure 
   *                                  we must prepare and do some initialisations
   *
   *              xc_atomadd update-- updating
   *
   *              xc_atomadd clean -- clean-up
   */	
int 
XC_AtomAddCmd(ClientData clientData, Tcl_Interp *interp,
	      int argc, const char *argv[])
{
  struct Togl *togl;
  char *csqn = (char*) Tcl_Alloc( sizeof(char) * 256 );
  GLdouble af,bf,cf; /* factors for vectors AF BF CF */

  if (argc < 3) {
    Tcl_SetResult(interp, "Usage: xc_atomadd <toglName> begin|update|clean ?AF? ?BF? ?CF?", TCL_STATIC);
    return TCL_ERROR;
  }

  /* find togl associated with toglName */
  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* if structure is not opened, just return silently */
  if ( !VPf.stropened ) return TCL_OK;

  /* so far just XC_3D is implemented; if not XC_3D return silently */
  if ( dimType != XC_3D ) return TCL_OK;

  /***************************************************************************/
  /* first "XC_ATOMADD <toglName> BEGIN" */
  if ( strcmp(argv[2],"begin") == 0 ) {
    int n;
    /* SO FAR ONLY 3D XC_ATOMADD IMPLEMENTATION */
    /* flag to know when we are in "atomadd mode" */
    VPf.atomadd = GL_TRUE;

    /* 3D ONLY: make basic vectors appear on the screen */    
    /* if xcr.lprimvec == GL_TRUE --> deal with primitiv vectors */
    /* if xcr.lconvvec == GL_TRUE --> deal with primitiv vectors */
    /* if xcr.dim  == 3 --> crystal
     *             == 2 --> SLAB; put C vector to (0.0,0.0,1.0)
     *		   == 1 --> polymer; put B,C vectors to 
     *                               (0.0,(1.0|0.0),(0.0|1.0))
     */
    n = xcr.dim;
    if ( xcr.celltype == XCR_PRIMCELL ) SetBasicVectors( vec.prim, n ); 
    if ( xcr.celltype == XCR_CONVCELL ) SetBasicVectors( vec.conv, n ); 

    /* default values for A B C fractions will be 0.3 0.3 0.3, that is what we
       will return */
    af = 0.3;
    bf = 0.3;
    cf = 0.3;
    /* now update and return current position */
    xcAddAtomUpdate(af, cf, bf, csqn);
    Tcl_SetResult(interp,csqn,TCL_DYNAMIC);
    /* now update display */
    xcDisplayAddAtomBox();
    glFlush();
    Togl_SwapBuffers(togl);
  }

  /***************************************************************************/
  /* now implement a "XC_ATOMADD <toglName> UPDATE %AF %BF %CF"; 
     %AF,%BF,%CF are fractions of AF BF CF crystal vectors */
  /* at the end we will calculate current position */
  else if ( strcmp(argv[2],"update") == 0 ) {

    if ( !VPf.atomadd ) {
      Tcl_SetResult(interp, "calling \"xc_atomadd <toglName> update\" command, before initialising atom-add process with \"xc_atomadd <toglName> begin\"", TCL_STATIC);
      return TCL_ERROR;
    }    

    /* argc must be 5 */
    if ( argc != 6 ) {
      Tcl_SetResult(interp, "wrong usage of \"xc_atomadd <toglName> update\" command, must be \"xc_addatom <toglName> update <AF> <BF> <CF>\"", TCL_STATIC);
      return TCL_ERROR;
    }
    /* argv[3],argv[4],argv[5] must be integers */
    if ( Tcl_GetDouble(interp,argv[3],&af) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"got invalid value for <AF>, while executing %s %s %s %s %s %s\n",argv[0],argv[1],argv[2],argv[3],argv[4],argv[5]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( Tcl_GetDouble(interp,argv[4],&bf) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"got invalid value for <BF>, while executing %s %s %s %s %s %s\n",argv[0],argv[1],argv[2],argv[3],argv[4],argv[5]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( Tcl_GetDouble(interp,argv[5],&cf) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"got invalid value for <CF>, while executing %s %s %s %s %s %s\n",argv[0],argv[1],argv[2],argv[3],argv[4],argv[5]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    /* now update and return current position */
    xcAddAtomUpdate(af, bf, cf, csqn);
    Tcl_SetResult(interp,csqn,TCL_DYNAMIC);
    /* update the display */
    Togl_PostRedisplay(togl);
  }

  /***************************************************************************/
  /* XC_ATOMADD <toglName> CLEAN */
  else if ( strcmp(argv[2],"clean") == 0 ) {
    if ( !VPf.atomadd ) {
      Tcl_SetResult(interp, "calling \"xc_atomadd <toglName> clean\" command, before initialising atom-add process with \"xc_atomadd <toglName> begin\"", TCL_STATIC);
      return TCL_ERROR;      
    }    
    VPf.atomadd = GL_FALSE;
    /* now update display */
    Togl_PostRedisplay(togl);
  }
  /***************************************************************************/
  /* unknown option */
  else {
    char rss[1024];
    snprintf(rss, sizeof(rss),"unknown option \"%s\" for \"xc_atomadd\" command, must be one of begin, update, clean",argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  return TCL_OK;
}


/*****************************************************************************/
static void 
SetBasicVectors( double vec[][4], int n )
{
  int i,j;

  /* for polymer & SLAB modify (B) and C vector to 
     (0.0,(1.0|0.0),(0.0|1.0)) */
  for (i = 0; i<=2; i++) {
    if ( i >= n ) {
      vec[i][i] = 1.0; /* modified is just part that should be 1.0 */
    }
    /* arrow will be 0.25 of vector length */
    for (j = 0; j<=2; j++) {
      renvec[i].coor[j][0] = 0.0;
      renvec[i].coor[j][1] = (1.0 - VECTOR_ARROWSIZE) * vec[i][j];
      renvec[i].coor[j][2] = vec[i][j];
    }

    renvec[i].vecthick = VECTOR_THICKF * rrod;
    renvec[i].arrthick = VECTOR_ARRTHICKF * renvec[i].vecthick;

    fprintf(stderr,"vec-thick: %f %f\n", renvec[i].vecthick, renvec[i].arrthick);

    /* convert cartesian coordinates to gluCylinder coor. */
    GetSelCylinderPar(renvec[i].coor[0][1] - renvec[i].coor[0][0],
		      renvec[i].coor[1][1] - renvec[i].coor[1][0],
		      renvec[i].coor[2][1] - renvec[i].coor[2][0],
		      &renvec[i].vecx, &renvec[i].vecy, &renvec[i].vecz,
		      &renvec[i].vecfi, &renvec[i].vecl);

    GetSelCylinderPar(renvec[i].coor[0][2] - renvec[i].coor[0][1], 
		      renvec[i].coor[1][2] - renvec[i].coor[1][1],
		      renvec[i].coor[2][2] - renvec[i].coor[2][1],
		      &renvec[i].arrx, &renvec[i].arry, &renvec[i].arrz,
		      &renvec[i].arrfi, &renvec[i].arrl);
  }

  mxx = mx;
  myy = my;
  mzz = mz;
  /**************************************************************************/
  /* maybe for MOLECULE it's better to draw vectors from MASS-CENTER,       */
  /* because of the presence of point-groups                                */
  /**************************************************************************/
  if ( xcr.dim == 0 ) {
    mxx = 0.0;
    myy = 0.0;
    mzz = 0.0;
  }
}


static void
RenderBasicVectors(void) {
  int i;

  glPushMatrix();
  glTranslated( -mxx, -myy, -mzz );
  LoadCageOrVecMaterial( GLPAR_PRIMVEC );
  for (i=0; i<3; i++) {
    xcSolidVector( renvec[i] );
  }
  LoadStructMaterial();    
  glPopMatrix();
}

/*****************************************************************************/
/* CALCULATE NEW CURRENT ADDATOM POSITION & RETURN IT */
static void
xcAddAtomUpdate(double af, double bf, double cf, char *addatomdata)
{  
  int i;

  /* points for cube-box edges */
  for(i=0; i<3; i++) {
    points[0][i] = 0.0;
    points[1][i] = af*renvec[0].coor[i][2];
    points[2][i] = af*renvec[0].coor[i][2] + bf*renvec[1].coor[i][2];
    points[3][i] = bf*renvec[1].coor[i][2];
    points[4][i] = cf*renvec[2].coor[i][2];
    points[5][i] = af*renvec[0].coor[i][2] + cf*renvec[2].coor[i][2];
    points[6][i] = points[5][i] + bf*renvec[1].coor[i][2];
    points[7][i] = bf*renvec[1].coor[i][2] + cf*renvec[2].coor[i][2];

    addatompos[i] = af*renvec[0].coor[i][2] + bf*renvec[1].coor[i][2] + cf*renvec[2].coor[i][2];
  }    

  sprintf(addatomdata,"%.10f   %.10f   %.10f", addatompos[0], addatompos[1], addatompos[2]);
}


void
xcDisplayAddAtomBox(void)
{
  double mxx, myy, mzz;

  mxx = mx;
  myy = my;
  mzz = mz;
  /**************************************************************************/
  /* maybe for MOLECULE it's better to draw vectors from MASS-CENTER,       */
  /* because of the presence of point-groups                                */
  /**************************************************************************/
  if ( xcr.dim == 0 ) {
    mxx = 0.0;
    myy = 0.0;
    mzz = 0.0;
  }

  glLineStipple(2, ADDLINE_PAT);
  glEnable(GL_LINE_STIPPLE);
  glCallList(tempDisable3Dlist);
  glLineWidth( ADDLINE_WIDTH );
  /* make a addatombox */
  /* box is made from this poits:
   *                              01,12,23,30
   *                              04,15,26,37
   *                              45,56,67,74
   */
  glTranslated( -mxx, -myy, -mzz ); /* (mx,my,mz) is origin shift vector */
  glBegin(GL_LINE_STRIP);
  glVertex3dv( points[0]);
  glVertex3dv( points[1]);
  glVertex3dv( points[2]);
  glVertex3dv( points[3]);
  glVertex3dv( points[0]);
  glEnd();
  glBegin(GL_LINES);
  glVertex3dv( points[0]);
  glVertex3dv( points[4]);
  glVertex3dv( points[1]);
  glVertex3dv( points[5]);
  glVertex3dv( points[2]);
  glVertex3dv( points[6]);
  glVertex3dv( points[3]);
  glVertex3dv( points[7]);
  glEnd();
  glBegin(GL_LINE_STRIP);
  glVertex3dv( points[4]);
  glVertex3dv( points[5]);
  glVertex3dv( points[6]);
  glVertex3dv( points[7]);
  glVertex3dv( points[4]);
  glEnd();
  /* make a point on atomadd-position */
  glPointSize(6);
  glBegin(GL_POINTS);
  glVertex3dv( addatompos );
  glEnd();
  glTranslated( mxx, myy, mzz );

  glCallList(tempEnable3Dlist);
  glDisable(GL_LINE_STIPPLE);

  RenderBasicVectors();
} 


/*
  static void 
  xcMakeAtomAddList(void) 
  {
  float color[3] = {1.0, 1.0, 1.0};
  GLUquadricObj *quadObj;
  quadObj = gluNewQuadric ();
  glNewList(AtomAddList, GL_COMPILE);
  glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, 
  color);
  gluQuadricDrawStyle (quadObj, GLU_FILL);
  gluQuadricNormals (quadObj, GLU_FLAT);
  gluSphere (quadObj, 0.1, 4, 4);
  glEndList();
  }
*/
