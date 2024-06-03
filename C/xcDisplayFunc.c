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
 * Source: $XCRYSDEN_TOPDIR/C/xcDisplayFunc.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj  
 * ------                                                                    *
 * Modified by Eric Verfaillie ericverfaillie@yahoo.fr EV                    *
 * may 2004                                                                  *
 * modifcations are near EV comments                                         *
 *****************************************************************************

*/

#include <togl.h> 
#include <stdio.h> 
#include <math.h> 
#include <stdlib.h>
#include <string.h> 
#include <GL/glu.h> 
#include "struct.h" 
#include "3D.h" 
#include "xcGLparam.h" 
#include "isosurf.h"
#include "molsurf.h"
#include "xcfunc.h"

#define APPR_SCREENSIZE 1024
#define SELCOL 0
#define SELLINE_WIDTH 2.0   /* width of "selline" */
#define SELLINE_PAT 0xAAAA  /* selline stipple pattern */

#define SPHERE_TESS_PREFACTOR 2
#define BOND_TESS_PREFACTOR   1


/*********************************************/
/*                 EV                        */
/*********************************************/
#include "anaglyph.h"   

GLdouble radiu;
CAMERA camera;
XYZ origin = {0.0,0.0,0.0};
/*********************************************/
/*                 EV                        */
/*********************************************/



extern struct Togl *mesa_togl; 

extern Options3D      is; 
extern RasterFontSize rf; 
extern PERSPECTIVE    persp;
extern AtomicLabel    *atomLabel, globalAtomLabel;
extern short          *do_not_display_atomlabel;
extern char           *element[];
extern GLuint         fontOffset;
extern GLfloat        AtCol_Ambient_by_Diffuse;
extern XCantialias    antialias;
extern ForceVector    FV;
extern XYZ_Attrib     xyz;
OrthoProj             ort; 

extern realTimeMove makeMovie;

/***************************************************************************** 
 * HERE LIST VARIABLES ARE DIFINED                                           */
GLuint PointList = -1; 
/*
  static GLuint SolidSpaceFillList = -1; 
  static GLuint WireSpaceFillList = -1; 
  static GLuint SolidBallList = -1; 
  static GLuint WireBallList = -1; 
  static GLuint SolidStickList = -1; 
  static GLuint WireStickList = -1; 
*/
/*
  static GLuint SolidFrameList = -1; 
  static GLuint WireFrameList = -1; 
  static GLuint LineFrameList = -1; 
*/
/*  static GLuint LabelBallList = -1;  */
/*  static GLuint LabelSpaceList = -1;  */

extern GLuint crdList; /* crdList is list for coor-sist.  */
extern GLuint xyplaneList; /* list for XY plane that goes with coor-sist  */

extern GLuint tempDisable3Dlist; 
extern GLuint tempEnable3Dlist; 
/* --- this is for selection-lists --- */ 
/*extern GLuint SphereSelList[MAXSEL];*/
/*extern GLuint LineSelList[MAXSEL-1];*/
/* --- this is for cell-adding-type-procedure of ATOMINSE Crystal command  */
/*extern GLuint BasicVectorsList, AtomAddList; */
/****************************************************************************/ 



/*****************************************************************************/

/* function prototypes  */
GLuint xcGenLists( GLsizei i );
void (*xcDisplay)(struct Togl *togl) = 0x0; 
void xcDisplayFunc( void (*Func)(struct Togl *togl) ); 
int xcToglDisplayFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv);
void xcDummyDisplay( struct Togl *togl ); 
void xcGenDispList(void); 
void xcWireFrame2D(struct Togl *togl); 
void xcMakePointList(void); 
void xcPointLine2D(struct Togl *togl); 
void xcBall2D(int iatom);
/*void xcMakeBallLists(void); 
  void xcMakeNewBallList( GLdouble *sizeArray , int natom ); */
void xcBallStick12D(struct Togl *togl); 
void xcPipeBall2D(struct Togl *togl); 
void xcBallStick22D(struct Togl *togl); 
void xcSpaceFill2D(struct Togl *togl);
void EnableOr2D_Or3D(void); 
static void _ambient_by_diffuse(GLfloat dst[4], GLfloat scale, GLfloat src[4]);
void xcRenderSolidBalls3D(void);
void xcRenderWireBalls3D(void);
void xcRenderBonds(GLenum type); 
void xcRenderSolidSpaceFills3D(void);
void xcRenderWireSpaceFills3D(void);
void xcRenderSolidFrame3D(void);
void xcRenderWireFrame3D(void);
void xcRenderLineFrame3D(void);
static void xcLineStippleFrame(int ii); 
static void MatInv33_44(double m[4][4], double minv[4][4], double *det); 
static void MultMat44_33Vec3(double mat[4][4], double vec[3], double res[3]); 
static void _makeAtomLabels3D(double *rad);
void xcMakeBallLabel3D(void); 
void xcMakeSpaceLabel3D(void); 
void xcDisplay3D(struct Togl *togl);
static void _xcDisplay3D(struct Togl *togl);
static void __xcDisplay3D(struct Togl *togl);
void xcMaybeDestroyLists(void); 
void xcMakeProjection2D(const char *mode); 
void xcMakeProjection3D(const char *mode); 
static void _xcMakeProjection(GLdouble size);
void xcClearScreen(struct Togl *togl); 
void xcChangeBackground(struct Togl *togl, GLclampf bckg[4]); 
void xcDisplayXYZ(void); 
static void goToPrevProj(void); 
void xcTurnXYZOff(struct Togl *togl); 

void xcWireSphere (GLdouble radius); 
void xcSolidSphere (GLdouble radius); 
void xcWireCylinder (GLdouble radius, GLdouble height); 
void xcSolidCylinder (GLdouble radius, GLdouble height); 
void xcSolidBond (GLdouble radius, GLdouble height, int bondFlag); 
void xcSolidCone (GLdouble baseradius, GLdouble topradius, GLdouble height); 
static int CalcTessFactor(void);
static void HandleDisplay(struct Togl *togl);
static void draw_scene(void);
static void CameraHome(void);
static void Normalise(XYZ *p);
/*static XYZ CalcNormal(XYZ p, XYZ p1, XYZ p2);*/

/* --- extern function prototypes ---  */
extern int MakeSticks1(int i,
		       GLdouble *x1, GLdouble *y1, GLdouble *z1, 
		       GLdouble *x2, GLdouble *y2, GLdouble *z2, 
		       GLdouble *x3, GLdouble *y3, GLdouble *z3, 
		       GLdouble *x4, GLdouble *y4, GLdouble *z4);
extern int MakeSticks2(int i, int col, int flag); 

/* --- readstrf.c ---*/ 
extern void FindMaxRad(void); 

/* --- xcSelect.c --- */
extern void xcRenderSelAtoms3D(void);
extern void xcRenderSelBonds3D(void);

/* --- xcAtomAdd.c --- */ 
extern void xcDisplayAddAtomBox(void); 

/* --- isorender.c --- */ 
extern void xcRenderIsosurf(int obj); 
extern void xcRenderColorplane(int obj); 
extern void xcRenderSurface(void); 

/* --- xcSuperCell.c --- */ 
extern void (*xcSuperCell)(void); 

/* --- xcIsoSpaceSel.c --- */ 
extern void IsoSpaceSel_Parallelogram(void); 
extern void IsoSpaceSel_3D(void); 

/* --- xcviewport.c --- */ 
extern void MaybeClipAndMakeProjection(void); 


/* --- cryNewContext.c ---  */
extern NEW_WIN_CONTEXT *FindWinContextByTogl(struct Togl *togl); 

/* xcdebug */
extern void xcdebug(const char *text);

/* xcFont.c */
extern void xcFont_PrintString (const char *s);

/* sgiAux.c */
extern GLuint findList1 (int lindex, GLdouble *paramArray, int size); 
extern GLuint makeModelPtr1 (int lindex, GLdouble *sizeArray, int count); 

GLuint 
xcGenLists( GLsizei i ) 
{ 
  return glGenLists (i);
  /*
    static GLuint nlist = 1; 
    GLuint n; 
    
    n = nlist; 
    nlist += i; 
    return n;
  */ 
} 



/*****************************************************************************/
void
xcDisplayFunc( void (*Func)(struct Togl *togl) ) 
{ 
  xcDisplay = Func; 
} 
/*****************************************************************************/


/*****************************************************************************/
/*void cryDisplayFunc( void (*DispFunc)(struct Togl *togl), struct Togl *togl );*/ 
int
xcToglDisplayFunc(ClientData clientData, Tcl_Interp *interp,
                  int objc, Tcl_Obj *const *objv)
{
  Togl *togl;
  /*static int toglIsInteractive = 0;*/

  if (objc != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, "pathName");
    return TCL_ERROR;
  }
  if (Togl_GetToglFromObj(interp, objv[1], &togl) != TCL_OK) {
    return TCL_ERROR;
  }
  
  if ( togl == mesa_togl ) {
    if (xcDisplay) (*xcDisplay)(togl); 
  } else { 
    NEW_WIN_CONTEXT *wc = FindWinContextByTogl( togl ); 
    /*cryDisplayFunc( wc->xcDisplay, togl );*/ 
    if (wc->xcDisplay) wc->xcDisplay( togl ); 
  }
  return TCL_OK;
} 
/* 
void cryDisplayFunc( void (*DispFunc)(struct Togl *togl), struct Togl *togl ) 
{ (*DispFunc)(togl); } */ 
/*****************************************************************************/


/*****************************************************************************/
void 
xcDummyDisplay( struct Togl *togl ) 
{   
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ); 
  glClearColor( bg[0], bg[1], bg[2], bg[3] ); 
  Togl_SwapBuffers(togl); 
  xcdebug("dummy display function"); 
} 



/* this routine generate display Lists */ 
void 
xcGenDispList(void) 
{ 
  /*int i; */

  /* generate all Display Lists; actually we need just to reserve numbers 
   * for lists */ 
  PointList = xcGenLists(1);   
  /*
    SolidSpaceFillList = xcGenLists(1); 
    WireSpaceFillList = xcGenLists(1); 
    SolidBallList = xcGenLists(1); 
    WireBallList = xcGenLists(1); 
    SolidStickList = xcGenLists(1); 
    WireStickList = xcGenLists(1); 
  */
  /*
    SolidFrameList = xcGenLists(1);
    WireFrameList = xcGenLists(1); 
    LineFrameList = xcGenLists(1); 
  */
/*    LabelBallList = xcGenLists(1); */ 
/*    LabelSpaceList = xcGenLists(1); */ 
  /*for(i=0; i<MAXSEL; i++) SphereSelList[i] = xcGenLists(1);*/
  /*for(i=0; i<MAXSEL-1; i++) LineSelList[i] = xcGenLists(1);*/
  /*BasicVectorsList = xcGenLists(1); */
  /*AtomAddList = xcGenLists(1); */
} 




/* xcWireFrame2D is xcDisplayFunc for WireFrames   
 * -------------------------------------------------------------- 
 *                xcWireFrame2D display a WireFrames 
 *               Here I'm dealing with 'coor' structure 
 */ 

void 
xcWireFrame2D(struct Togl *togl) 
{ 
  int i; 

  glClear(GL_COLOR_BUFFER_BIT); 
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glLineWidth( (GLfloat) VPf.WFlinewidth ); 
  for(i = 1; i <= tmp_nobjects; i++) { 
    glBegin(GL_LINES); 
       if ( VPf.unibond) glColor3fv(unibondCol);
       else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] ); 
       glVertex3d( (coor + *(iwksp + i))->x1, 
		   (coor + *(iwksp + i))->y1,
		   (coor + *(iwksp + i))->z1); 
       glVertex3d( (coor + *(iwksp + i))->x2, 
		   (coor + *(iwksp + i))->y2,
		   (coor + *(iwksp + i))->z2); 
    glEnd(); 
  } 

  /* check if coordinate sistem should be displayed */ 
  if ( VPf.xyzOn ) xcDisplayXYZ(); 

  glFlush(); 
  Togl_SwapBuffers(togl); 

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
} 


/* make list for point */ 
void 
xcMakePointList(void) 
{  
  GLint i, nstep; 
  GLdouble sine, cosine, dstep, r; 
  /* double cos(), sin();  */

  if ( glIsList(PointList) ) glDeleteLists( PointList, 1 ); 

  r = VPf.PLradius / VPf.VPfactor; 
  nstep = 12; 
  dstep = (GLdouble) nstep; 

  /* PointList list was generarted in xcGenDispList */ 
  glNewList( PointList, GL_COMPILE ); 
    glBegin(GL_POLYGON); 
      for(i=0; i < nstep; i++) 
        { 
	  cosine = r * cos((GLdouble) i * 2.0 * PI / dstep); 
	  sine   = r * sin((GLdouble) i * 2.0 * PI / dstep); 
	  glVertex2d( cosine, sine ); 
	} 
    glEnd(); 
  glEndList(); 
} 


/* xcPointLine2D is xcDisplayFunc for PointLines 
 * -------------------------------------------------------------- 
 *                xcPointLine2D display a PointLines 
 *              Here I'm dealing with 'coor' structure 
 */ 
void 
xcPointLine2D(struct Togl *togl) 
{ 
  int i; 

  /*printf("In xcPointLine2D\n",NULL); 
  fflush(stdout);*/ 
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);
 
  glLineWidth( (GLfloat) VPf.PLlinewidth ); 
  /*glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);*/
  for(i = 1; i <= tmp_nobjects; i++) {  
    if ( (coor + *(iwksp + i))->flag == BOND ) 
      { 
	glBegin(GL_LINES); 
	   /* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */ 
           if ( VPf.unibond) glColor3fv(unibondCol);
           else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] ); 
	   glVertex3d( (coor + *(iwksp + i))->x1, 
		       (coor + *(iwksp + i))->y1,
		       (coor + *(iwksp + i))->z1 ); 
	   glVertex3d( (coor + *(iwksp + i))->x2, 
		       (coor + *(iwksp + i))->y2,
		       (coor + *(iwksp + i))->z2 ); 
	   glEnd(); 
      } 
    if ( (coor + *(iwksp + i))->flag == ATOM ) 
      {	 
	/* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */ 
	glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] ); 
	glPushMatrix();
	  glTranslated( (coor + *(iwksp + i))->x1,  
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 ); 
	  glCallList( PointList );	 
	glPopMatrix();	 
      } 
  } 

  /* check if coordinate sistem should be displayed */ 
  if ( VPf.xyzOn ) xcDisplayXYZ(); 

  glFlush(); 
  Togl_SwapBuffers(togl); 

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
} 


void xcBall2D( int iatom ) {
  GLint     nstep, i, atmn;
  GLdouble  sine, cosine;
  GLdouble  sizeArray[3];
  GLuint    displayList[2];

  atmn   = (coor + iatom)->nat; 
  nstep = 3 * CalcTessFactor();
  if ( nstep < 8 )  nstep = 8; 

  sizeArray[0]   = rball[atmn];
  sizeArray[1]   = (GLdouble) nstep; 
  sizeArray[2]   = (GLdouble) VPf.OUTlinewidth;
  displayList[0] = findList1 (BALL, sizeArray, 3);
  displayList[1] = findList1 (OUTLINEBALL, sizeArray, 3);

  if ( displayList[0] == 0 ) {  
    glNewList( makeModelPtr1 (BALL, sizeArray, 3),
	       GL_COMPILE_AND_EXECUTE ); 
      glBegin(GL_POLYGON);
      for(i=0; i < nstep; i++) {
	cosine = rball[atmn] * cos((GLdouble) i * 2.0 * PI / (GLdouble) nstep);
	sine   = rball[atmn] * sin((GLdouble) i * 2.0 * PI / (GLdouble) nstep);
	glVertex2d( cosine, sine );
      }
      glEnd();
    glEndList();

    glNewList( makeModelPtr1 (OUTLINEBALL, sizeArray, 3),
	       GL_COMPILE_AND_EXECUTE );			      
      glColor3f (0.0, 0.0, 0.0); 
      glBegin(GL_LINE_LOOP); 
      for(i=0; i < nstep; i++) { 
	cosine = rball[atmn] * cos((GLdouble) i * 2.0 * PI / (GLdouble) nstep); 
	sine   = rball[atmn] * sin((GLdouble) i * 2.0 * PI / (GLdouble) nstep); 
	glVertex2d( cosine, sine ); 
      } 
      glEnd(); 
    glEndList();
  } else {
    /*glColor3fv (atm.col[(coor + iatom)->sqn]);*/
    glCallList (displayList[0]);
    glCallList (displayList[1]);
  }
}
void xcSmallBall2D( int iatom ) {
  GLint     nstep, i; 
  GLdouble  sine, cosine;
  GLdouble  sizeArray[3];
  GLuint    displayList[2];

  nstep = 3 * CalcTessFactor();
  if ( nstep < 8 )  nstep = 8; 

  sizeArray[0]   = rball[1]; /* size of hydrogen atom */
  sizeArray[1]   = (GLdouble) nstep; 
  sizeArray[2]   = (GLdouble) VPf.OUTlinewidth;
  displayList[0] = findList1 (BALL, sizeArray, 3);
  displayList[1] = findList1 (OUTLINEBALL, sizeArray, 3);

  if ( displayList[0] == 0 ) {  
    /*glColor3fv (atm.col[(coor + iatom)->sqn]);*/
    glNewList( makeModelPtr1 (BALL, sizeArray, 3),
	       GL_COMPILE_AND_EXECUTE ); 
      glBegin(GL_POLYGON);
      for(i=0; i < nstep; i++) {
	cosine = rball[1] * cos((GLdouble) i * 2.0 * PI / (GLdouble) nstep);
	sine   = rball[1] * sin((GLdouble) i * 2.0 * PI / (GLdouble) nstep);
	glVertex2d( cosine, sine );
      }
      glEnd();
    glEndList();

    glNewList( makeModelPtr1 (OUTLINEBALL, sizeArray, 3),
	       GL_COMPILE_AND_EXECUTE );			      
      glColor3f (0.0, 0.0, 0.0); 
      glBegin(GL_LINE_LOOP); 
      for(i=0; i < nstep; i++) { 
	cosine = rball[1] * cos((GLdouble) i * 2.0 * PI / (GLdouble) nstep); 
	sine   = rball[1] * sin((GLdouble) i * 2.0 * PI / (GLdouble) nstep); 
	glVertex2d( cosine, sine ); 
      } 
      glEnd(); 
    glEndList();
  } else {
    /*glColor3fv (atm.col[(coor + iatom)->sqn]);*/
    glCallList (displayList[0]);
    glCallList (displayList[1]);
  }
}
void xcBigBall2D( int iatom ) {
  GLint     nstep, i, atmn; 
  GLdouble  sine, cosine;
  GLdouble  sizeArray[3];
  GLuint    displayList[2];

  atmn   = (coor + iatom)->nat; 
  nstep = 6 * CalcTessFactor();
  if ( nstep < 8 )  nstep = 8; 

  sizeArray[0]   = atrad[atmn]; /* size of hydrogen atom */
  sizeArray[1]   = (GLdouble) nstep; 
  sizeArray[2]   = (GLdouble) VPf.OUTlinewidth;
  displayList[0] = findList1 (BALL, sizeArray, 3);
  displayList[1] = findList1 (OUTLINEBALL, sizeArray, 3);

  if ( displayList[0] == 0 ) {  
    /*glColor3fv (atm.col[(coor + iatom)->sqn]);*/
    glNewList( makeModelPtr1 (BALL, sizeArray, 3),
	       GL_COMPILE_AND_EXECUTE ); 
      glBegin(GL_POLYGON);
      for(i=0; i < nstep; i++) {
	cosine = atrad[atmn] * cos((GLdouble) i * 2.0 * PI / (GLdouble) nstep);
	sine   = atrad[atmn] * sin((GLdouble) i * 2.0 * PI / (GLdouble) nstep);
	glVertex2d( cosine, sine );
      }
      glEnd();
    glEndList();

    glNewList( makeModelPtr1 (OUTLINEBALL, sizeArray, 3),
	       GL_COMPILE_AND_EXECUTE );			      
      glColor3f (0.0, 0.0, 0.0); 
      glBegin(GL_LINE_LOOP); 
      for(i=0; i < nstep; i++) { 
	cosine = atrad[atmn] * cos((GLdouble) i * 2.0 * PI / (GLdouble) nstep); 
	sine   = atrad[atmn] * sin((GLdouble) i * 2.0 * PI / (GLdouble) nstep); 
	glVertex2d( cosine, sine ); 
      } 
      glEnd(); 
    glEndList();
  } else {
    /*glColor3fv (atm.col[(coor + iatom)->sqn]);*/
    glCallList (displayList[0]);
    glCallList (displayList[1]);
  }
}


/* xcBallStick12D is xcDisplayFunc for Ballstick1, 
 * -----------------------------------------------
 *             xcBallStick12D display a BallSticks1
 *            Here I'm dealing with 'coor' structure
 */
void
xcBallStick12D(struct Togl *togl)
{
  register int i;
  GLdouble x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glLineWidth( VPf.OUTlinewidth );
  glEnable( GL_CULL_FACE );
  glCullFace( GL_BACK );
  for(i = 1; i <= tmp_nobjects; i++) {
    if ( (coor + *(iwksp + i))->flag == ATOM ) {
      glPushMatrix();
      glTranslated( (coor + *(iwksp + i))->x1,
		    (coor + *(iwksp + i))->y1,
		    (coor + *(iwksp + i))->z1 );
      glColor3fv (atm.col[(coor + *(iwksp + i))->sqn]);
      xcBall2D (*(iwksp + i));
      glPopMatrix();
    }
    else if ( (coor + *(iwksp + i))->flag == BOND ) {
      MakeSticks1(*(iwksp + i), 
		  &x1, &y1, &z1,   &x2, &y2, &z2, 
		  &x3, &y3, &z3,   &x4, &y4, &z4);
      glBegin(GL_POLYGON);
        if ( VPf.unibond) glColor3fv(unibondCol);
        else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] ); 
	glVertex3d( x1, y1, z1);
	glVertex3d( x2, y2, z2);
	glVertex3d( x3, y3, z3);
	glVertex3d( x4, y4, z4);
      glEnd();

      glBegin(GL_LINE_LOOP);
        glColor3f( 0.0, 0.0, 0.0 );
	glVertex3d( x1, y1, z1);
	glVertex3d( x2, y2, z2);
	glVertex3d( x3, y3, z3);
	glVertex3d( x4, y4, z4);
      glEnd();	      
    }
  }
  glDisable( GL_CULL_FACE );

  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void xcPipeBall2D(struct Togl *togl)
{
  register int i;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glLineWidth( VPf.OUTlinewidth );
  glEnable( GL_CULL_FACE );
  glCullFace( GL_BACK );
  for(i = 1; i <= tmp_nobjects; i++) {
    if ( (coor + *(iwksp + i))->flag == ATOM ) {
      glLineWidth( VPf.OUTlinewidth );
      glPushMatrix();
      glTranslated( (coor + *(iwksp + i))->x1,
		    (coor + *(iwksp + i))->y1, 
		    (coor + *(iwksp + i))->z1 );
      glColor3fv (atm.col[(coor + *(iwksp + i))->sqn]);
      xcSmallBall2D (*(iwksp + i));
      glPopMatrix();

      if ( VPf.labelsOn == 1 ) {
	makeAtomLabel2D ( *(iwksp + i), CENTERED_LABEL, CENTERED_LABEL );
      }
    }
    else if ( (coor + *(iwksp + i))->flag == SELATOM ) {
      glPushMatrix();
      glTranslated( (coor + *(iwksp + i))->x1,
		    (coor + *(iwksp + i))->y1, 
		    (coor + *(iwksp + i))->z1 );
      glColor3fv( atcol[SELCOL] );
      xcSmallBall2D (*(iwksp + i));
      /*glCallList( (coor + *(iwksp + i))->list2 );*/
      glPopMatrix();

      if ( VPf.labelsOn == 1 ) {
	makeAtomLabel2D ( *(iwksp + i), CENTERED_LABEL, CENTERED_LABEL );
      }
    }
    else if ( (coor + *(iwksp + i))->flag == BOND ) {	  
      glLineWidth( VPf.OUTlinewidth );
      /* MakeSticks2(*(iwksp + i), (coor + *(iwksp + i))->nat); */
      MakeSticks2(*(iwksp + i), (coor + *(iwksp + i))->sqn, PIPEBALL);
    }
    else if ( (coor + *(iwksp + i))->flag == SELBOND ) {
      glLineWidth( VPf.OUTlinewidth );
      MakeSticks2(*(iwksp + i), 0, PIPEBALL);
    }
    else if ( (coor + *(iwksp + i))->flag == SELLINE ) {
      glLineStipple( 2, SELLINE_PAT );
      glEnable(GL_LINE_STIPPLE);
      glLineWidth( SELLINE_WIDTH ); 
      glColor3fv( atcol[SELCOL] );	  
      glBegin(GL_LINES);
        glVertex3d( (coor + *(iwksp + i))->x1, 
		    (coor + *(iwksp + i))->y1,
		    (coor + *(iwksp + i))->z1 );
	glVertex3d( (coor + *(iwksp + i))->x2, 
		    (coor + *(iwksp + i))->y2,
		    (coor + *(iwksp + i))->z2 );
      glEnd(); 
      glDisable(GL_LINE_STIPPLE);
    }
    else if ( (coor + *(iwksp + i))->flag == FRAME ) {
      if ( (coor + *(iwksp + i))->nat == 2 ) {
	glLineStipple(2, FRAME_PAT);
	glEnable(GL_LINE_STIPPLE);
      }
      glLineWidth( (GLfloat) VPf.framewidth );  
      glBegin(GL_LINES);
        glColor3fv( framecol );
	glVertex3d( (coor + *(iwksp + i))->x1, 
		    (coor + *(iwksp + i))->y1,
		    (coor + *(iwksp + i))->z1 );
	glVertex3d( (coor + *(iwksp + i))->x2, 
		    (coor + *(iwksp + i))->y2,
		    (coor + *(iwksp + i))->z2 );
      glEnd(); 
      glDisable(GL_LINE_STIPPLE);
    }	
  }
  glDisable( GL_CULL_FACE );

  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}

/* xcBallStick22D is xcDisplayFunc for Ballstick2, 
 * -----------------------------------------------
 *             xcBallStick12D display a BallSticks2
 *            Here I'm dealing with 'coor' structure
 */
void
xcBallStick22D(struct Togl *togl)
{
  register int i;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glLineWidth( VPf.OUTlinewidth );
  glEnable( GL_CULL_FACE );
  glCullFace( GL_BACK );
  for(i = 1; i <= tmp_nobjects; i++) {
    if ( (coor + *(iwksp + i))->flag == ATOM ) {
      glPushMatrix();
      glTranslated( (coor + *(iwksp + i))->x1,
		    (coor + *(iwksp + i))->y1, 
		    (coor + *(iwksp + i))->z1 );
      glColor3fv (atm.col[(coor + *(iwksp + i))->sqn]);
      xcBall2D (*(iwksp + i));
      glPopMatrix();
    }
    else if ( (coor + *(iwksp + i))->flag == BOND ) {	  
      /* MakeSticks2(*(iwksp + i), (coor + *(iwksp + i))->nat); */
      MakeSticks2(*(iwksp + i), (coor + *(iwksp + i))->sqn, BALL);
    }
  }
  glDisable( GL_CULL_FACE );

  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void xcSpaceFill2D(struct Togl *togl)
{
  register int i;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glLineWidth( VPf.OUTlinewidth );
  glEnable( GL_CULL_FACE );
  glCullFace( GL_BACK );
  for(i = 1; i <= tmp_nobjects; i++) {
    if ( (coor + *(iwksp + i))->flag == ATOM ) {
      glLineWidth( VPf.OUTlinewidth );
      glPushMatrix();
      glTranslated( (coor + *(iwksp + i))->x1,
		    (coor + *(iwksp + i))->y1, 
		    (coor + *(iwksp + i))->z1 );
      glColor3fv (atm.col[(coor + *(iwksp + i))->sqn]);
      xcBigBall2D (*(iwksp + i));
      glPopMatrix();
      if ( VPf.labelsOn == 1 ) {
	makeAtomLabel2D ( *(iwksp + i), CENTERED_LABEL, CENTERED_LABEL );
      }
    }
    else if ( (coor + *(iwksp + i))->flag == SELATOM ) {
      glPushMatrix();
      glTranslated( (coor + *(iwksp + i))->x1,
		    (coor + *(iwksp + i))->y1, 
		    (coor + *(iwksp + i))->z1 );
      glColor3fv( atcol[SELCOL] );
      xcBigBall2D (*(iwksp + i));
      glPopMatrix();

      if ( VPf.labelsOn == 1 ) {
	makeAtomLabel2D ( *(iwksp + i), CENTERED_LABEL, CENTERED_LABEL );
      }
    }
    else if ( (coor + *(iwksp + i))->flag == SELLINE ) {
      glLineStipple( 2, SELLINE_PAT );
      glEnable(GL_LINE_STIPPLE);
      glLineWidth( SELLINE_WIDTH*2 ); 
      glColor3fv( atcol[SELCOL] );	  
      glBegin(GL_LINES);
        glVertex3d( (coor + *(iwksp + i))->x1, 
		    (coor + *(iwksp + i))->y1,
		    (coor + *(iwksp + i))->z1 );
	glVertex3d( (coor + *(iwksp + i))->x2, 
		    (coor + *(iwksp + i))->y2,
		    (coor + *(iwksp + i))->z2 );
      glEnd(); 
      glDisable(GL_LINE_STIPPLE);
    }
    else if ( (coor + *(iwksp + i))->flag == FRAME ) {
      if ( (coor + *(iwksp + i))->nat == 2 ) {
	glLineStipple(2, FRAME_PAT);
	glEnable(GL_LINE_STIPPLE);
      }
      glLineWidth( (GLfloat) VPf.framewidth );  
      glBegin(GL_LINES);
        glColor3fv( framecol );
	glVertex3d( (coor + *(iwksp + i))->x1, 
		    (coor + *(iwksp + i))->y1,
		    (coor + *(iwksp + i))->z1 );
	glVertex3d( (coor + *(iwksp + i))->x2, 
		    (coor + *(iwksp + i))->y2, 
		    (coor + *(iwksp + i))->z2 );
      glEnd(); 
      glDisable(GL_LINE_STIPPLE);
    }	
  }
  glDisable( GL_CULL_FACE );

  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


/*****************************************************************************/
/* Do we have 2D or 3D ?????? */
void 
EnableOr2D_Or3D(void)
{
  register int i;
  GLenum LIGHT[8] = { GL_LIGHT0,
		      GL_LIGHT1,
		      GL_LIGHT2,
		      GL_LIGHT3,
		      GL_LIGHT4,
		      GL_LIGHT5,
		      GL_LIGHT6,
		      GL_LIGHT7 };

  if ( dimType == XC_2D ) 
    {
      glDisable( GL_DEPTH_TEST );      
      glDisable( GL_LIGHTING );
      glDisable( GL_DITHER );
      glShadeModel( GL_FLAT );
      for (i=0; i<MAX_NUMBER_OF_LIGHTS; i++) 
	glDisable( LIGHT[i] );
    }
  else if (dimType == XC_3D )
    {
      /*printf("XC_3D\n",NULL);*/
      if (is.smooth) glShadeModel( GL_SMOOTH );
      else glShadeModel( GL_FLAT );
      glEnable( GL_DITHER );
      glEnable( GL_DEPTH_TEST );      
      glEnable( GL_LIGHTING );
      glEnable( GL_LIGHT0 );
      for (i=0; i<MAX_NUMBER_OF_LIGHTS; i++) 
	if ( light[i].isenabled )
	  glEnable( LIGHT[i] );
 
      /* turn-off 2D antialiasing */
      xcAntiAlias (GL_FALSE);
    }
}


/********************************************  
 *   xcPrimitives (sphere, cylinder, bonds)
 ********************************************/
void 
xcWireSphere (GLdouble radius)
{
  GLint nstep;
  GLUquadricObj *quadObj;
  GLdouble      sizeArray[3];
  GLuint        displayList;
  
  nstep = SPHERE_TESS_PREFACTOR * CalcTessFactor();
  if ( tr.b1motion || tr.b2motion || tr.shiftB1motion ) nstep /= 4.0;

  /* nstep = 4 * (GLint) ( radius * VPf.VPfactor / tessFactor ); */
  /*nstep = (GLint) (2 * APPR_SCREENSIZE / (tessFactor * ort.size));*/
  
  /* nstep shouldn't be lower than 8 and greater than XX */
  if ( nstep < 6 ) nstep = 6; 
  /* if ( nstep > 16 ) nstep = 16; */

  sizeArray[0]  = radius;
  sizeArray[1]  = (double) nstep;
  sizeArray[2]  = (double) VPf.WF3Dlinewidth;
  displayList   = findList1 (SPHEREWIRE, sizeArray, 3);

  if ( displayList == 0 ) {
    glNewList( makeModelPtr1 (SPHEREWIRE, sizeArray, 3),
	       GL_COMPILE_AND_EXECUTE );
      glLineWidth (VPf.WF3Dlinewidth);
      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_LINE);
      gluQuadricNormals (quadObj, GLU_SMOOTH);
      gluSphere (quadObj, radius, 2*nstep, nstep);
    glEndList();
  } else {
    glCallList(displayList);
  }
}


void 
xcSolidSphere (GLdouble radius)
{
  GLint         nstep;
  GLUquadricObj *quadObj;
  GLdouble      sizeArray[2];
  GLuint        displayList;

  nstep = SPHERE_TESS_PREFACTOR * CalcTessFactor();
  if ( tr.b1motion || tr.b2motion || tr.shiftB1motion ) nstep /= 4.0;

  /*if ( tr.b1motion || tr.b2motion || tr.shiftB1motion ) tessFactor *= 4.0;*/
  /* nstep = 4 * (GLint) ( radius * VPf.VPfactor / tessFactor ); */
  /*nstep = (GLint) (2 * APPR_SCREENSIZE / (tessFactor * ort.size));*/

  /* nstep shouldn't be lower than 8 and gretaer than XX */
  if ( nstep < 6 )  nstep = 6; 
  /*if ( nstep > 16 ) nstep = 16;*/

  sizeArray[0]	= radius;
  sizeArray[1]  = (double) nstep;
  displayList   = findList1 (SPHERESOLID, sizeArray, 2);

  if ( displayList == 0 ) {
    glNewList( makeModelPtr1 (SPHERESOLID, sizeArray, 2),
	       GL_COMPILE_AND_EXECUTE );			      
      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_FILL);
      gluQuadricNormals (quadObj, GLU_SMOOTH);
      gluSphere (quadObj, radius, 2*nstep, nstep);
    glEndList();
  } else {
    glCallList(displayList);
  }
}


void 
xcWireCylinder (GLdouble radius, GLdouble height)
{
  GLint nstep;
  GLUquadricObj *quadObj;
  GLdouble      sizeArray[4];
  GLuint        displayList;

  nstep = (int) (BOND_TESS_PREFACTOR * (float)CalcTessFactor());
  if ( tr.b1motion || tr.b2motion || tr.shiftB1motion ) nstep /= 4.0;
  /*if ( tr.b1motion || tr.b2motion || tr.shiftB1motion ) tessFactor *= 4.0;*/
  /* nstep = 4 * (GLint) ( radius * VPf.VPfactor / tessFactor ); */
  /*nstep = (GLint) (1.3 * APPR_SCREENSIZE / (tessFactor * ort.size));*/

  /* nstep shouldn't be lower than 6 and gretaer than XX */
  if ( nstep < 6 ) nstep = 6; 
  /*if ( nstep > 10 ) nstep = 10;*/

  sizeArray[0]  = radius;
  sizeArray[1]  = height;
  sizeArray[2]  = (double) nstep;
  sizeArray[3]  = (double) VPf.WF3Dlinewidth;
  displayList   = findList1 (CYLINDERWIRE, sizeArray, 4);
  
  if ( displayList == 0 ) {
    glNewList( makeModelPtr1 (CYLINDERWIRE, sizeArray, 4),
	       GL_COMPILE_AND_EXECUTE );
      glLineWidth (VPf.WF3Dlinewidth);			      
      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_LINE);
      gluQuadricNormals (quadObj, GLU_SMOOTH);
      gluCylinder (quadObj, radius, radius, height, nstep, 1);
    glEndList();
  } else {
    glCallList(displayList);
  }
}

void 
xcSolidCylinder (GLdouble radius, GLdouble height)
{
  GLint nstep;
  GLUquadricObj *quadObj;
  GLdouble      sizeArray[3];
  GLuint        displayList;

  nstep = (int) (BOND_TESS_PREFACTOR * (float)CalcTessFactor());
  if ( tr.b1motion || tr.b2motion || tr.shiftB1motion ) nstep /= 4.0;   
  /* nstep = 4 * (GLint) ( radius * VPf.VPfactor / tessFactor ); */
  /*nstep = (GLint) (1.3 * APPR_SCREENSIZE / (tessFactor * ort.size));*/

  /* nstep shouldn't be lower than 6 and gretaer than XX */
  if ( nstep < 6 ) nstep = 6; 
  /* if ( nstep > 10 ) nstep = 10; */

  sizeArray[0]  = radius;
  sizeArray[1]  = height;
  sizeArray[2]  = (double) nstep;
  displayList   = findList1 (CYLINDERSOLID, sizeArray, 3);

  if ( displayList == 0 ) {
    glNewList( makeModelPtr1 (CYLINDERSOLID, sizeArray, 3),
	       GL_COMPILE_AND_EXECUTE );		      
      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_FILL);
      gluQuadricNormals (quadObj, GLU_SMOOTH);
      gluCylinder (quadObj, radius, radius, height, nstep, 1);
 
      /* lower-end disk */
      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_FILL);
      gluQuadricNormals (quadObj, GLU_FLAT);
      gluQuadricOrientation( quadObj, GLU_INSIDE );
      gluDisk( quadObj, 0.0, radius, nstep, 1 );

      /* upper end disk */
      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_FILL);
      gluQuadricNormals (quadObj, GLU_FLAT);
      gluQuadricOrientation( quadObj, GLU_OUTSIDE );
      glPushMatrix();
        glTranslated(0.0, 0.0, height);
	gluDisk( quadObj, 0.0, radius, nstep, 1 );
      glPopMatrix();
    glEndList();
  } else {
    glCallList(displayList);
  }
}


void 
xcSolidBond (GLdouble radius, GLdouble height, int bondFlag)
{
  GLint nstep;
  GLUquadricObj *quadObj;
  GLdouble      sizeArray[4];
  GLuint        displayList;

  
  nstep = (int) (BOND_TESS_PREFACTOR * (float) CalcTessFactor());
  if ( tr.b1motion || tr.b2motion || tr.shiftB1motion ) nstep /= 4.0;

  /* nstep = 4 * (GLint) ( radius * VPf.VPfactor / tessFactor ); */
  /*nstep = (GLint) (1.3 * APPR_SCREENSIZE / (tessFactor * ort.size));*/

  /* nstep shouldn't be lower than 6 and gretaer than XX */
  if ( nstep < 6 ) nstep = 6;

  sizeArray[0]  = radius;
  sizeArray[1]  = height;
  sizeArray[2]  = (double) nstep;
  sizeArray[3]  = (double) bondFlag;
  displayList   = findList1 (BONDSOLID, sizeArray, 4);

  if ( displayList == 0 ) {
    glNewList( makeModelPtr1 (BONDSOLID, sizeArray, 4),
	       GL_COMPILE_AND_EXECUTE );		      
      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_FILL);
      gluQuadricNormals (quadObj, GLU_SMOOTH);
      gluCylinder (quadObj, radius, radius, height, nstep, 1);

      if ( (bondFlag == BOND_ATOM_TO_MIDBOND) ||
	   (bondFlag == BOND_ATOM_TO_ATOM) ) {
	/* lower-end sphere */
	quadObj = gluNewQuadric ();
	gluQuadricDrawStyle (quadObj, GLU_FILL);
	gluQuadricNormals (quadObj, GLU_SMOOTH);
	gluQuadricOrientation( quadObj, GLU_OUTSIDE );
	gluSphere (quadObj, radius, nstep, nstep);
      }

      if ( (bondFlag == BOND_ATOM_TO_ATOM ) ||
	   (bondFlag == BOND_MIDBOND_TO_ATOM) ) {
	/* upper-end sphere */
	quadObj = gluNewQuadric ();
	gluQuadricDrawStyle (quadObj, GLU_FILL);
	gluQuadricNormals (quadObj, GLU_SMOOTH);
	gluQuadricOrientation( quadObj, GLU_OUTSIDE );
	glPushMatrix();
          glTranslated(0.0, 0.0, height);
          gluSphere (quadObj, radius, nstep, nstep);
	glPopMatrix();
      }
    glEndList();
  } else {
    glCallList(displayList);
  }   
}


/*  Render solid Cone.  
 */
void 
xcSolidCone (GLdouble baseradius, GLdouble topradius, GLdouble height)
{
  GLint nstep;
  GLUquadricObj *quadObj;
  GLdouble      sizeArray[4];
  GLuint        displayList;
  /*double tessFactor = VPf.tessFactor;*/

  nstep = (int) (BOND_TESS_PREFACTOR * (float) CalcTessFactor());
  if ( tr.b1motion || tr.b2motion || tr.shiftB1motion ) nstep /= 4.0;

  /*if ( tr.b1motion || tr.b2motion || tr.shiftB1motion ) tessFactor *= 4.0;*/
  /* nstep = 4 * (GLint) ( baseradius * VPf.VPfactor / tessFactor ); */
  /*nstep = (GLint) (1.3 * APPR_SCREENSIZE / (tessFactor * ort.size));*/

  /* nstep shouldn't be lower than 6 and gretaer than XX */
  if ( nstep < 6 ) nstep = 6; 

  sizeArray[0]  = baseradius;
  sizeArray[1]  = topradius;
  sizeArray[2]  = height;
  sizeArray[3]  = (double) nstep;
  displayList   = findList1 (CONESOLID, sizeArray, 4);

  if ( displayList == 0 ) {
    glNewList( makeModelPtr1 (CONESOLID, sizeArray, 4),
	       GL_COMPILE_AND_EXECUTE );		      
      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_FILL);
      gluQuadricNormals (quadObj, GLU_SMOOTH);
      gluCylinder (quadObj, baseradius, topradius, height, nstep, 1);

      /* lower end disk */
      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_FILL);
      gluQuadricNormals (quadObj, GLU_FLAT);
      gluQuadricOrientation( quadObj, GLU_INSIDE );
      gluDisk( quadObj, 0.0, baseradius, nstep, 1 );
    glEndList();
  } else {
    glCallList(displayList);
  }   
}

static void _ambient_by_diffuse(GLfloat dst[4], GLfloat scale, GLfloat src[4]) {
  COPY_AND_SCALE_V(4, dst, scale, src);
  CLAMP_V(4, 1.0, dst);
}

#define WIRE    GL_FALSE
#define SOLID   GL_TRUE
void
xcRenderSolidBalls3D(void) {
  register int i;
  GLfloat ambient[4];

  for(i = 1; i <= natoms; i++) {
    glPushMatrix();
      glTranslated( *(xat + i), *(yat + i), *(zat + i) );
      if ( AtCol_Ambient_by_Diffuse == 1.0) 
	glMaterialfv ( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, atm.col[i] );
      else {
	_ambient_by_diffuse(ambient, AtCol_Ambient_by_Diffuse, atm.col[i]);
	glMaterialfv (GL_FRONT, GL_AMBIENT, ambient);
	glMaterialfv (GL_FRONT, GL_DIFFUSE, atm.col[i]);
      }

      /* new */
      /*glMaterialfv( GL_FRONT, GL_SPECULAR, atm.col[i]);*/
      if (is.pipemode) {
	xcSolidSphere( rball[1] );
      } else {
	xcSolidSphere( rball[ *(nat + i) ] );
      }
    glPopMatrix();
  }
}  


void
xcRenderWireBalls3D(void) {
  register int i;
  GLfloat ambient[4];

  for(i = 1; i <= natoms; i++) {
    glPushMatrix();
      glTranslated( *(xat + i), *(yat + i), *(zat + i) );

      if ( AtCol_Ambient_by_Diffuse == 1.0) 
	glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, atm.col[i] );
      else {
	_ambient_by_diffuse(ambient, AtCol_Ambient_by_Diffuse, atm.col[i]);
	glMaterialfv (GL_FRONT, GL_AMBIENT, ambient);
	glMaterialfv (GL_FRONT, GL_DIFFUSE, atm.col[i]);
      }
      
      if (is.pipemode) {
	xcWireSphere( rball[1] );
      } else {
	xcWireSphere( rball[ *(nat + i) ] );
      }
    glPopMatrix();
  }
}


void
xcRenderBonds(GLenum type)
{
  register int i;
  GLfloat ambient[4];

  if ( VPf.unibond) {
    if ( AtCol_Ambient_by_Diffuse == 1.0) 
      glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, unibondCol);
    else {
      _ambient_by_diffuse(ambient, AtCol_Ambient_by_Diffuse, unibondCol);
      glMaterialfv (GL_FRONT, GL_AMBIENT, ambient);
      glMaterialfv (GL_FRONT, GL_DIFFUSE, unibondCol);
    }
  }
  
  if ( type == GL_FILL ) {
    /* hack for white bonds !!! */
    
    /*fprintf(stderr,"Bonds GL-coor:\n----------\n");*/
    for(i = 1; i <= nbonds; i++) {
      glPushMatrix();
      /*
        fprintf(stderr,"-----     %f %f %f   %f %f %f %f\n",
		(coor3D + i)->x1,
		(coor3D + i)->y1,
		(coor3D + i)->z1,
		(coor3D + i)->fibond,
		(coor3D + i)->xrvb,
		(coor3D + i)->yrvb,
		(coor3D + i)->zrvb );
      */
	glTranslated( (coor3D + i)->x1,
		      (coor3D + i)->y1,
		      (coor3D + i)->z1 );
	glRotated( (coor3D + i)->fibond,
		   (coor3D + i)->xrvb,
		   (coor3D + i)->yrvb,
		   (coor3D + i)->zrvb );
	if (!VPf.unibond) {
	  if ( AtCol_Ambient_by_Diffuse == 1.0)
	    glMaterialfv( GL_FRONT, 
			  GL_AMBIENT_AND_DIFFUSE, atm.col[(coor3D + i)->sqn] );
	  else {
	    _ambient_by_diffuse(ambient, AtCol_Ambient_by_Diffuse, atm.col[(coor3D + i)->sqn]);
	    glMaterialfv (GL_FRONT, GL_AMBIENT, ambient);
	    glMaterialfv (GL_FRONT, GL_DIFFUSE, atm.col[(coor3D + i)->sqn]);
	  }
	}
	xcSolidBond( rrod, (coor3D + i)->bondl, (coor3D + i)->bondFlag );
	/*xcSolidCylinder( rrod, (coor3D + i)->bondl );*/
	glPopMatrix();
    }
  } else {
    /* GL_LINE */    
    for(i = 1; i <= nbonds; i++) {
      glPushMatrix();
        glTranslated( (coor3D + i)->x1,
		      (coor3D + i)->y1,
		      (coor3D + i)->z1 );
	glRotated( (coor3D + i)->fibond,
		   (coor3D + i)->xrvb,
		   (coor3D + i)->yrvb,
		   (coor3D + i)->zrvb );
	if (!VPf.unibond) {
	  if ( AtCol_Ambient_by_Diffuse == 1.0)
	    glMaterialfv( GL_FRONT, 
			  GL_AMBIENT_AND_DIFFUSE, atm.col[(coor3D + i)->sqn] );
	  else {
	    _ambient_by_diffuse(ambient, AtCol_Ambient_by_Diffuse, atm.col[(coor3D + i)->sqn]);
	    glMaterialfv (GL_FRONT, GL_AMBIENT, ambient);
	    glMaterialfv (GL_FRONT, GL_DIFFUSE, atm.col[(coor3D + i)->sqn]);
	  }
	}
	xcWireCylinder( rrod, (coor3D + i)->bondl );
      glPopMatrix();
    }
  }
}


void
xcRenderSolidSpaceFills3D(void)
{
  register int i;
  GLfloat ambient[4];

  for(i = 1; i <= natoms; i++) {
    glPushMatrix();
      glTranslated( *(xat + i), *(yat + i), *(zat + i) );      
      if ( AtCol_Ambient_by_Diffuse == 1.0)
	glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, atm.col[i] );
      else {
	_ambient_by_diffuse(ambient, AtCol_Ambient_by_Diffuse, atm.col[i]);
	glMaterialfv (GL_FRONT, GL_AMBIENT, ambient);
	glMaterialfv (GL_FRONT, GL_DIFFUSE, atm.col[i]);
      }

      /* new */
      /*glMaterialfv( GL_FRONT, GL_SPECULAR, atm.col[i]);*/

      xcSolidSphere( atrad[ *(nat + i) ] );	  
    glPopMatrix();
  }
}


void
xcRenderWireSpaceFills3D(void)
{
  register int i;
  GLfloat ambient[4];

  for(i = 1; i <= natoms; i++) {
    glPushMatrix();
      glTranslated( *(xat + i), *(yat + i), *(zat + i) );
      
      if ( AtCol_Ambient_by_Diffuse == 1.0)
	glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, atm.col[i] );
      else {
	_ambient_by_diffuse(ambient, AtCol_Ambient_by_Diffuse, atm.col[i]);
	glMaterialfv (GL_FRONT, GL_AMBIENT, ambient);
	glMaterialfv (GL_FRONT, GL_DIFFUSE, atm.col[i]);
      }
      xcWireSphere( atrad[ *(nat + i) ] );	  
    glPopMatrix();
  }
}


void 
xcRenderSolidFrame3D(void) {
  register int i;
  GLfloat ambient[4];

  if ( AtCol_Ambient_by_Diffuse == 1.0)
    glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, framecol );
  else {
    _ambient_by_diffuse(ambient, AtCol_Ambient_by_Diffuse, framecol );
    glMaterialfv (GL_FRONT, GL_AMBIENT, ambient);
    glMaterialfv (GL_FRONT, GL_DIFFUSE, framecol);
  }
  
  for(i = nbonds + 1; i <= nframes + nbonds; i++) {
    /* if frametype == 1 -> draw a xcSolidBond frame else,
     * draw a line_stipple frame
     */
    if ( (coor3D + i)->nat == 1 ) {
      glPushMatrix();
        glTranslated( (coor3D + i)->x1,
		      (coor3D + i)->y1,
		      (coor3D + i)->z1 );
	glRotated( (coor3D + i)->fibond,
		   (coor3D + i)->xrvb,
		   (coor3D + i)->yrvb,
		   (coor3D + i)->zrvb );
	xcSolidBond ( rframe, (coor3D + i)->bondl, BOND_ATOM_TO_ATOM );
      glPopMatrix();
    } 
    else xcLineStippleFrame(i - nbonds);	
  }
}


void
xcRenderWireFrame3D(void) {
  register int i;
  GLfloat ambient[4];

  if ( AtCol_Ambient_by_Diffuse == 1.0)
    glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, framecol );
  else {
    _ambient_by_diffuse(ambient, AtCol_Ambient_by_Diffuse, framecol );
    glMaterialfv (GL_FRONT, GL_AMBIENT, ambient);
    glMaterialfv (GL_FRONT, GL_DIFFUSE, framecol);
  }

  for(i = nbonds + 1; i <= nframes + nbonds; i++) {
    if ( (coor3D + i)->nat == 1 ) {
      glPushMatrix();
        glTranslated( (coor3D + i)->x1,
		      (coor3D + i)->y1,
		      (coor3D + i)->z1 );
	glRotated( (coor3D + i)->fibond,
		   (coor3D + i)->xrvb,
		   (coor3D + i)->yrvb,
		   (coor3D + i)->zrvb );
	xcWireCylinder ( rframe, (coor3D + i)->bondl );
      glPopMatrix();
    } 
    else xcLineStippleFrame(i - nbonds);
  }
}


void
xcRenderLineFrame3D(void) {
  register int i;

  glColor3fv( framecol );
  /* temporary disable lighting, dihtering */
  glCallList( tempDisable3Dlist );
  for(i = 1; i <= nframes; i++) {
    if ( *(frametype + i) != 1 ) {
      glLineStipple(2, FRAME_PAT);
      glEnable(GL_LINE_STIPPLE);
    }
    glCallList(tempDisable3Dlist);
    glLineWidth( (GLfloat) VPf.framewidth );
    glBegin(GL_LINES);
      glVertex3d( *(xframe + i), *(yframe + i), *(zframe + i) );
      glVertex3d( *(xframe2 + i), *(yframe2 + i), *(zframe2 + i) );
    glEnd();
    glDisable(GL_LINE_STIPPLE);
    glCallList(tempEnable3Dlist);
  }
  /* now again enable lighting, dihtering */
  glCallList( tempEnable3Dlist );
  glEndList();
}    


static void 
xcLineStippleFrame(int ii)
{  
  glLineStipple(2, FRAME_PAT);
  glEnable(GL_LINE_STIPPLE);
  glCallList( tempDisable3Dlist );
  glColor3fv( framecol );
  glLineWidth( (GLfloat) VPf.framewidth );
  glBegin(GL_LINES);
    glVertex3d( *(xframe + ii), *(yframe + ii), *(zframe + ii) );
    glVertex3d( *(xframe2 + ii), *(yframe2 + ii), *(zframe2 + ii) );
  glEnd();
  glCallList( tempEnable3Dlist );
  glDisable(GL_LINE_STIPPLE);
}


static void
MatInv33_44(double m[4][4], double minv[4][4], double *det)
{
  double f0, f1, f2, idet;
  register int i, j;

  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      minv[i][j] = m[i][j];
  
  f0 = m[1][1]*m[2][2] - m[1][2]*m[2][1];
  f1 = m[1][2]*m[2][0] - m[1][0]*m[2][2];  
  f2 = m[1][0]*m[2][1] - m[1][1]*m[2][0];

  *det = m[0][0]*f0 + m[0][1]*f1 + m[0][2]*f2;

  if (*det < 1.0e-10) {
    /* determinant is very small; return the minv=m */
    return;
  }

  idet = 1.0 / *det;

  minv[0][0] = f0 * idet;
  minv[1][0] = f1 * idet;
  minv[2][0] = f2 * idet;

  f0 = m[0][0] * idet;
  f1 = m[0][1] * idet;
  f2 = m[0][2] * idet;

  minv[0][1] = f2*m[2][1] - f1*m[2][2];
  minv[1][1] = f0*m[2][2] - f2*m[2][0];
  minv[2][1] = f1*m[2][0] - f0*m[2][1];
  minv[0][2] = f1*m[1][2] - f2*m[1][1];
  minv[1][2] = f2*m[1][0] - f0*m[1][2];
  minv[2][2] = f0*m[1][1] - f1*m[1][0];
}


static void
MultMat44_33Vec3(double mat[4][4], double vec[3], double res[3])
{
  int i;
  /* mat is !!!COLUMN!!!, not row major mode (due to OGL) */
  for(i=0; i<3; i++)
    res[i] = mat[0][i] * vec[0] + mat[1][i] * vec[1] + mat[2][i] * vec[2];
}


static void _makeAtomLabels3D(double *rad) {
  register int i;
  double       dummy, minv[4][4], v[3], z[3], sizeF;
  
  sizeF = 2.0 * VPf.VPfactor;

  /* disable lighting & dihtering ... */
  glCallList( tempDisable3Dlist );
  
  for(i = 1; i <= natoms; i++) {
    
    z[2] = 1.01 * rad[*(nat + i)];
    
    if ( !globalAtomLabel.base && ! atomLabel[i].base ) {
      char    *label;
      GLfloat *color;
      /* 
	 render default raster font 
      */      
      if (atomLabel[i].base) {
	label = atomLabel[i].label;
	color = atomLabel[i].bright_color;
	if (atomLabel[i].tkfont) {
	  z[0] = (double) ( -Tk_TextWidth(atomLabel[i].tkfont, label, strlen(label)) ) / sizeF;
	} else {
	  z[0] = (double) (-atomLabel[i].width * (int)strlen(label)) / sizeF;
	}
	z[1] = (double) (-atomLabel[i].height) / sizeF;	
      } else {
	label = element[*(nat + i)];
	color = globalAtomLabel.bright_color;
	z[0]  = rf.w2 * (double)strlen(label); 
	z[1]  = rf.h2;
      }

      MatInv33_44((double (*)[4]) vec.crdvec, minv, &dummy);
      MultMat44_33Vec3(minv, z, v);
      glColor3fv (color);
      glRasterPos3d ( *(xat + i) + v[0],
		      *(yat + i) + v[1],
		      *(zat + i) + v[2] );

      if (atomLabel[i].base && atomLabel[i].do_display \
	  && ! do_not_display_atomlabel[i]) {
	/* custom label */
	glListBase(atomLabel[i].base);
	xcFont_PrintString (label);
      } else {
	/* default label */
	if (globalAtomLabel.do_display && !do_not_display_atomlabel[i])
	  glCallList( atomlabelOffset + *(nat + i) );
      }
    } else {
      /*
	render one of new Togl fonts 
      */
      
      if (atomLabel[i].base) {
	/* new */
	glListBase (atomLabel[i].base); /* t.k. */

	if (atomLabel[i].tkfont) {
	  z[0] = (double) ( -Tk_TextWidth(atomLabel[i].tkfont, 
					  atomLabel[i].label, 
					  strlen(atomLabel[i].label)) ) / sizeF;
	} else {
	  z[0] = (double) (-atomLabel[i].width * (int)strlen(atomLabel[i].label)) / sizeF;
	}
	z[1] = (double) (-atomLabel[i].height) / sizeF;	

	MatInv33_44((double (*)[4]) vec.crdvec, minv, &dummy);
	MultMat44_33Vec3(minv, z, v);
	glColor3fv (atomLabel[i].bright_color);	
	glRasterPos3d( *(xat + i) + v[0],
		       *(yat + i) + v[1],
		       *(zat + i) + v[2]);
	if (atomLabel[i].do_display && !do_not_display_atomlabel[i])
	  xcFont_PrintString( atomLabel[i].label );	
      } else {
	/* old */
	glListBase (globalAtomLabel.base); /* t.k. */	
	if (globalAtomLabel.tkfont) {
	  z[0] = (double) ( -Tk_TextWidth(globalAtomLabel.tkfont, 
					  element[*(nat + i)], 
					  strlen(element[*(nat + i)])) ) / sizeF;
	} else {
	  z[0] = ((double) -globalAtomLabel.width * strlen(element[*(nat + i)])) / sizeF;
	}

	z[1] = ((double) -globalAtomLabel.height) / sizeF;	
	z[2] = 1.01 * rad[*(nat + i)];
	MatInv33_44((double (*)[4]) vec.crdvec, minv, &dummy);
	MultMat44_33Vec3(minv, z, v);
	glColor3fv (globalAtomLabel.bright_color);
	glRasterPos3d( *(xat + i) + v[0],
		       *(yat + i) + v[1],
		       *(zat + i) + v[2]);
	if (globalAtomLabel.do_display && !do_not_display_atomlabel[i])
	  xcFont_PrintString( element[*(nat + i)] );
      }
    }
  }
  
  /* now again enable lighting, dihtering */
  glCallList( tempEnable3Dlist );
}


void
xcMakeBallLabel3D(void)
{
  _makeAtomLabels3D(rball);
}

void
xcMakeSpaceLabel3D(void)
{
  _makeAtomLabels3D(atrad);
}


/*****************************************************************************/
/* when I enter this procedure I just Made a Display, that means that 
 * everything must be made before that
 */
void
xcDisplay3D(struct Togl *togl)
{

  if ( ! VPf.antialias ) 
    {
      /* 
	 NO ANTIALIASING 
      */
      _xcDisplay3D(togl);
    } 
  else 
    {
      /* 
	 ANTIALIASING 
      */
      GLint viewport[4];
      GLfloat sx, sy;
      GLfloat scale, dx, dy;
      int min, max, count, i, j;
      enum {
	XORG, YORG, WID, HT
      }; 

      glGetIntegerv(GL_VIEWPORT, viewport);
      
      sx = 2 * ort.maxx / viewport[WID];
      sy = 2 * ort.maxy / viewport[WID];

      min    = -antialias.degree;
      max    = -min + 1;
      count  = -min + max;
      count *= count;

      /* uniform scaling, less than one pixel wide */
      scale = -antialias.offset / min;

      glClear(GL_ACCUM_BUFFER_BIT);

      for (j = min; j < max; j++) {
	for (i = min; i < max; i++) 
	  {
	    dx = sx * scale * i;
	    dy = sy * scale * j;
	    glMatrixMode(GL_PROJECTION);
	    glLoadIdentity();
	    
	    if ( ! VPf.perspective ) 
	      {
		glOrtho(ort.minx + dx, ort.maxx + dx, 
			ort.miny + dy, ort.maxy + dy, 
			ort.minz+ort.size, ort.maxz+ort.size);
	      } 
	    else 
	      {
		glFrustum (ort.minx + dx, ort.maxx + dx, 
			   ort.miny + dy, ort.maxy + dy, 
			   persp.near, persp.far);
	      }
	    glMatrixMode(GL_MODELVIEW);
	    _xcDisplay3D(togl);
	    glAccum(GL_ACCUM, 1.0 / (GLfloat)count);
	  }
      }
      glAccum(GL_RETURN, 1.0);
    }

  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


static void _xcDisplay3D(struct Togl *togl) {
  /* /\* testing ... *\/ */
  /*   GLboolean can_do_stereo; */
  /*   glGetBooleanv(GL_STEREO, &can_do_stereo); */
  /*   is.stereo = GL_TRUE; */
  /* testing-end */

  if (is.stereo)
    {
      /***************/
      /* STEREO MODE */
      /***************/

      /* draw left eye image */

      glDrawBuffer(GL_BACK_LEFT);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      /* INSERT HERE ROTATION MATRIX */
      glLoadIdentity();
      if (VPf.perspective) 
	glTranslated (0.0, 0.0, persp.shiftZ);
      else
	glTranslated(0.0, 0.0, -ort.size);
      glRotated(1.0, 0.0, 1.0, 0.0);
      glMultMatrixd (vec.crdvec);

      __xcDisplay3D(togl);

      /* draw right eye image */

      glDrawBuffer(GL_BACK_RIGHT);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      /* INSERT HERE ROTATION MATRIX */
      glLoadIdentity();
      if (VPf.perspective) 
	glTranslated (0.0, 0.0, persp.shiftZ);
      else
	glTranslated(0.0, 0.0, -ort.size);
      glRotated(-1.0, 0.0, 1.0, 0.0);
      glMultMatrixd (vec.crdvec);

      __xcDisplay3D(togl);
    }
  else
    {
      /*******************/
      /* non STEREO MODE */
      /*******************/

      /* draw single image */
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      /* INSERT HERE ROTATION MATRIX */
      glLoadIdentity();
      if (VPf.perspective) 
	glTranslated (0.0, 0.0, persp.shiftZ);
      else
	glTranslated(0.0, 0.0, -ort.size);
      glMultMatrixd (vec.crdvec);

      __xcDisplay3D(togl);
    }
}


static void __xcDisplay3D(struct Togl *togl) {
  int i;

  /* fog */
  xcFog (togl, VPf.fog, VPf.perspective);
  
  if (is.anaglyph)
    {
      /**************************/
      /* DISPLAY with ANAGLYPHs */
      /**************************/

      /* EV display anaglyph */
      /* anaglyph is defined as wire is.solid=GL_FALSE */    
      glEnable(GL_COLOR_MATERIAL);
        CameraHome();
        HandleDisplay(togl);
      glDisable(GL_COLOR_MATERIAL);
    }
  else
    {
      /*****************************/
      /* DISPLAY without ANAGLYPHs */
      /*****************************/

      /* display FORCES */
      if (VPf.force && xcr.lforce) {
	xcRenderVectorForces();
      }
      
      /* display H-bonds */
      if (VPf.Hbond) xcRenderHbonds3D();        

      /* if (is.anaglyph) xcTurnXYZOff(togl); */
      
      tmp_nobjects = 0;
      if (is.solid) {
	/* display SOLID Lists */
	if (is.spacefillmode) {
	  xcRenderSolidSpaceFills3D();
	  tmp_nobjects += natoms;
	}
	if (is.ballmode)      {
	  xcRenderSolidBalls3D();
	  tmp_nobjects += natoms;
	}
	if (is.stickmode)     {
	  xcRenderBonds(GL_FILL);
	  tmp_nobjects += nbonds;
	}
	if (!is.lineframe && VPf.framesOn) {
	  xcRenderSolidFrame3D();
	  tmp_nobjects += nframes;
	}
      } else {
	if (is.spacefillmode) {
	  xcRenderWireSpaceFills3D();
	  tmp_nobjects += natoms;
	}
	if (is.ballmode)      {
	  xcRenderWireBalls3D();
	  tmp_nobjects += natoms;
	}
	if (is.stickmode)     {
	  xcRenderBonds(GL_LINE);
	  /* glCallList(WireStickList);*/
	  tmp_nobjects += nbonds;
	}
	if (!is.lineframe && VPf.framesOn) {
	  xcRenderWireFrame3D();
	  tmp_nobjects += nframes;
	}
      }
  
      if (is.lineframe && VPf.framesOn) {
	/* set color & width for frames */
	glLineWidth( VPf.framewidth );      
	xcRenderLineFrame3D();
	tmp_nobjects += nframes;
      }
    
      for (i=0; i<MAX_ISOOBJECTS; i++) {
	/* THIS IS NEEDED TO RENDER ISOSURFACE */
	if ( VPf.isosurf[i] ) xcRenderIsosurf(i);
	/* THIS IS NEEDED TO RENDER COLORPLANE */
	if ( VPf.colorplane[i] ) xcRenderColorplane(i);
	/* THIS IS NEEDED TO RENDER ISOLINE */
	if ( VPf.isoline[i] ) xcRenderIsoLine2D(i);
      }
      
      /* IsoSpaceSelection */
      if ( VPf.isospacesel2D ) IsoSpaceSel_Parallelogram();
      if ( VPf.isospacesel3D ) IsoSpaceSel_3D();
      
      if ( VPf.labelsOn ) {
	if ( is.ballmode || is.stickmode ) xcMakeBallLabel3D();
	if ( is.spacefillmode ) xcMakeSpaceLabel3D();
      }    
      
      /* t.k: this is temporal:
       * VORONOI & WIGNER-SEITZ
       */
      if ( xcr.lprimwigner && VPf.wignerseitz ) xcRenderVoronoi();
      
      /*
       * SUPERCELL OPTION
       */
      if ( VPf.supercell) (*xcSuperCell)();
      
      /*
       * NEW SURFACE OPTIONS 
       */
      if ( VPf.nsurface > 0 ) xcRenderSurface();
    
      /* are we in selection mode ??? */
      if ( VPf.selection ) {
	xcRenderSelAtoms3D();
	xcRenderSelBonds3D();
	/*for(i=0; i<MAXSEL-1; i++) {*/
	  /*if ( glIsList(SphereSelList[i]) ) glCallList(SphereSelList[i]);*/
	  /*if ( glIsList(LineSelList[i]) ) glCallList(LineSelList[i]);*/
	/*}*/
	/*if ( glIsList(SphereSelList[MAXSEL-1]) ) glCallList(SphereSelList[3]);*/
      }
      
      /* are we in atomadd mode ??? */
      if ( VPf.atomadd ) {
	xcDisplayAddAtomBox();	
      }
      
      /* check if coordinate sistem should be displayed */
      if ( VPf.xyzOn ) xcDisplayXYZ();
    }
    
  glFlush();
      
  printf("   ----------------\n");
  printf("ReMakeStr:: tmp_nobjects  = %d\n",tmp_nobjects);
  printf("ReMakeStr:: nbonds        = %d\n",nbonds);
  printf("ReMakeStr:: natoms        = %d\n",natoms);
  printf("ReMakeStr:: nframes       = %d\n",nframes);
  printf("ReMakeStr:: nobjects      = %d\n",nobjects);
  printf("--------------------------------------------------\n\n");
  fflush(stdout);
}


/* function delete 3D lists */
/* within 3D modes a whole structure is a list, so when I close a structure 
 * and open a new one or when I modify some radius I must delete 3D List and
 * create a new ones;
 */
/*
void
xcMaybeDestroyLists(void)
{  
  int i;
  
  for(i=0; i<25; i++)
    lists[i] = NULL;
}
*/

void
xcMakeProjection2D(const char *mode)
{
  GLdouble size;
    
  FindMaxRad();
  /* 
     this is too much complicated and the drawback is that the size of 
     displayed structure changes when chaning mode --> simplify !
  */

  /* if ( strcmp(mode,"WF") == 0 || strcmp(mode,"PL") == 0 ) */
  /*   /\* 1.1 is just some additional offset *\/ */
  /*   size = (MVf.structsize + VPf.PLradius / VPf.VPfactor) * 1.1;  */
  /* else if ( mode[0] == 'B' || strcmp(mode,"PB") == 0 ) { */
  /*   size = (MVf.structsize + max.r * VPf.ballf) * 1.1; */
  /* } */
  /* else if ( strcmp(mode,"SF") == 0 ) { */
  /*   size = (MVf.structsize + max.r * VPf.atradf) * 1.1; */
  /* } */

  /* this is much simpler */
  /*size = MVf.structsize * (size + max.r) * 1.1;*/
  /*size = MVf.structsize * (1.0 + max.r) * 1.1;*/
  size = (MVf.structsize + 2.0*max.r) * 1.05;

  _xcMakeProjection(size);
}


void
xcMakeProjection3D(const char *mode)
{
  GLdouble size;

/*   GLenum LIGHT[8] = { GL_LIGHT0, */
/* 		      GL_LIGHT1, */
/* 		      GL_LIGHT2, */
/* 		      GL_LIGHT3, */
/* 		      GL_LIGHT4, */
/* 		      GL_LIGHT5, */
/* 		      GL_LIGHT6, */
/* 		      GL_LIGHT7 }; */

  FindMaxRad();
 /* 
     this is too much complicated and the drawback is that the size of 
     displayed structure changes when chaning mode --> simplify !
  */
  
  /* size = MVf.structsize; */
  /* if ( strcmp(mode,"sticks") == 0 || strcmp(mode,"balls") == 0 ) */
  /*   size = (size + max.r * VPf.ballf) * 1.1; */
  /* else if ( strcmp(mode,"space") == 0 ) */
  /*   size = (size + max.r) * 1.1; */

  /* this is much simpler */
  /*size = MVf.structsize * (size + max.r) * 1.1;*/
  /*size = MVf.structsize * (1.0 + max.r) * 1.1;*/
  size = (MVf.structsize + 2.0*max.r) * 1.05;
  
  _xcMakeProjection(size);
}
 
static void _xcMakeProjection(GLdouble size) {
  register int i;
  int isoobject = 0;
  MOL_SURF     *m;
  ISOSURFACE   *iso;
  float _size;
  double structsize, fulls_vs_structs, ratio;
  static double old_fulls_vs_structs = 1.0;

  structsize = size;
  
  /* if isoXXXXX is rendered, than it may be greater than the structure,
   * take this into account
   */
  for (i=ISOOBJ_BASE; i<MAX_ISOOBJECTS; i++)
    if ( VPf.isosurf[i] || VPf.colorplane[i] || VPf.isoline[i] ) {
      isoobject = 1;
      break;
    }
  if ( isoobject && MVf.isosize > size ) size = MVf.isosize;

  /* are we in Iso_Space_Selection mode ??? */
  if ( (VPf.isospacesel2D || VPf.isospacesel3D) && 
       MVf.isospaceselsize > size ) 
    size = MVf.isospaceselsize;

  /* do we specify to render Wigner-Seitz cells ??? */
  if ( VPf.wignerseitz && MVf.wignerseitzsize > size )
    size = MVf.wignerseitzsize;

  /* this is temporal
     if ( xcr.lprimwigner ) size = MVf.structsize + wsp.max;
  */
  
  /* size of forces */
 
  size = FV.max_size > size ? FV.max_size : size;

  /* size of molecular surfaces */

  for (i=0; i<VPf.nsurface; i++) 
    {
      if (!VPf.surface[i]) continue;
      
      m   = (MOL_SURF *) VPf.surfacePtr[i];
      iso = FindIsoSurf( m->isosurf_index );
      if (!iso) {
	continue;
      }
      _size = distfv (m->size);
      size = _size > size ? _size : size;
    }

  if ( VPf.stropened && structsize > 0.0 ) {
    fulls_vs_structs     = size / structsize;
    ratio                = fulls_vs_structs / old_fulls_vs_structs;
    old_fulls_vs_structs = fulls_vs_structs;
    tr.zoom *= ratio;
  }

  /* assign glOrtho dimensions to OrhoProj ort */
  ort.size = (size > 1.0 ? size : 1.0);

  /* MaybeClipAndMakeProjection(); */
  xcViewPort();

  /* 
   * take care about some minimium light-settings
   */
/*   size = ort.maxx; */
/*   if ( ort.maxy > size ) */
/*     size = ort.maxy; */
  
/*   for (i=0; i< MAX_NUMBER_OF_LIGHTS; i++) { */
/*     if ( light[i].isenabled ) { */
/*       light[i].position[0] = size * light[i].fract_position[0]; */
/*       light[i].position[1] = size * light[i].fract_position[1]; */
/*       light[i].position[2] = size * light[i].fract_position[2]; */
/*       light[i].position[3] = light[i].fract_position[3]; */
/*       glLightfv(LIGHT[i], GL_POSITION, light[i].position); */
/*       glEnable( LIGHT[i] ); */
/*     } */
/*   } */

  /* incoming: */
  LoadLights();
}


void
xcClearScreen(struct Togl *togl)
{
  /* make a black polygon on whole window -> so structure will disapear */
  glClear(GL_COLOR_BUFFER_BIT);
  glClearColor( bg[0], bg[1], bg[2], bg[3] );
  glFlush();
  Togl_SwapBuffers( togl );
}


void
xcChangeBackground(struct Togl *togl, GLclampf bckg[4])
{
  glClear(GL_COLOR_BUFFER_BIT);
  glClearColor( bckg[0], bckg[1], bckg[2], bckg[3] );
  glFlush();
  Togl_PostRedisplay( togl );
}


void
xcDisplayXYZ(void)
{
  GLint    shadetype, fog; /* we will query shade model */
  GLdouble shift, f;  

  /* disable Depth-cuing */
  glGetIntegerv( GL_FOG, &fog);
  if (fog) glDisable (GL_FOG);

  /*glLoadMatrixd( vec.crdvec );*/
  glViewport( 0, 0, CRDS_SIZE, CRDS_SIZE);

  /* prepare for matrices for coordinate sistem */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (VPf.perspective) {    
    /* PERSPECTIVE projection */
    f = 1.5 + 2.4 / persp.shiftZ;
    shift = -2.0*persp.shiftZ;
  } else {
    /* ORTHOGRAPHIC projection */
    f = 1.5 + 2.4 / ort.size;
    shift = f*ort.size;
  }
  glOrtho( -1.2, 1.2, -1.2, 1.2, -shift-1.2, 1.2 );


  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslated (0.0, 0.0, shift);
  glMultMatrixd (vec.crdvec );
          
  /* ========================================== */
  /*        now display coor sistem             */
  
  /* if we are in XC_2D mode clear DETPTH_BUFFER */
  if ( dimType == XC_2D ) {
    glEnable(GL_DEPTH_TEST);  
    glClear( GL_DEPTH_BUFFER_BIT );
  }
  /* if we are in 3D turn off */
  if ( dimType == XC_3D ) {
    /* when displaying coord-sist in 3D mode, 
     * we must Disable some things
     */
    glCallList( tempDisable3Dlist );   

    glGetIntegerv( GL_SHADE_MODEL, &shadetype);
 
    /* go to FLAT shades, but before remember what was the shade model */
    if ( shadetype == GL_SMOOTH ) glShadeModel( GL_FLAT);      
  }

  /* display XYZ axis */
  glColor4fv( xyz.axis_color );
  glCallList (crdList);  
  
  glColor4fv( xyz.xyplane_color );
  glCallList (xyplaneList);

  /*glLoadIdentity();*/

  if ( dimType == XC_2D ) glDisable(GL_DEPTH_TEST);  
  if ( dimType == XC_3D ) {
    /* it's better idea to make a glList for this, but so far ... */
    glCallList( tempEnable3Dlist );
    if ( shadetype == GL_SMOOTH )
      glShadeModel( GL_SMOOTH );
  }  
  if (fog) glEnable (GL_FOG);

  goToPrevProj();
}


static void
goToPrevProj(void)
{
  xcViewPort();

  /* now go back to previous GL_PROJECTION & GL_MODELVIEW */
/*   glMatrixMode(GL_PROJECTION); */
/*   glLoadIdentity(); */
/*   if ( dimType == XC_2D ) { */
/*     gluOrtho2D(ort.minx, ort.maxx, ort.miny, ort.maxy); */
/*   } else  */
/*     glOrtho(ort.minx, ort.maxx, ort.miny, ort.maxy, ort.minz, ort.maxz); */
/*   glMatrixMode(GL_MODELVIEW); */
/*   glViewport( (GLint) VPf.x, (GLint) VPf.y,  */
/* 	      (GLsizei) VPf.sizex, (GLsizei) VPf.sizey); */
/*   if ( dimType == XC_2D ) glLoadIdentity();  */
}
    
  
void
xcTurnXYZOff(struct Togl *togl)
{
  /* turn VPf.xyzOn off */
  VPf.xyzOn = 0;
  Togl_PostRedisplay (togl);

/*   if ( VPf.stropened ) { */
/*   GLfloat size = 1.1 * CRDS_SIZE; */
/*     /\* on place where the coor-sist is desplayed, display black rectangle *\/ */
/*     /\* prepere for matrices for coordinate sistem *\/ */
/*     glMatrixMode(GL_PROJECTION); */
/*     glOrtho( -1.2, 1.2, -1.2, 1.2, -1.2, 1.2 ); */
/*     glMatrixMode(GL_MODELVIEW); */
/*     glViewport( 0.0, 0.0, size, size); */
/*     glLoadIdentity(); */
    
/*     glColor3fv( bg ); */
/*     /\* if we are in 3D mode we must go to "pseudo" 2D mode *\/ */
/*     if (dimType == XC_3D) { */
/*       glCallList(tempDisable3Dlist); */
/*       glDisable(GL_DEPTH_TEST); */
/*     } */
/*     glBegin(GL_POLYGON); */
/*       glVertex2f(-size, size); */
/*       glVertex2f( size, size); */
/*       glVertex2f( size,-size); */
/*       glVertex2f(-size,-size); */
/*     glEnd(); */
/*     /\* if 3D -> going back *\/ */
/*     if (dimType == XC_3D) { */
/*       glCallList(tempEnable3Dlist); */
/*       glEnable(GL_DEPTH_TEST); */
/*     } */
/*     goToPrevProj(); */
/*     glFlush(); */
/*     Togl_SwapBuffers(togl); */
/*   }   */
}


static int CalcTessFactor(void) {
  if ( natoms < 10 ) {
    return (int) VPf.tessFactor;
  } else if ( natoms < 30 ) {
    return (int) (VPf.tessFactor / 1.5);
  } else if ( natoms < 100 ) {
    return (int ) (VPf.tessFactor / 2);
  } else {
    return (int) (VPf.tessFactor / 3);
  }
}




/********************************************************************
 ********************************************************************
                              EV   ANAGLYPH PART
 ********************************************************************
 ********************************************************************/

/*
  This is the basic display callback routine
  It creates the geometry, lighting, and viewing position
*/

static void 
HandleDisplay(struct Togl *togl)
{
  XYZ right,focus;
   
  /* Determine the focal point */
  Normalise(&camera.vd);
  focus.x = camera.vp.x + camera.focallength * camera.vd.x;
  focus.y = camera.vp.y + camera.focallength * camera.vd.y;
  focus.z = camera.vp.z + camera.focallength * camera.vd.z;

  /* Derive the the "right" vector */
  CROSSPROD(camera.vd,camera.vu,right);
  Normalise(&right);
  right.x *= camera.eyesep / 2.0;
  right.y *= camera.eyesep / 2.0;
  right.z *= camera.eyesep / 2.0;

  /* Set the buffer for writing and reading */
  glDrawBuffer(GL_BACK);
  glReadBuffer(GL_BACK);
   
  /* Clear things */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  
   
   
  /* Set projection */
  xcViewPort();
   
   
  /* Left eye filter red */
  glColorMask(GL_TRUE,GL_FALSE,GL_FALSE,GL_TRUE);
   
   
  /* Create the model */
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if (VPf.perspective)
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);
  glMultMatrixd (vec.crdvec);
   
  gluLookAt(camera.vp.x - right.x,
	    camera.vp.y - right.y,
	    camera.vp.z - right.z,
	    focus.x,focus.y,focus.z,
	    camera.vu.x,camera.vu.y,camera.vu.z);
   
  draw_scene();
  glFlush();
   

  glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
   
  glDrawBuffer(GL_BACK);
   
  glClear(GL_DEPTH_BUFFER_BIT);
   
  xcViewPort();

  /* Right eye filter */
  
  glColorMask(GL_FALSE,GL_FALSE,GL_TRUE,GL_TRUE);
   
   
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if (VPf.perspective)
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);
  glMultMatrixd (vec.crdvec); 
  gluLookAt(camera.vp.x + right.x,
	    camera.vp.y + right.y,
	    camera.vp.z + right.z,
	    focus.x,focus.y,focus.z,
	    camera.vu.x,camera.vu.y,camera.vu.z);
   
  draw_scene();
  glFlush();
  glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
   
   
  /* Let's look at it */
  /*Togl_SwapBuffers(togl);*/
   
}


static void
draw_scene() 
{
  int i;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix();
  /*  Here are the options set in _xcDisplay3D for wire or solid */ 
  /* display FORCES */
  if (VPf.force && xcr.lforce) {
    xcRenderVectorForces();
  }
  
  /* display H-bonds */
  if (VPf.Hbond) xcRenderHbonds3D();

  tmp_nobjects = 0;

  if (is.solid) 
    {
      /* display SOLID Lists */
      if (is.spacefillmode) {
	xcRenderSolidSpaceFills3D();
	tmp_nobjects += natoms;
      }
      if (is.ballmode)      {
	xcRenderSolidBalls3D();
	tmp_nobjects += natoms;
      }
      if (is.stickmode)     {
	xcRenderBonds(GL_FILL);
	tmp_nobjects += nbonds;
      }
      if (!is.lineframe && VPf.framesOn) {
	xcRenderSolidFrame3D();
	tmp_nobjects += nframes;
      }
    } 
  else 
    {
      if (is.spacefillmode) {
	xcRenderWireSpaceFills3D();
	tmp_nobjects += natoms;
      }
      if (is.ballmode)      {
	xcRenderWireBalls3D();
	tmp_nobjects += natoms;
      }
      if (is.stickmode)     {
	xcRenderBonds(GL_LINE);
	/* glCallList(WireStickList);*/
	tmp_nobjects += nbonds;
      }
      if (!is.lineframe && VPf.framesOn) {
	xcRenderWireFrame3D();
	tmp_nobjects += nframes;
      }
    }
  
  if (is.lineframe && VPf.framesOn) {
    /* set color & width for frames */
    glLineWidth( VPf.framewidth );      
    xcRenderLineFrame3D();
    tmp_nobjects += nframes;
  }
 
  for (i=0; i<MAX_ISOOBJECTS; i++) {
    /* THIS IS NEEDED TO RENDER ISOSURFACE */
    if ( VPf.isosurf[i] ) xcRenderIsosurf(i);
    /* THIS IS NEEDED TO RENDER COLORPLANE */
    if ( VPf.colorplane[i] ) xcRenderColorplane(i);
    /* THIS IS NEEDED TO RENDER ISOLINE */
    if ( VPf.isoline[i] ) xcRenderIsoLine2D(i);
  }

  /* IsoSpaceSelection */
  if ( VPf.isospacesel2D ) IsoSpaceSel_Parallelogram();
  if ( VPf.isospacesel3D ) IsoSpaceSel_3D();

  if ( VPf.labelsOn ) {
    if ( is.ballmode || is.stickmode ) xcMakeBallLabel3D();
    if ( is.spacefillmode ) xcMakeSpaceLabel3D();
  }    

  /* t.k: this is temporal:
   * VORONOI & WIGNER-SEITZ
   */
  if ( xcr.lprimwigner && VPf.wignerseitz ) xcRenderVoronoi();

  /*
   * SUPERCELL OPTION
   */
  if ( VPf.supercell) (*xcSuperCell)();

  /*
   * NEW SURFACE OPTIONS 
   */
  if ( VPf.nsurface > 0 ) xcRenderSurface();
  
  /* are we in selection mode ??? */
  if ( VPf.selection ) {
    xcRenderSelAtoms3D();
    xcRenderSelBonds3D();
    /*for(i=0; i<MAXSEL-1; i++) {*/
      /*if ( glIsList(SphereSelList[i]) ) glCallList(SphereSelList[i]);*/
      /*if ( glIsList(LineSelList[i]) ) glCallList(LineSelList[i]);*/
    /*}*/
    /*if ( glIsList(SphereSelList[MAXSEL-1]) ) glCallList(SphereSelList[3]);*/
  }

  /* are we in atomadd mode ??? */
  if ( VPf.atomadd ) {
    xcDisplayAddAtomBox();    
    
    /* this one is not displayed because it cause problem */
    /* this is the only bug i have seen */

    /* check if coordinate sistem should be displayed */
    /*  if ( VPf.xyzOn ) xcDisplayXYZ(); */
    
    
  }
  
  glPopMatrix();
}


/*
  Move the camera to the home position
  Or to a predefined stereo configuration
  The model is assumed to be in a 10x10x10 cube
  Centered at the origin
*/
static void 
CameraHome()
{
  camera.aperture = 60;

   
  camera.pr = origin;
  camera.vd.x = 1;
  camera.vd.y = 0;
  camera.vd.z = 0;

  camera.vu.x = 0;
  camera.vu.y = 1;
  camera.vu.z = 0;
  
  camera.vp.x = 0;
  camera.vp.y = 0;
  camera.vp.z = 0;

  /*for the following value :*/
  /*    5 -> the scene semm to be on the screen */
  /*   15 -> the scene semm to be out of the screen */
  /*   10 -> the scene semm to be in the screen */
   
   
  camera.focallength = 5;
   
   
  /* Non stressful stereo setting */
  camera.eyesep = camera.focallength / 30.0;
   
}

/*
  Normalise a vector
*/
static void 
Normalise(XYZ *p)
{
  double length;

  length = sqrt(p->x * p->x + p->y * p->y + p->z * p->z);
  if (length != 0) {
    p->x /= length;
    p->y /= length;
    p->z /= length;
  } else {
    p->x = 0;
    p->y = 0;
    p->z = 0;
  }
}

/*
  Calculate the unit normal at p given two other points
  p1,p2 on the surface. The normal points in the direction
  of p1 crossproduct p2
*/
/* static XYZ  */
/* CalcNormal(XYZ p,XYZ p1,XYZ p2) */
/* { */
/*   XYZ n,pa,pb; */

/*   pa.x = p1.x - p.x; */
/*   pa.y = p1.y - p.y; */
/*   pa.z = p1.z - p.z; */
/*   pb.x = p2.x - p.x; */
/*   pb.y = p2.y - p.y; */
/*   pb.z = p2.z - p.z; */
/*   Normalise(&pa); */
/*   Normalise(&pb); */

/*   n.x = pa.y * pb.z - pa.z * pb.y; */
/*   n.y = pa.z * pb.x - pa.x * pb.z; */
/*   n.z = pa.x * pb.y - pa.y * pb.x; */
/*   Normalise(&n); */

/*   return(n); */
/* } */




/********************************************************************
 ********************************************************************
                              EV
 ********************************************************************
 ********************************************************************/
