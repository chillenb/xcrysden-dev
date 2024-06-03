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
* Source: $XCRYSDEN_TOPDIR/C/xcSelect.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

extern void xcdebug(const char *text);
#include <togl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tk.h>
#include "struct.h"
#include "displayfunc.h"
#include "xcfunc.h"

#define APPR_SCREENSIZE 1024

#define PIXTOL 5  /* tolerance in screen pixels for selection */
#define SELTOL -1.0   /* tolerance for selection; this is just a flag,
		       * because tolerance is either radius (rball, atrad)
		       * or rodr
		       */
#define SELCOL 0  /* index of selected color */
#define SELLINE_WIDTH 2.0   /* width of "selline" */
#define SELLINE_PAT 0xAAAA  /* selline stipple pattern */
extern OrthoProj ort; 
extern Options3D is;
extern char *element[];
extern GLfloat DefAtCol[MAXNAT + 1][3];
extern GLuint PointList;
extern GLuint tempDisable3Dlist, tempEnable3Dlist; /* this is for switching 
						      on/off lighting */
extern RasterFontSize rf;
extern PERSPECTIVE    persp;
extern realTimeMove makeMovie;

/* this is where data of selected atoms go */
typedef struct {
  int id;     /* ID of selected atom */
  int sqn;    /* SQN of selected atom */
} SelAtom;
SelAtom asel[MAXSEL + 1];       /* asel is counted from 1..MAXSEL */

typedef struct {
  int nat;
  float xat;
  float yat;
} AtomsCoor;  
AtomsCoor *acrd; /* here we don't have zat, because of 'hpsort_index1' 
		  * function, which
		  * don't allow usage of struct-variables and also 
		  * for selections only xat & yat are important
		  */
double *acrd_z;   /* in above struct, we don't have "Z-orientation"
		     variable, this is one */

int old_nobjects; /* it tells how many objects were there before
		     selection process begins */
int *zindex; /* here will go an Z-orientation order of atoms */

float tol;  /* X and Y tolerance for selections */
int nsel = 0;   /* number of selected atoms */
int nline = 0;  /* id of selected line */
GLboolean isWF; /* WF and PL uses the same Sel. Disp Function, this is a flag
		   for distinguishing between this two modes */

/* --- this is for selection-lists --- */
/*GLuint SphereSelList[MAXSEL];*/
/*GLuint LineSelList[MAXSEL - 1];*/

/* --- temporary variables needed for debugging */
typedef struct {
  float xx;
  float yy;
} SELECT_DEBUG;
SELECT_DEBUG sel_debug;

/*****************************************************************************/
/* --- function prototypes --- */
int XC_SelectCmd(ClientData clientData, Tcl_Interp *interp, 
		 int argc, const char *argv[]);
static float dist2(float x, float y);
static void xcUpdateCoorf( double mat[][4], float *x, float *y, double *z);
void (*old_xcDisplay)(struct Togl *togl) = 0x0;
static void old_xcDisplayFunc( void (*Func)(struct Togl *togl) );
static void MapToSelDispFunc(void);
int xcSelectSqn(int x, int y);
static void xcSelectAtom(int sqn, char *atomdata);
/*
  static void xcMakeSelBallList(int sqn);
  static void xcMakeNewSelBallList( GLdouble *sizeArray , int natom );
*/
static void  xcSelectLine(int sqn1, int sqn2);
void GetSelCylinderPar(double x21, double y21, double z21, 
		       double *xrvb, double *yrvb, double *zrvb, 
		       double *fibond, double *bondl);
void xcPointLineSel(struct Togl *togl);
void xcBallStick1Sel(struct Togl *togl);
void xcBallStick2Sel(struct Togl *togl);
void xcRenderSelAtoms3D(void);
void xcRenderSelBonds3D(void);
int XC_DeselectCmd(ClientData clientData, Tcl_Interp *interp,
		   int argc, const char *argv[]);
static void PushOneCoor(int i, int ii);
/*static void xcUpdateSelLists(int id);*/
/*static void SelectDebug(struct Togl *togl);*/

/* ================================== */
/* --- extern function prototypes --- */
/* --- xcDisplayFunc2.c --- */
void xcUpdateCoor( double mat[][4], double *x, double *y, double *z);
/* --- xcDisplayFunc.c --- */
extern void xcPipeBall2D(struct Togl *togl);
extern void xcSpaceFill2D(struct Togl *togl);
extern void (*xcDisplay)(struct Togl *togl); 
extern void xcMakePointList(void);
extern void xcDisplayFunc( void (*Func)(struct Togl *togl) );
extern void xcSolidSphere(GLdouble radius);
extern void xcSolidBond(GLdouble radius, GLdouble height, int bondFlag);
extern GLuint findList1(int lindex, GLdouble *paramArray, int size);
GLuint makeModelPtr1(int lindex, GLdouble *sizeArray, int count);
/* extern void xcMakeNewBallList( GLdouble *sizeArray , int natom ); */
extern void xcDisplayXYZ(void);
/* --- auxilproc.c ---*/
extern double dist6(double x1, double x2, 
		    double y1, double y2, double z1, double z2);
extern double dist3(double x, double y, double z);
extern int MakeSticks2(int i, int col, int flag);

/* --- xcviewport.c --- */
extern void Screen2Model_Coords(int xs, int ys, float *xm, float *ym);


/* XC_SelectCmd --> inplementation of 'xc_select' custom Tcl command
 * ---------------------
 * Usage: xc_select <toglName> begin|sqn|atom|line|clean|flush \
 *                             ?value1? ?value2? .....
 *		
 *              xc_select begin -- before we begin selection we must prepare 
 *                                 and do some initialisations
 *
 *              xc_select sqn %x %y -- a sequential number of selected atom
 *                                     is returned (this number is an ID)
 *              
 *              xc_select atom <sqn> - <sqn> atom is painted to selection color
 *
 *              xc_select line <sqn1> <sqn2> - a line between atom <sqn1> and
 *                                             atom <sqn2> is made
 *
 *              xc_select clean  -- when we are done with selection -> 
 *                                  to clean up
 *
 *              xc_select flush  -- flush(display) the selection changes
 *
 *              xc_select finish -- finish(display) the selection changes
 *                                  for 2D is the same as xc_select flush,
 *                                  but for 3D it does (*xcDisplay)(), whereas
 *                                  xc_select flush just swap the buffers
 */	
int 
XC_SelectCmd(ClientData clientData, Tcl_Interp *interp,
	     int argc, const char *argv[])
{
  struct Togl *togl;
  int x,y; /* mesaWin coord of Button-1 */

  if (argc < 3) {
    Tcl_SetResult(interp, "Usage: xc_select <toglName> sqn|atom|line|clean|flush|finish ?value1? ?value2?", TCL_STATIC);
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

  /*
   * xc_select BEGIN 
   */
  if ( strcmp(argv[2],"begin") == 0 ) {
    /* now do some initializations !!!! */
    old_nobjects = tmp_nobjects;     /* remember number of tmp_nobjects */
    old_xcDisplayFunc(xcDisplay);    /* remember display function that was on
				      * before selection proces begin
				      */   
    /* flag to know when we are in "selection mode" */
    VPf.selection = GL_TRUE;

    /* now map from xcDisplay to function that has properties of xcDisplay +
       stuff for selection */
    MapToSelDispFunc();
    /* now malloc acrd structure and zindex */
    acrd = (AtomsCoor *) malloc(sizeof(AtomsCoor) * (natoms + 1));
    acrd_z = (double *) malloc(sizeof(double) * (natoms + 1));
    zindex = (int *) malloc(sizeof(int) * (natoms + 1));
    /* set nsel to 0 */
    nsel = 0;

    /******************************************************************/
    /* there is a strange efect due to reconfiguring of Mesa window. 
       Because we want to change cursor of Mesa window for selection mode, 
       this cause display of Mesa window to vanish. To fix that just call 
       tkSwapBuffers(); */
    /* tkSwapBuffers(); */
  }

  /*
   * xc_select SQN %x %y
   *
   * %x %y are toglName Button-1 coord. 
   */
  else if ( strcmp(argv[2],"sqn") == 0 ) {
    int sqn;
    char *csqn = Tcl_Alloc( sizeof(char) * 10 );
    /* argc must be 5 */
    if ( argc != 5 ) {
      Tcl_SetResult(interp, "wrong usage of \"xc_select <toglName> sqn\" command, must be \"xc_select <toglName> sqn <x> <y>\"", TCL_STATIC);
      return TCL_ERROR;
    }
    /* agrv[3],argv[4] must be integers */
    if ( Tcl_GetInt(interp,argv[3],&x) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"got invalid value for <x>, while executing %s %s %s %s %s\n",argv[0],argv[1],argv[2],argv[3],argv[4]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( Tcl_GetInt(interp,argv[4],&y) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"got invalid value for <y>, while executing %s %s %s %s %s\n",argv[0],argv[1],argv[2],argv[3],argv[4]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }  
    /* if there were no "xc_selct begin" than ERROR */
    if ( !VPf.selection ) {
      Tcl_SetResult(interp, "calling \"xc_select sqn\" command, before initialising selection process with \"xc_select begin\"", TCL_STATIC);
      return TCL_ERROR;
    }
    /* now find sequential number of atom, if no atom find, retun 0 */
    sqn = xcSelectSqn(x,y);
    /* BEGIN -- debug only */
    /* SelectDebug( togl ); */
    /* END   -- debug only */
    sprintf(csqn,"%d",sqn);
    Tcl_SetResult(interp,csqn,TCL_DYNAMIC);
  }

  /*
   * xc_select ATOM <sqn>
   */
  else if ( strcmp(argv[2],"atom") == 0 ) {
    int sqn;
    char *atomdata = Tcl_Alloc(sizeof(char) * 80 );
    if ( argc != 4 ) {
      Tcl_SetResult(interp, "wrong usage of \"xc_select <toglName> atom\" command, must be \"xc_select <toglName> atom <sqn>\"", TCL_STATIC);
      return TCL_ERROR;
    }
    /* agrv[3] must be integer */
    if ( Tcl_GetInt(interp,argv[3],&sqn) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"got invalid value for <sqn>, while executing %s %s %s %s\n",argv[0],argv[1],argv[2],argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    /* sqn must be in [1,natoms], else nothing good */
    if ( sqn < 1 || sqn > natoms ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"got invalid value for <sqn>, while executing %s %s %s %s; must be from 1 to %d",argv[0],argv[1],argv[2],argv[3],natoms);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    /* if there were no "xc_selct begin" than ERROR */
    if ( !VPf.selection ) {
      Tcl_SetResult(interp, "calling \"xc_select atom\" command, before initialising selection process with \"xc_select begin\"", TCL_STATIC);
      return TCL_ERROR;
    }

    xcSelectAtom(sqn, atomdata);
    Tcl_SetResult(interp,atomdata,TCL_DYNAMIC);
  }

  /***************************************************************************/
  /* implementation of "XC_SELECT <toglName> LINE <sqn1> <sqn2> */
  else if ( strcmp(argv[2],"line") == 0 ) {
    int sqn1, sqn2;
    if ( argc != 5 ) {
      Tcl_SetResult(interp, "wrong usage of \"xc_select <toglName> line\" command, must be \"xc_select <toglName> line <sqn1> <sqn2>\"", TCL_STATIC);
      return TCL_ERROR;
    }
    /* agrv[3],argv[4] must be integers */
    if ( Tcl_GetInt(interp,argv[3],&sqn1) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"got invalid value for <sqn1>, while executing %s %s %s %s %s\n",argv[0],argv[1],argv[2],argv[3],argv[4]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( Tcl_GetInt(interp,argv[4],&sqn2) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"got invalid value for <sqn2>, while executing %s %s %s %s %s\n",argv[0],argv[1],argv[2],argv[3],argv[4]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    /* sqn1,sqn2 must be in [1,natoms], else nothing good */
    if ( (sqn1 < 1 && sqn1 > natoms) || (sqn2 < 1 && sqn2 > natoms) ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"got invalid value for <sqn1> or <sqn2>, while executing %s %s %s %s %s; must be from 1 to %d",argv[0],argv[1],argv[2],argv[3],argv[4],natoms);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( sqn1 == sqn2 ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"<sqn1> and <sqn2> are the same; while executing %s %s %s %s %s",argv[0],argv[1],argv[2],argv[3],argv[4]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    /* when calling XC_SELECT LINE, at least two atoms must be selected !!!! */
    if ( nsel < 2 ) {
      Tcl_SetResult(interp, "at least two atoms must be selected, before calling \"xc_select line <sqn1> <sqn2>\"", TCL_STATIC);
      return TCL_ERROR;
    }
    /* if there were no "xc_selct begin" than ERROR */
    if ( !VPf.selection ) {
      Tcl_SetResult(interp, "calling \"xc_select line\" command, before initialising selection process with \"xc_select begin\"", TCL_STATIC);
      return TCL_ERROR;
    }

    xcSelectLine(sqn1, sqn2);
  }

  /***************************************************************************/
  /* implementaion of "XC_SELECT <toglName> CLEAN"; to clean */
  else if ( strcmp(argv[2],"clean") == 0 ) {  
    /* if there were no "xc_select <toglName> begin" than ERROR */
    if ( !VPf.selection ) {
      Tcl_SetResult(interp, "calling \"xc_select clean\" command, before initialising selection process with \"xc_select begin\"", TCL_STATIC);
      return TCL_ERROR;
    }    

    /* cleaning up */
    tmp_nobjects = old_nobjects;   
    xcDisplayFunc(old_xcDisplay);  /* going to DisplayFunc that was used before
				      selections */

    /* delete all lists */
    /*glDeleteLists( SphereSelList[0], MAXSEL );*/
    /*glDeleteLists( LineSelList[0], MAXSEL - 1);*/

    VPf.selection = GL_FALSE;
    /* now free acrd structure, zindex & acrd_z */
    free(acrd);
    free(acrd_z);
    free(zindex);
    /* now Z-reorientate */
    if ( dimType == XC_2D) hpsort_index1(tmp_nobjects, zorient, iwksp);
    /* now update display */
    Togl_PostRedisplay(togl);
  }

  /***************************************************************************/
  /* XC_SELECT FLUSH */
  else if ( strcmp(argv[2],"flush") == 0 ) {
    if ( dimType == XC_2D ) Togl_PostRedisplay(togl); 
    if ( dimType == XC_3D ) {
      glFlush();
      Togl_SwapBuffers(togl);
    }
  }

  /***************************************************************************/
  /* XC_SELECT FINISH */
  else if ( strcmp(argv[2],"finish") == 0 ) Togl_PostRedisplay(togl); 

  /***************************************************************************/
  /* unknown option */
  else {
    char rss[1024];
    snprintf(rss, sizeof(rss),"unknown option \"%s\" for \"xc_select\" command, must be one of begin, sqn, atom, line, clean, flush, finish",argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  return TCL_OK;
}


/*****************************************************************************/
static float
dist2(float x, float y)
{
  /*return( sqrtf(x * x + y * y) );*/
  return( sqrt(x * x + y * y) );
}

static void
xcUpdateCoorf( double mat[][4], float *x, float *y, double *z)
{
  float x1, y1;
  double z1;

  x1 = mat[0][0] * *x + mat[1][0] * *y + mat[2][0] * *z;
  y1 = mat[0][1] * *x + mat[1][1] * *y + mat[2][1] * *z;
  z1 = mat[0][2] * *x + mat[1][2] * *y + mat[2][2] * *z;
  *x = x1;
  *y = y1;
  *z = z1;
}


/*****************************************************************************/
static void
old_xcDisplayFunc( void (*Func)(struct Togl *togl) )
{
  old_xcDisplay = Func;
}


/*****************************************************************************/
static void
MapToSelDispFunc(void)
{
  /* now find out what is the current display function, and according to that
   * assign new display function, that will have all necessary stuff for 
   * selections
   *
   * now we must also add some tolerance to xx,yy; that means if we are close
   * enough to atom, atom will be selected;
   * if 2D --> if WF or PL tolerance will be PIXTOL pixels
   *
   * if 3D --> "sticks" tolerance is rod radius
   *
   * else tolerance is either rball or atrad
   */

  if (dimType == XC_2D) {
    /* WF and PL display functions will have common selection function, because
     * selection atom & line are visual the same for this two modes
     */
    if ( xcDisplay == xcWireFrame2D ||
	 xcDisplay == xcWireFrameL2D ||
	 xcDisplay == xcWireFrameF2D ||
	 xcDisplay == xcWireFrameFL2D ) {
      /* make PointList */
      xcMakePointList();
      isWF = 1;
      xcDisplayFunc(xcPointLineSel);    
      tol = (float) PIXTOL / VPf.scrAnX;
    }
    else if ( xcDisplay == xcPointLine2D ||
	      xcDisplay == xcPointLineL2D ||
	      xcDisplay == xcPointLineF2D ||
	      xcDisplay == xcPointLineFL2D ) {
      isWF = 0;
      xcDisplayFunc(xcPointLineSel);    
      tol = (float) PIXTOL / VPf.scrAnX;      
    }
    else if ( xcDisplay == xcBallStick12D ||
	      xcDisplay == xcBallStick1L2D ||
	      xcDisplay == xcBallStick1F2D ||
	      xcDisplay == xcBallStick1FL2D ) {
      xcDisplayFunc(xcBallStick1Sel);
      tol = SELTOL;      
    }
    else if ( xcDisplay == xcBallStick22D ||
	      xcDisplay == xcBallStick2L2D ||
	      xcDisplay == xcBallStick2F2D ||
	      xcDisplay == xcBallStick2FL2D ) {
      xcDisplayFunc(xcBallStick2Sel);
      tol = SELTOL;
    }
    else if ( xcDisplay == xcPipeBall2D ) {
      xcDisplayFunc(xcPipeBall2D);
      tol = SELTOL;
    }
    else if ( xcDisplay == xcSpaceFill2D ) {
      xcDisplayFunc(xcSpaceFill2D);
      tol = SELTOL;
    }
  } 
  else if ( dimType == XC_3D ) { 
    if (is.stickmode && !is.ballmode) 
      tol = 1.0 * rrod;
    else tol = SELTOL;      
  }
}



/*****************************************************************************/
/* this function find an atom on (x,y) position, if there is one; 
 * always find the topmost one;
 * x,y are mesaWin Button-1 coord.
 */
int
xcSelectSqn(int x, int y)
{
  float xx,yy;  /* coordinates of Button-1 in Angstroms */
  float dis, tol_;
  int i;

  /* transform screen coordinates to model coordinates */
  Screen2Model_Coords( x, y, &xx, &yy);
  /* TEST--TEST */
  sel_debug.xx = xx;
  sel_debug.yy = yy;
  printf("xcSelectSqn> x=%d, y=%d;   xx=%f, yy=%f\n",x, y, xx,yy);
  fflush(stdout);
  /* TEST--END */
  /* now write acrd struct and rotate according to vec.crdmajor */
  for(i=1; i<=natoms; i++)
    {
      (acrd + i)->nat = *(nat + i);
      (acrd + i)->xat = (float) *(xat + i);
      (acrd + i)->yat = (float) *(yat + i);
      *(acrd_z + i) = *(zat + i);	
      xcUpdateCoorf( vec.crdmajor, &(acrd + i)->xat, &(acrd + i)->yat, 
		     (acrd_z + i) );      
    }
  /* now orientate according to Z */
  hpsort_index1(natoms, acrd_z, zindex);

  /* now try to find an atom, if there is no atom, just return 0 */
  for(i=natoms; i>=1; i--)
    {

      tol_ = tol;

      /* now take care of tolerance */
      if ( tol == SELTOL ) {
	/* first take care of XC_2D modes */
	if ( dimType == XC_2D ) {
	  /* check if balls are >= PIXTOL */
	  if ( (rball[ nat[zindex[i]] ] * VPf.scrAnX) >= PIXTOL ) 
	    tol_ = rball[ nat[zindex[i]] ];
	  else tol_ = (float) PIXTOL / VPf.scrAnX;
	}
	else if ( dimType == XC_3D ) {
	  /* are we in SF mode */
	  if (is.spacefillmode) tol_ = atrad[ nat[zindex[i]] ];
	  else if (is.ballmode) {
	    if ( (rball[ nat[zindex[i]] ] * VPf.scrAnX) >= PIXTOL ) 
	      tol_ = rball[ nat[zindex[i]] ];
	    else tol_ = (float) PIXTOL / VPf.scrAnX;
	  }

	}
      }

      dis = dist2( xx - (acrd + zindex[i])->xat, 
		   yy - (acrd + zindex[i])->yat );
      if ( dis <= tol_ ) {
	fprintf(stderr, "DEBUG> tol_ = %f, dis = %f; VPf.scrAnX = %f\n", tol_, dis, VPf.scrAnX);
	/* WE FIND AN ATOM */
	return zindex[i];
      }
    }

  /* if we come here we didn't find an atom */
  return 0;
}


/*****************************************************************************/
/* this function select an atom and return it's data
 */
static void xcSelectAtom(int sqn, char *atomdata)
{  
  int i;
  /*char *atomdata = (char *) malloc(sizeof(char) * 80);*/

  nsel++; /* number of selected atoms is increased */

  /* put information of selected atom to 'asel' structure */
  if ( dimType == XC_2D ) (asel + nsel)->id = tmp_nobjects + 1;
  (asel + nsel)->sqn = sqn;

  /* if we are in 2D make "selatom" and redraw, if in 3D just draw "selatom" */
  if ( dimType == XC_2D ) {
    i = ++tmp_nobjects; /* increase TMP_NOBJECTS for one */    
    /* add "selatom" to coor structure */
    (coor + i)->flag = SELATOM;
    (coor + i)->sqn  = sqn;
    (coor + i)->nat  = *(nat + sqn);
    (coor + i)->x1   = *(xat + sqn);
    (coor + i)->y1   = *(yat + sqn);
    (coor + i)->z1   = *(zat + sqn);
    xcUpdateCoor( vec.crdmajor, 
		  &(coor + i)->x1, &(coor + i)->y1, &(coor + i)->z1 );

    /* add a little to Z-coor, because SELATOM  must go over the ATOM */
    *(zorient + i)  = (coor + i)->z1 + Z_OFFSET((coor + i)->flag);
    hpsort_index1(tmp_nobjects, zorient, iwksp);

    /* if BS1 or BS2 mode */
    /*if ( xcDisplay != xcPointLineSel ) {            
      xcMakeSelBallList(i);
      }*/

    /* now update a display */
    /* Togl_PostRedisplay(togl); I THINK THIS IS NOT NEEDED */
  }
  else if ( dimType == XC_3D ) {
    ;
  }    

  /* now return atomic data */
  sprintf(atomdata," %-3d %-3s %-3d  %+14.9f %+14.9f %+14.9f", sqn, 
	  element[*(nat + sqn)], *(nat + sqn), *(xat + sqn) + mx, 
	  *(yat + sqn) + my, *(zat + sqn) + mz);
}


/*
  static void xcMakeSelBallList(int sqn) {
  GLint     nstep, i, nat; 
  GLdouble  sine, cosine;
  GLdouble  sizeArray[3];
  GLuint    displayList[2];

  nat   = (coor + sqn)->nat; 
  nstep = 3 * CalcTessFactor();
  if ( nstep < 8 )  nstep = 8; 

  sizeArray[0]   = rball[nat];
  sizeArray[1]   = (GLdouble) nstep; 
  sizeArray[2]   = (GLdouble) VPf.OUTlinewidth;
  displayList[0] = findList1 (BALL, sizeArray, 3);
  displayList[1] = findList1 (OUTLINEBALL, sizeArray, 3);

  if ( displayList[0] == 0 ) {  
  glNewList( makeModelPtr1 (BALL, sizeArray, 3),
  GL_COMPILE_AND_EXECUTE ); 
  glBegin(GL_POLYGON);
  for(i=0; i < nstep; i++) {
  cosine = rball[nat] * cos((GLdouble) i * 2.0 * PI / (GLdouble) nstep);
  sine   = rball[nat] * sin((GLdouble) i * 2.0 * PI / (GLdouble) nstep);
  glVertex2d( cosine, sine );
  }
  glEnd();
  glEndList();

  glNewList( makeModelPtr1 (OUTLINEBALL, sizeArray, 3),
  GL_COMPILE_AND_EXECUTE );			      
  glColor3f (0.0, 0.0, 0.0); 
  glBegin(GL_LINE_LOOP); 
  for(i=0; i < nstep; i++) { 
  cosine = rball[nat] * cos((GLdouble) i * 2.0 * PI / (GLdouble) nstep); 
  sine   = rball[nat] * sin((GLdouble) i * 2.0 * PI / (GLdouble) nstep); 
  glVertex2d( cosine, sine ); 
  } 
  glEnd(); 
  glEndList();
  } else {
  glCallList (displayList[0]);
  glCallList (displayList[1]);
  }
  }
*/
/* / */

/* static void */
/* xcMakeSelBallList(int sqn) */
/* { */
/*   GLuint displayList; */
/*   GLdouble *sizeArray, *tmp; */

/*   sizeArray = (GLdouble *) malloc (sizeof (GLdouble) * 4);  */
/*   tmp = sizeArray; */
/*   *tmp++ = rball[(coor + sqn)->nat]; */

/*   displayList = findList1 (BALL, sizeArray, 1); */
/*   if (displayList == 0) xcMakeNewSelBallList( sizeArray, sqn ); */
/*   else { */
/*     (coor + sqn)->list2 = displayList; */
/*     free (sizeArray); */
/*   } */
/* } */


/**/
/* static void */
/* xcMakeNewSelBallList( GLdouble *sizeArray , int natom ) */
/* { */
/*   GLint nstep, i; */
/*   int nab = (coor + natom)->nat; */
/*   GLdouble sine, cosine, dstep; */

/*   nstep = (GLint) (3 * APPR_SCREENSIZE / (VPf.tessFactor * ort.size)); */

/*   if ( nstep < 8 ) nstep = 8; */
/*   dstep = (GLdouble) nstep; */
/*   (coor + natom)->list1 = makeModelPtr1 (BALL, sizeArray, 4); */
/*   glNewList( (coor + natom)->list2, GL_COMPILE ); */
/*     glBegin(GL_POLYGON); */
/*       for(i=0; i < nstep; i++) */
/*         { */
/* 	  cosine = rball[nab] * cos((GLdouble) i * 2.0 * PI / dstep); */
/* 	  sine = rball[nab] * sin((GLdouble) i * 2.0 * PI / dstep); */
/* 	  glVertex2d( cosine, sine ); */
/* 	} */
/*     glEnd(); */
/*   glEndList(); */
/* } */
/**/ 

static void  
xcSelectLine(int sqn1, int sqn2)
{
  double dis, len, fdis;
  int i;

  printf("\n\nsqn1=%d, sqn2=%d\n",sqn1,sqn2);
  /* now, if between sqn1 & sqn2 is chemmical bond, display line as a bond,
     else display as a dotted line */
  dis = dist6(*(xat + sqn1), *(xat + sqn2),
	      *(yat + sqn1), *(yat + sqn2),
	      *(zat + sqn1), *(zat + sqn2));
  len = rcov[*(nat + sqn1)] + rcov[*(nat + sqn2)];

  if ( dimType == XC_2D ) {
    /* tmp_nobjects must be increased for one, because we have new object */
    i = ++tmp_nobjects;  
    printf("lineID=%d\n",i);
    if ( dis < len ) {
      fdis = rcov[*(nat + sqn1)] / len;
      /* bond to ATOM sqn1 */
      (coor + i)->flag = SELBOND;
      (coor + i)->nat = nat[sqn1];
      (coor + i)->bondend = BOND_ATOM_TO_MIDBOND;
      (coor + i)->x1 = xat[sqn1];
      (coor + i)->y1 = yat[sqn1];
      (coor + i)->z1 = zat[sqn1];
      (coor + i)->x2 = xat[sqn1] + fdis * (xat[sqn2] - xat[sqn1]);
      (coor + i)->y2 = yat[sqn1] + fdis * (yat[sqn2] - yat[sqn1]);
      (coor + i)->z2 = zat[sqn1] + fdis * (zat[sqn2] - zat[sqn1]);
      xcUpdateCoor( vec.crdmajor, &(coor + i)->x1, &(coor + i)->y1, 
		    &(coor + i)->z1 );
      xcUpdateCoor( vec.crdmajor, &(coor + i)->x2, &(coor + i)->y2, 
		    &(coor + i)->z2 );
      *(zorient + i)  = 0.5 * ((coor + i)->z1 + (coor + i)->z2) + 
	Z_OFFSET((coor + i)->flag);

      /* bond to ATOM sqn2 */
      i = ++tmp_nobjects;
      printf("lineID=%d\n",i);
      (coor + i)->flag = SELBOND;
      (coor + i)->nat = nat[sqn2];
      (coor + i)->bondend = BOND_MIDBOND_TO_ATOM;
      (coor + i)->x2 = xat[sqn2];
      (coor + i)->y2 = yat[sqn2];
      (coor + i)->z2 = zat[sqn2];
      (coor + i)->x1 = xat[sqn2] + (1 - fdis) * (xat[sqn1] - xat[sqn2]);
      (coor + i)->y1 = yat[sqn2] + (1 - fdis) * (yat[sqn1] - yat[sqn2]);
      (coor + i)->z1 = zat[sqn2] + (1 - fdis) * (zat[sqn1] - zat[sqn2]);      
    }
    else {
      (coor + i)->flag = SELLINE;
      (coor + i)->nat = 0;  /* 0 is for SELCOL */
      (coor + i)->x1 = xat[sqn1];
      (coor + i)->y1 = yat[sqn1];
      (coor + i)->z1 = zat[sqn1];
      (coor + i)->x2 = xat[sqn2];
      (coor + i)->y2 = yat[sqn2];
      (coor + i)->z2 = zat[sqn2];
    }

    xcUpdateCoor( vec.crdmajor, &(coor + i)->x1, &(coor + i)->y1, 
		  &(coor + i)->z1 );
    xcUpdateCoor( vec.crdmajor, &(coor + i)->x2, &(coor + i)->y2, 
		  &(coor + i)->z2 );

    /* if ATOM and BOND and FRAME and SEL BOND are on the same plane, 
     * ATOM must be drawn first, than BOND, than FRAME, and than SELBOND,
     * so we will add very little to z-coord, but more than for FRAME
     */
    *(zorient + i)  = 0.5 * ((coor + i)->z1 + (coor + i)->z2) + Z_OFFSET((coor + i)->flag); 
    hpsort_index1( tmp_nobjects, zorient, iwksp);
    /* now update a display */
    /* (*xcDisplay)(); */
  }
  else if ( dimType == XC_3D ) {
    ;
  }    
}


void 
GetSelCylinderPar(double x21, double y21, double z21, 
		  double *xrvb, double *yrvb, double *zrvb, 
		  double *fibond, double *bondl)
{  
  double mintol = -1.0 + 1.0e-8;

  /* rotated vector is cross product between (0,0,1) & vector submited 
     to GetCylinderPar */
  *bondl = dist3(x21, y21, z21);  

  /* normalize vectror (x21,y21,z21) */
  normalizepv(&x21, &y21, &z21);
  if (z21 > mintol) {
    *xrvb = -y21;
    *yrvb = x21;
    *zrvb = 0.0;
  } else {
    *xrvb = 1.0;
    *yrvb = 0.0;
    *zrvb = 0.0;
  }

  if (ABS(z21) > 1.0000000001) {
    fprintf(stderr,"++ Z-component of normalized vector greater than 1.0\n");
    fprintf(stderr,"++ Vector: %f %f %f\n", x21, y21, z21);
    return;
  }
  *fibond = RAD2DEG * acos(z21);
}


/*****************************************************************************/
/*                      Selection Display Functions                          */
/*****************************************************************************/

/* display function for WF & PL 2D's modes for selections */
void
xcPointLineSel(struct Togl *togl)
{
  int i;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  for(i = 1; i <= tmp_nobjects; i++) { 
    if ( (coor + *(iwksp + i))->flag == BOND )
      {
	glLineWidth( (GLfloat) VPf.PLlinewidth );  
	glBegin(GL_LINES);
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
	/* display point only if old_xcDisplay is one of PL disp. functions */
	if ( isWF == 0 ) {
	  /* first display point */
	  glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	  glPushMatrix();
	  glTranslated( (coor + *(iwksp + i))->x1, 
			(coor + *(iwksp + i))->y1, 
			(coor + *(iwksp + i))->z1 );
	  glCallList( PointList );	
	  glPopMatrix();	
	}

	/* display labels only if VPf.lablesOn  */
	if ( VPf.labelsOn ) {
	  makeAtomLabel2D ( *(iwksp + i), rf.wp2, rf.hp2 );
	}
      }
    /* display frames only if VPf.framessOn */
    if ( VPf.framesOn ) {
      if ( (coor + *(iwksp + i))->flag == FRAME )
	{
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

    if ( (coor + *(iwksp + i))->flag == SELATOM )
      {
	glColor3fv( atcol[SELCOL] );
	glPushMatrix();
	glTranslated( (coor + *(iwksp + i))->x1, 
		      (coor + *(iwksp + i))->y1, 
		      (coor + *(iwksp + i))->z1 );
	glCallList( PointList );	
	glPopMatrix();
      }

    if ( (coor + *(iwksp + i))->flag == SELLINE ) 
      {
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
    if (  (coor + *(iwksp + i))->flag == SELBOND ) {
      {
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
      }
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


void
xcBallStick1Sel(struct Togl *togl)
{
  int i;
  GLdouble x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  glLineWidth( VPf.OUTlinewidth );
  glEnable( GL_CULL_FACE );
  glCullFace( GL_BACK );

  for(i = 1; i <= tmp_nobjects; i++) 
    {
      if ( (coor + *(iwksp + i))->flag == BOND )
	{
	  MakeSticks1(*(iwksp + i), 
		      &x1, &y1, &z1,   &x2, &y2, &z2, 
		      &x3, &y3, &z3,   &x4, &y4, &z4);
	  glColor3f( 0.0, 0.0, 0.0 );	    

	  if ( VPf.unibond) glColor3fv(unibondCol);
	  else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	  glBegin(GL_POLYGON);
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
      else if ( (coor + *(iwksp + i))->flag == ATOM ) 
	{
	  /* display ball */
	  glPushMatrix();
	  glTranslated( (coor + *(iwksp + i))->x1,
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	  glColor3fv (atm.col[(coor + *(iwksp + i))->sqn]);
	  xcBall2D (*(iwksp + i));
	  glPopMatrix();

	  if (VPf.labelsOn) {
	    makeAtomLabel2D ( *(iwksp + i), CENTERED_LABEL, CENTERED_LABEL);
	  }            
	}
      else if ( (coor + *(iwksp + i))->flag == SELATOM ) 
	{
	  /* display selected ball */
	  glPushMatrix();
	  glTranslated( (coor + *(iwksp + i))->x1,
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	  glColor3fv( atcol[SELCOL] );
	  xcBall2D (*(iwksp + i));
	  /*glCallList( (coor + *(iwksp + i))->list2 );*/
	  glPopMatrix();

	  if (VPf.labelsOn) {
	    makeAtomLabel2D ( *(iwksp + i), CENTERED_LABEL, CENTERED_LABEL);
	  }            
	}

      /* display frames only if VPf.framesOn */
      if ( VPf.framesOn ) {
	if ( (coor + *(iwksp + i))->flag == FRAME )
	  {
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

      if ( (coor + *(iwksp + i))->flag == SELBOND )
	{
	  MakeSticks1(*(iwksp + i), 
		      &x1, &y1, &z1,   &x2, &y2, &z2, 
		      &x3, &y3, &z3,   &x4, &y4, &z4);
	  glColor3fv( atcol[SELCOL] );
	  glBegin(GL_POLYGON);
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

      if ( (coor + *(iwksp + i))->flag == SELLINE ) 
	{
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


void
xcBallStick2Sel(struct Togl *togl)
{
  int i;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  glLineWidth( VPf.OUTlinewidth );
  glEnable( GL_CULL_FACE );
  glCullFace( GL_BACK );

  for(i = 1; i <= tmp_nobjects; i++) 
    {
      if ( (coor + *(iwksp + i))->flag == BOND )
	{	  
	  glLineWidth( VPf.OUTlinewidth );  
	  /* MakeSticks2(*(iwksp + i), (coor + *(iwksp + i))->nat); */
	  MakeSticks2(*(iwksp + i), (coor + *(iwksp + i))->sqn, BALL);
	}
      else if ( (coor + *(iwksp + i))->flag == ATOM )
	{
	  /* display ball */
	  glPushMatrix();
	  glTranslated( (coor + *(iwksp + i))->x1,
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	  glColor3fv (atm.col[(coor + *(iwksp + i))->sqn]);
	  xcBall2D (*(iwksp + i));
	  glPopMatrix();

	  if (VPf.labelsOn) {
	    makeAtomLabel2D ( *(iwksp + i), CENTERED_LABEL, CENTERED_LABEL);
	  } 
	}
      else if ( (coor + *(iwksp + i))->flag == SELATOM ) 
	{
	  glPushMatrix();
	  glTranslated( (coor + *(iwksp + i))->x1,
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	  glColor3fv( atcol[SELCOL] );
	  xcBall2D (*(iwksp + i));
	  /*glCallList( (coor + *(iwksp + i))->list2 );*/
	  glPopMatrix();
	  if (VPf.labelsOn) {
	    makeAtomLabel2D ( *(iwksp + i), CENTERED_LABEL, CENTERED_LABEL);
	  } 
	}

      /* display frames only if VPf.framesOn */
      if ( VPf.framesOn ) {
	if ( (coor + *(iwksp + i))->flag == FRAME )
	  {
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

      if ( (coor + *(iwksp + i))->flag == SELBOND )
	{
	  glLineWidth( VPf.OUTlinewidth );
	  MakeSticks2(*(iwksp + i), 0, BALL);
	}

      if ( (coor + *(iwksp + i))->flag == SELLINE ) 
	{
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


void
xcRenderSelAtoms3D(void) {
  register int i, ind;
  register double r;

  glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, atcol[SELCOL] );

  for (i=1; i<=nsel; i++) {
    ind = asel[i].sqn;
    if(is.spacefillmode) r = atrad[*(nat + ind)];
    else if(is.ballmode) r = rball[*(nat + ind)];
    else r =  1.5 * rrod;
    glPushMatrix();
    glTranslated( xat[ind], yat[ind], zat[ind] );
    xcSolidSphere( 1.009 * r );
    glPopMatrix(); 
  }
}


void
xcRenderSelBonds3D(void) {
  register int    i, i1, i2;
  register double dis, len;
  double xrvb, yrvb, zrvb, fibond, bondl;

  glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, 
		atcol[SELCOL] );

  for (i=1; i<nsel; i++) {
    i1 = asel[i].sqn;
    i2 = asel[i+1].sqn;

    dis = dist6( xat[i1], xat[i2],  yat[i1], yat[i2],  zat[i1], zat[i2]);
    len = rcov[ nat[i1] ] + rcov[ nat[i2] ];

    if (dis < len) {
      GetSelCylinderPar(xat[i2]-xat[i1], yat[i2]-yat[i1], zat[i2]-zat[i1],
			&xrvb, &yrvb, &zrvb, &fibond, &bondl);
      glPushMatrix();
      glTranslated( *(xat + i1), *(yat + i1), *(zat + i1) );
      glRotated( fibond, xrvb, yrvb, zrvb );
      xcSolidBond( 1.2 * rrod, bondl, BOND_ATOM_TO_ATOM );
      glPopMatrix();
    } else {
      glCallList(tempDisable3Dlist);
      glEnable(GL_LINE_STIPPLE);
      glColor3fv( framecol );
      glLineWidth( SELLINE_WIDTH );
      glBegin(GL_LINES);
      glVertex3d( *(xat + i1), *(yat + i1), *(zat + i1) );
      glVertex3d( *(xat + i2), *(yat + i2), *(zat + i2) );
      glEnd();
      glDisable(GL_LINE_STIPPLE);   
      glCallList(tempEnable3Dlist);
    }
  }
}


/*****************************************************************************/
/* XC_DeselectCmd --> inplementation of 'xc_deselect' custom Tcl command
 * ---------------------
 * Usage: xc_deselect <toglName> atom|line ?sqn? 
 *		
 *              xc_deselect atom -- deselect atom <sqn>
 *
 *              xc_deselect line -- deselect all lines
 */	
int 
XC_DeselectCmd(ClientData clientData, Tcl_Interp *interp,
	       int argc, const char *argv[])
{
  struct Togl *togl;
  int sqn; /* mesaWin coord of Button-1 */
  int i, id, ii, j, jj;

  if (argc < 3) {
    Tcl_SetResult(interp, "Usage: xc_deselect <toglName> atom|line ?sqn?", TCL_STATIC);
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

  /***************************************************************************/
  /* first "XC_DESELECT <toglName> ATOM <SQN>" 
   * this command deselect <sqn>th atom !!! 
   */
  if ( strcmp(argv[2],"atom") == 0 ) {
    if (argc != 4) {
      Tcl_SetResult(interp, "Usage: xc_deselect <toglName> atom <sqn>", TCL_STATIC);
      return TCL_ERROR;
    }
    /* agrv[3] must be integer */
    if ( Tcl_GetInt(interp,argv[3],&sqn) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"got invalid value for <sqn>, while executing %s %s %s %s\n",argv[0],argv[1],argv[2],argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    /* first find an atom for deselection */
    id = -1;
    for(i=1; i<=nsel; i++) {
      printf("id=%d, sqn=%d\n",(asel + i)->id,(asel + i)->sqn);
      if ( (asel + i)->sqn == sqn ) {
	printf("DESELECTED:: ID=%d, SQN=%d\n",(asel + i)->id,(asel + i)->sqn); 
	if ( dimType == XC_2D ) id = (asel + i)->id;
	if ( dimType == XC_3D ) {
	  id = i;
	  printf("3D id=%d\n",id);
	  fflush(stdout);
	}
	/* atom (asel + i), will be deleted, that's (asel + i + 1) will 
	   became (asel + i) under condition that  (asel + i) is not the
	   last one */
	for(j=i; j<nsel; j++) {
	  if ( dimType == XC_2D ) (asel + j)->id = (asel + j + 1)->id - 1;
	  (asel + j)->sqn = (asel + j + 1)->sqn;
	}
	break;
      }
    }    
    /* if id is still -1, than we try to deselect an atom, 
     * that was never selected
     */
    if ( id == -1 ) {
      Tcl_SetResult(interp, "attempt to deselect an atom, that was never selected", TCL_STATIC);
      return TCL_ERROR;
    }
    /* now I must delete sqn ATOM from "coor" and than rewrite the rest 
     * of the "coor" 
     */
    if ( dimType == XC_2D ) {
      printf("tmp_nobjects=%d\n",tmp_nobjects);
      for(i=id; i<tmp_nobjects; i++) {
	ii = i + 1;
	printf("exc:: %d %d\n",i,ii);	
	PushOneCoor(i,ii);
      }
      /* now decrease tmp_nobjects */
      --tmp_nobjects;  
      printf("tmp_nobjects=%d\n",tmp_nobjects);
      fflush(stdout);    
      hpsort_index1(tmp_nobjects, zorient, iwksp);
    }
    else if ( dimType == XC_3D ) {
      printf("before deleting:: id=%d\n",id);
      fflush(stdout);
      /*glDeleteLists(SphereSelList[id-1], nsel - id + 1);*/
      /*xcUpdateSelLists(id);*/
    }      
    /* if we come here, than atom was deselected; decrease nsel */
    --nsel;
  }

  /***************************************************************************/
  /* now "XC_DESELECT <toglName> LINE" 
   * this command deselect all lines
   */
  else if ( strcmp(argv[2],"line") == 0 ) {
    if (argc != 3) {
      Tcl_SetResult(interp, "Usage: xc_deselect <toglName> line", TCL_STATIC);
      return TCL_ERROR;
    }

    if ( dimType == XC_2D ) {
      for(i=old_nobjects; i<=tmp_nobjects; i++) 
	printf("i=%d, flag=%d\n",i,(coor + i)->flag);
      for(i=old_nobjects; i<=tmp_nobjects; i++) {
	printf("i=%d, flag=%d\n",i,(coor + i)->flag);
	if ( (coor + i)->flag == SELLINE || (coor + i)->flag == SELBOND ) {
	  for(j=i; j<tmp_nobjects; j++) { 
	    /* if j+1-objects is SELATOM, then all asel->id's elements
	       that have greater id than j must be decreased for one */
	    if ( (coor + j + 1)->flag == SELATOM )
	      for(jj=1; jj<=nsel; jj++)
		if ((asel + jj)->id > j) (asel + jj)->id -= 1;
	    /* now also Push coor for one down */
	    PushOneCoor(j,j + 1);	
	  }    
	  --i;
	  --tmp_nobjects;	      
	}
      }
      hpsort_index1(tmp_nobjects, zorient, iwksp);
    }
    else if ( dimType == XC_3D ) {
      xcdebug("delete LINES 3D\n");
      /*       for(i=0; i<MAXSEL-1; i++) a */
      /* 	if ( glIsList(LineSelList[i]) ) glDeleteLists(LineSelList[i],1); */
      nline=0;
    }
  }

  else {
    char rss[1024];
    snprintf(rss, sizeof(rss),"unknown option \"%s\" for \"xc_deselect\" command, must be one of atom, line",argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* this is just for now */
  /* (*xcDisplay)(); */

  return TCL_OK;
}


/***************************************************************************/  
/* this function push i-th coor away, and then ii-th goes to i-th place    */
static void 
PushOneCoor(int i, int ii)
{

  (coor + i)->flag  = (coor + ii)->flag;
  (coor + i)->list1 = (coor + ii)->list1;
  (coor + i)->list2 = (coor + ii)->list2;
  (coor + i)->nat   = (coor + ii)->nat;
  (coor + i)->x1    = (coor + ii)->x1;
  (coor + i)->y1    = (coor + ii)->y1;
  (coor + i)->z1    = (coor + ii)->z1;
  (coor + i)->x2    = (coor + ii)->x2;
  (coor + i)->y2    = (coor + ii)->y2;
  (coor + i)->z2    = (coor + ii)->z2;
  *(zorient + i)    = *(zorient + ii);
}


/* when we delete an atom we must update a SphereSelList's lists 
 * from deleted atom on 
 */
/* static void */
/* xcUpdateSelLists(int id) */
/* { */
/*   int sqn, i; */
/*   double r; */

/*   for(i=id-1; i<nsel-1; i++) {     */
/*     sqn = (asel + i + 1)->sqn; */
/*     if (is.spacefillmode) r = atrad[*(nat + sqn)]; */
/*     else if (is.ballmode) r = rball[*(nat + sqn)]; */
/*     else r =  1.5 * rrod; */

/*     glNewList(SphereSelList[i], GL_COMPILE); */
/*       printf("UPD> SphereSelList[i]=%d",SphereSelList[i]); */
/*       fflush(stdout); */
/*       glTranslated( *(xat + sqn), *(yat + sqn), *(zat + sqn) ); */
/*       glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE,  */
/* 		    atcol[SELCOL] ); */
/*       xcSolidSphere( 1.1 * r ); */
/*       glTranslated( -*(xat + sqn), -*(yat + sqn), -*(zat + sqn) ); */
/*     glEndList(); */
/*   } */

/*   /\*  (*xcDisplay)(); *\/ */
/* } */

/*
  static void
  SelectDebug(struct Togl *togl) 
  {
  glColor3f(1.0, 1.0, 0.0);
  if ( dimType == XC_3D ) 
  glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, atcol[0] );
  glPointSize(20.0);
  glBegin(GL_POINTS);
  glVertex2f(0.0, 0.0);
  glVertex2f(sel_debug.xx, sel_debug.yy);    
  glEnd();
  Togl_SwapBuffers(togl);
  }
*/
