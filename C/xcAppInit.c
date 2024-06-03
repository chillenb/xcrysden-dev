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
* Source: $XCRYSDEN_TOPDIR/C/xcAppInit.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
* ------                                                                    *
* Modified by Eric Verfaillie ericverfaillie@yahoo.fr EV                    *
* may 2004                                                                  *
* modifcations are near EV comments                                         *
*****************************************************************************

*/

#include <togl.h>
#include <tk.h>
#include <tcl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "struct.h"
#include "system.h"
#include "memory.h"
#include "xcfunc.h"

#if defined(TOGL_X11)
#   include "GL/glx.h"
#endif

#define XSF_OPEN   0
#define XSF_UPDATE 1

/* on some systems this is needed; I don't know why ????? */
/* #define Tk_MainEx Tk_Main */

extern AtomicLabel *atomLabel, globalAtomLabel;
extern short *do_not_display_atomlabel;
extern XCfog fog;
extern GLfloat def_fog_color[4];
extern GLfloat def_xyz_axis_color[4];
extern GLfloat def_xyz_xyplane_color[4];
extern XYZ_Attrib xyz;
extern ForceVector FV;
extern struct Togl *mesa_togl;


int togl_exists;

/* name of photo image to be used for printing */
const char *printImage = "dump2ppm.82617";

/*===========================================================================*/
/*                THIS IS FOR 3D object (Spheres, Cylinders)                 */
/*===========================================================================*/

Options3D is = { GL_FALSE, GL_FALSE, GL_FALSE, GL_TRUE, GL_TRUE, GL_TRUE };


/* ====================================================== */
/* ---------------function prototypes-------------------- */
/* ====================================================== */

/* implemented Tcl/Tk commands; this file */
int Tcl_AppInit(Tcl_Interp *interp);
int Xcrys_Init( Tcl_Interp *interp );

int XC_OpenStrCmd(ClientData clientData,Tcl_Interp *interp,
		  int argc, const char *argv[]);

int XC_CloseStrCmd(ClientData clientData,Tcl_Interp *interp,
		   int argc, const char *argv[]);

int XC_DisplayMode2DCmd(ClientData clientData, Tcl_Interp *interp,
			int argc, const char *argv[]);

int XC_RotateCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);

int XC_TranslateCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);


int XC_DisplayMode3DCmd(ClientData clientData, Tcl_Interp *interp,
			int argc, const char *argv[]);

int XC_DrawStyle3DCmd(ClientData clientData, Tcl_Interp *interp,
		      int argc, const char *argv[]);

int XC_ShadeModel3DCmd(ClientData clientData, Tcl_Interp *interp,
		       int argc, const char *argv[]);

int XC_PointSizeCmd(ClientData clientData, Tcl_Interp *interp,
		    int argc, const char *argv[]);

int XC_ResetVarCmd(ClientData clientData, Tcl_Interp *interp,
		   int argc, const char *argv[]);

int XC_NewValueCmd(ClientData clientData, Tcl_Interp *interp,
		   int argc, const char *argv[]);

int XC_OldAtmColCmd(ClientData clientData, Tcl_Interp *interp,
                    int argc, const char *argv[]);

int XC_GetDefaultCmd(ClientData clientData, Tcl_Interp *interp,
		     int argc, const char *argv[]);

int XC_GetValueCmd(ClientData clientData, Tcl_Interp *interp,
		   int argc, const char *argv[]);

/* int XC_MesaContextCmd(ClientData clientData, Tcl_Interp *interp,
   int argc, const char *argv[]); */

int XC_UpdateStrCmd(ClientData clientData, Tcl_Interp *interp,
		    int argc, const char *argv[]);

int XC_DisplayCmd(ClientData clientData, Tcl_Interp *interp,
		  int argc, const char *argv[]);

int XC_SwapBufferCmd(ClientData clientData, Tcl_Interp *interp,
		     int argc, const char *argv[]);

int XC_RotationMatrixCmd(ClientData clientData, Tcl_Interp *interp,
			 int argc, const char *argv[]);
int XC_TranslParamCmd(ClientData clientData, Tcl_Interp *interp,
		      int argc, const char *argv[]);

/* extern int XC_MesaWinCmd(ClientData clientData, Tcl_Interp *interp, 
   int argc, const char *argv[]);*/

extern int XC_SelectCmd(ClientData clientData, Tcl_Interp *interp, 
			int argc, const char *argv[]);

extern int XC_DeselectCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, const char *argv[]);

extern int XC_AtomAddCmd(ClientData clientData, Tcl_Interp *interp,
			 int argc, const char *argv[]);

/* implemented xc_iso* Tcl/Tk commands; xcIsoSurf.c */
extern int XC_IsoCmd(ClientData clientData, Tcl_Interp *interp,
		     int argc, const char *argv[]);

extern int XC_IsostackCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, const char *argv[]);

extern int XC_IsosignCmd(ClientData clientData, Tcl_Interp *interp,
			 int argc, const char *argv[]);

extern int XC_IsofilesCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, const char *argv[]);

extern int XC_IsopointsCmd(ClientData clientData, Tcl_Interp *interp,
			   int argc, const char *argv[]);

extern int XC_IsodataCmd(ClientData clientData, Tcl_Interp *interp,
			 int argc, const char *argv[]);

extern int XC_IsosurfCmd(ClientData clientData, Tcl_Interp *interp,
			 int argc, const char *argv[]);

extern int XC_IsoplaneCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, const char *argv[]);

extern int XC_IsoexpandCmd(ClientData clientData, Tcl_Interp *interp,
			   int argc, const char *argv[]);

extern int XC_SetGLparamCmd(ClientData clientData, Tcl_Interp *interp,
			    int argc, const char *argv[]);

extern int XC_GetGLparamCmd(ClientData clientData, Tcl_Interp *interp,
			    int argc, const char *argv[]);

extern int XC_BzCmd(ClientData clientData,Tcl_Interp *interp,
		    int argc, const char *argv[]);

extern int XC_SuperCellCmd(ClientData clientData,Tcl_Interp *interp,
			   int argc, const char *argv[]);

extern int XC_IsoSpaceSelCmd(ClientData clientData, Tcl_Interp *interp,
			     int argc, const char *argv[]);

extern int XC_WignerCmd(ClientData clientData,Tcl_Interp *interp,
			int argc, const char *argv[]);

extern int XC_IsoDataGridCmd(ClientData clientData, Tcl_Interp *interp,
			     int argc, const char *argv[]);

extern int XC_ColorschemeCmd(ClientData clientData, Tcl_Interp *interp,
			     int argc, const char *argv[]);

extern int XC_ReadXSFCmd(ClientData clientData, Tcl_Interp *interp,
			 int argc, const char *argv[]);

extern int XC_ReadBandXSFCmd(ClientData clientData, Tcl_Interp *interp,
			     int argc, const char *argv[]);

extern int XC_GridValueCmd(ClientData clientData, Tcl_Interp *interp,
			   int argc, const char *argv[]);

extern int XC_F3toI4Cmd(ClientData clientData, Tcl_Interp *interp,
			int argc, const char *argv[]);

extern int XC_FractCoorCmd(ClientData clientData, Tcl_Interp *interp,
			   int argc, const char *argv[]);

extern int XC_MolSurfCmd(ClientData clientData, Tcl_Interp *interp,
			 int argc, const char *argv[]);

extern int XC_MolSurfRegCmd(ClientData clientData, Tcl_Interp *interp,
			    int argc, const char *argv[]);
extern int XC_MolSurfConfigCmd(ClientData clientData, Tcl_Interp *interp,
			       int argc, const char *argv[]);

extern int CRY_DispFuncCmd(ClientData clientData,Tcl_Interp *interp,
			   int argc, const char *argv[]);

extern int CRY_DispFuncMultiFSCmd(ClientData clientData,Tcl_Interp *interp,
				  int argc, const char *argv[]);

extern int CRY_SurfRegCmd(ClientData clientData,Tcl_Interp *interp,
			  int argc, const char *argv[]);

extern int CRY_SurfConfigCmd(ClientData clientData,Tcl_Interp *interp,
			     int argc, const char *argv[]);

extern int CRY_SurfCmd(ClientData clientData,Tcl_Interp *interp,
		       int argc, const char *argv[]);

extern int XC_ForcesCmd(ClientData clientData,Tcl_Interp *interp, int argc, const char *argv[]);

extern int XC_HBondsCmd(ClientData clientData,Tcl_Interp *interp, int argc, const char *argv[]);

extern int XC_StereoCmd(ClientData clientData,Tcl_Interp *interp, int argc, const char *argv[]);

/* END:: implemented xc_iso* Tcl/Tk commands; xcIsoSurf.c */


/* ppmPrintTogl.c */
extern int CRY_Dump2PpmCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
/* togl_ppm.c */
extern int Togl_DumpToPpmFile(struct Togl *togl, const char *filename);
  
/* realTimeMovie.c */
extern int CRY_RealTimeMovieCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);

/* gl2psPrintTogl.c */
extern int CRY_gl2psPrintToglCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);

/* cryTogl */
extern int CRY_ToglZoomCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);

/* still this file */
int xcRotateX(double fi);
int xcRotateY(double fi);
int xcRotateZ(double fi);
int xcRotateXY(double fiX, double fiY);
void FreeAllVariables(void);

/* custom event handler */
/* extern void XC_MesaEvent(ClientData clientData, XEvent *eventPtr);*/

/* --- auxilary functions; auxilproc.c --- */
extern void Rotate(double *x, double *y, double *z, double cosfi, double sinfi);
/* --- xcballstick.c --- */
extern void MakeArcPoints(void);

/* --- readstrf.c --- */
extern int ReadStructFile(FILE *fp, const char *file, int format, int mode );
/* --- viewport.c --- */
extern void xcViewPort(void);

/* --- loadValues --- */
extern void LoadDefaultValues(void);
extern void ResetVar(struct Togl *togl, int var);
extern void ReloadVars(void);
extern void LoadNewValue(struct Togl *togl, int var, double value1, 
			 double value2, double value3, double value4);
extern double GetDefault(int var);
extern double GetValue(int var);
extern void LoadOldAtomicColors(void);

/* --- xcMesaContext.c --- */
/* GLenum xcCreateMesaContext(Display *xcDisplay, Window xcMesaWin, 
   int xcscreen);*/
/* void xcDestroyMesaContext(void); */

/* --- 3D.c --- */
extern void From2Dto3D(void);
extern void GetOldMat(void);
extern void GetMajorMat(void);
extern void GetRotX(GLdouble fiX);
extern void GetRotY(GLdouble fiY);
extern void GetRotZ(GLdouble fiZ);
extern void LoadIdentity(GLdouble matrix[][4]);

/* --- vectors.c --- */
extern void vecMatToVec(double major[4][4], double vec[16]);
extern void VecRotateX(double cosfi, double sinfi);
extern void VecRotateY(double cosfi, double sinfi);
extern void VecRotateZ(double cosfi, double sinfi);
extern void VecRotateXY(double cosfiX, double sinfiX, 
			double cosfiY, double sinfiY);
extern void  VecRotTo_XY(void);
extern void  VecRotTo_XZ(void);
extern void  VecRotTo_YZ(void);
extern void  VecRotTo_AB(void);
extern void  VecRotTo_AC(void);
extern void  VecRotTo_BC(void);

/* --- xcLabels.c --- */
extern void makeAtomLabels(void);
extern void makeXYZLabels(void);
extern void makeRasterFont(void);
extern void makeTemp3D2DList(void);
extern void makeCrdList(void);

/*****************************************************************************/
/*                          MESA  --  MESA  --  MESA                         */
/*****************************************************************************/
/* --- xcDisplayFunc.c --- */
extern void (*xcDisplay)(struct Togl *togl);
extern void xcGenDispList(void);
extern void xcDisplayFunc( void (*Func)(struct Togl *togl) );
extern int xcToglDisplayFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv);
extern void xcDummyDisplay( struct Togl *togl );
extern void xcWireFrame2D(struct Togl *togl);
extern void xcPointLine2D(struct Togl *togl);
extern void xcBallStick12D(struct Togl *togl);
extern void xcBallStick22D(struct Togl *togl);
extern void xcDisplay3D(struct Togl *togl);
extern void EnableOr2D_Or3D(void);
extern void xcMakeProjection2D(const char *mode);
extern void xcMakeProjection3D(const char *mode);
extern void xcClearScreen(struct Togl *togl);
/* extern void xcMakeBallLists(void); */
extern void xcMakePointList();
/*extern void xcMaybeDestroyLists(void);*/

/* --- xcDisplayFunc2.c --- */
extern void xcAssignDisplayFunc(const char *dispMode);
extern void UpdateDispFunc(void);
extern void RewriteCoor(GLenum type);

/* --- lighting.c --- */
extern void LoadStructMaterial(void);
extern void LoadLights(void);
extern void LoadBlendfunc_And_Frontface(void);

/* --- cells.c --- */
extern void CellTypes(void);

/* --- xcTogl.c --- */
extern int xcToglCreateFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv);
extern int xcToglReshapeFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv);
extern int xcToglDestroyFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv);
extern int xcToglTimerFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv);
extern int XC_B1MotionCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_ShiftB1MotionCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_B2MotionCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_ButtonReleaseCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);

/* --- xcFont.c --- */
extern int XC_SetFontCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_SetAtomLabelCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_ClearAtomLabelCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_QueryFontCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);


/* --- sInfo.c --- */
extern void Set_sInfoArray( Tcl_Interp *interp );

/* --- datagrid.c --- */
extern void CloseGridList(void);

/* xcdebug.c */
extern void xcdebug(const char *text);
extern void xcErrDebug(const char *text);
extern void breakpoint(const char *text);

/* --- pia.c ---*/
/*extern void Pia(void);*/

#if defined _LINUXALPHA || defined _DEC_CUSTOM_FPE_HANDLER
/* --- signal.c --- */
extern void xcFPEHandler(void);
#endif

/* cryNewContext.c */
extern NEW_WIN_CONTEXT *FindWinContextByTogl(struct Togl *togl);



#ifndef WIN32 
static int 
XErrorFunction(Display *display,XErrorEvent *event)
{
  char buffer[512];

  buffer[0]='\0';
  XGetErrorText(display,event->error_code,buffer,sizeof(buffer)-1);
  buffer[511]='\0';
  fprintf(stderr,"X Error of failed request: %s\n",buffer);
  fprintf(stderr,"    Major opcode of failed request: %d\n",event->request_code);
  switch(event->error_code) {
  case BadValue:
    fprintf(stderr,"    Value in failed request: 0x%lx\n",event->resourceid);
    break;
  case BadAtom:
    fprintf(stderr,"    AtomID in failed request: 0x%lx\n",event->resourceid);
    break;
  default:
    if (event->error_code != BadRequest)
      fprintf(stderr,"    ResourceID in failed request: 0x%lx\n",event->resourceid);
    break;
  }
  fprintf(stderr,"    Serial number of failed request: %ld\n",event->serial);
  return 0;
}

/*****************************************************************************/
int
main(int argc, char *argv[])  
{  
  /* create signal handler for FPE for LinuxAlpha */
#ifdef _LINUXALPHA
  fprintf(stderr,"Note: *** custom SIG_FPE signal handler is used for Linux-Alpha platform\n", NULL);
  xcFPEHandler();
#endif
#ifdef _DEC_CUSTOM_FPE_HANDLER
  fprintf(stderr,"Note: *** custom SIG_FPE signal handler will be used !!!\n", NULL);
  xcFPEHandler();
#endif

  XSetErrorHandler(XErrorFunction);
  Tk_Main(argc, argv, Tcl_AppInit);
  exit(0);
}


/*
 * Tcl_AppInit is called from Tcl_Main
 * after the Tcl interpreter has been created,
 * and before the script file
 * or interactive command loop is entered.
 */
int
Tcl_AppInit(Tcl_Interp *interp) {
  /*
   * Initialize packages
   * Tcl_Init sets up the Tcl library facility.
   */
  if (Tcl_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }
  if (Tk_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }  

  if (Xcrys_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }

  /* HERE:
   * Define startup filename if any. This file is read in
   * case the program is run interactively.
   */

  return TCL_OK;
}

#endif /* not WIN32 */


#ifdef WIN32
__declspec(dllexport)
#endif
int 
Xcrys_Init( Tcl_Interp *interp )
{
  char *dir = NULL;

  /*
   * Define XCRYSDEN application-specific commands here.
   * 
   */
  int i;
  struct CUSTOM_COM {
    Tcl_CmdProc *proc;
    char        *name;
  } custom_command[] = { 
    { (Tcl_CmdProc *)XC_ForcesCmd,         "xc_forces"},
    { (Tcl_CmdProc *)XC_GridValueCmd,      "xc_gridvalue"},
    { (Tcl_CmdProc *)CRY_DispFuncCmd,      "cry_dispfunc"},
    { (Tcl_CmdProc *)CRY_DispFuncMultiFSCmd, "cry_dispfuncmultiFS"},
    { (Tcl_CmdProc *)CRY_SurfRegCmd,	    "cry_surfreg"},     
    { (Tcl_CmdProc *)CRY_SurfConfigCmd,    "cry_surfconfig"},		 
    { (Tcl_CmdProc *)CRY_SurfCmd,          "cry_surf"},
    { (Tcl_CmdProc *)XC_RotationMatrixCmd, "xc_rotationmatrix"},
    { (Tcl_CmdProc *)XC_TranslParamCmd,    "xc_translparam"},
    { (Tcl_CmdProc *)XC_WriteXSFCmd,       "xc_writeXSF" },
    { (Tcl_CmdProc *)XC_HBondsCmd,         "xc_hbonds"},
    { (Tcl_CmdProc *)XC_WriteBandXSFCmd,   "xc_writebandXSF"},
    { (Tcl_CmdProc *)XC_StereoCmd,         "xc_stereo"},
    { NULL, "" }
  };

  struct TOGL_CMD {
    Tcl_CmdProc *proc;
    char        *name;
  } togl_cmd[] = {    
    { (Tcl_CmdProc *)XC_RotateCb,          "xc_rotate" },       
    { (Tcl_CmdProc *)XC_TranslateCb,       "xc_translate" },    
    { (Tcl_CmdProc *)XC_B1MotionCb,        "xc_B1motion" },     
    { (Tcl_CmdProc *)XC_ShiftB1MotionCb,   "xc_ShiftB1motion" },
    { (Tcl_CmdProc *)XC_B2MotionCb,        "xc_B2motion" },     
    { (Tcl_CmdProc *)XC_ButtonReleaseCb,   "xc_Brelease" },     
    { (Tcl_CmdProc *)CRY_gl2psPrintToglCb, "cry_gl2psPrintTogl" },
    { (Tcl_CmdProc *)CRY_Dump2PpmCb,       "cry_dump2ppm" },
    { (Tcl_CmdProc *)CRY_ToglZoomCb,       "cry_toglzoom" },
    { (Tcl_CmdProc *)CRY_RealTimeMovieCb,  "xc_realtimemovie" },
    { (Tcl_CmdProc *)XC_SetFontCb,         "xc_setfont" },
    { (Tcl_CmdProc *)XC_SetAtomLabelCb,    "xc_setatomlabel" },
    { (Tcl_CmdProc *)XC_ClearAtomLabelCb,  "xc_clearatomlabel" },
    { (Tcl_CmdProc *)XC_QueryFontCb,       "xc_queryfont" },    
    { NULL, "" }
  };

  /* some XC_initializations */
  LoadDefaultValues();
  CellTypes();
  xcr.prim_nat   = 0;
  xcr.prim_coor  = 0;
  xcr.prim_fcoor = 0;
  xcr.prim_forc  = 0; 
  xcr.conv_nat   = 0;
  xcr.conv_fcoor = 0;

  xc_system.pid         = getpid();
  xc_system.scratch_dir = (char *) malloc ( sizeof(char) * 2048 );

  atm.sqn = atm.sqnat = NULL;
  atm.col = NULL;

  /* $system(SCRDIR) = $XCRYSDEN_SCRATCH/xc_$$ */
  /* get some minimum of environmental data */
  if ( (dir = getenv("XCRYSDEN_SCRATCH")) != NULL )
    sprintf(xc_system.scratch_dir,"%s/xc_%d/", dir, xc_system.pid);

  MakeArcPoints();

  /*
   * TOGL stuff ...
   */ 
  if (Togl_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }  
  /* Togl Setup & Initialization functions */
  Tcl_CreateObjCommand(interp, "togl_create",  xcToglCreateFunc,  (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateObjCommand(interp, "togl_display", xcToglDisplayFunc, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateObjCommand(interp, "togl_reshape", xcToglReshapeFunc, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateObjCommand(interp, "togl_destroy", xcToglDestroyFunc, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
  Tcl_CreateObjCommand(interp, "togl_timer",   xcToglTimerFunc,   (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  

  /*
   * XCRYSDEN application-specific Tcl/Tk commands
   */

  Tcl_CreateCommand(interp, "xc_openstr", (Tcl_CmdProc *)XC_OpenStrCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_closestr", (Tcl_CmdProc *)XC_CloseStrCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_displayMode2D", (Tcl_CmdProc *)XC_DisplayMode2DCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_displayMode3D", (Tcl_CmdProc *)XC_DisplayMode3DCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_drawStyle3D", (Tcl_CmdProc *)XC_DrawStyle3DCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_shadeModel3D", (Tcl_CmdProc *)XC_ShadeModel3DCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_pointSize", (Tcl_CmdProc *)XC_PointSizeCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  /*
    Tcl_CreateCommand(interp, "xc_mesawin", (Tcl_CmdProc *)XC_MesaWinCmd, 
    (ClientData)Tk_MainWindow(interp), 
    (Tcl_CmdDeleteProc *)NULL);

    Tcl_CreateCommand(interp, "xc_mesacontext", (Tcl_CmdProc *)XC_MesaContextCmd, 
    (ClientData)Tk_MainWindow(interp), 
    (Tcl_CmdDeleteProc *)NULL);
  */

  Tcl_CreateCommand(interp, "xc_resetvar", (Tcl_CmdProc *)XC_ResetVarCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_newvalue", (Tcl_CmdProc *)XC_NewValueCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_oldatmcol", (Tcl_CmdProc *)XC_OldAtmColCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_getdefault", (Tcl_CmdProc *)XC_GetDefaultCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_getvalue", (Tcl_CmdProc *)XC_GetValueCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_updatestr", (Tcl_CmdProc *)XC_UpdateStrCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_swapbuffer", (Tcl_CmdProc *)XC_SwapBufferCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_display", (Tcl_CmdProc *)XC_DisplayCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_select", (Tcl_CmdProc *)XC_SelectCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_deselect", (Tcl_CmdProc *)XC_DeselectCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_atomadd", (Tcl_CmdProc *)XC_AtomAddCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_iso", (Tcl_CmdProc *)XC_IsoCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_isostack", (Tcl_CmdProc *)XC_IsostackCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_isosign", (Tcl_CmdProc *)XC_IsosignCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_isofiles", (Tcl_CmdProc *)XC_IsofilesCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_isopoints", (Tcl_CmdProc *)XC_IsopointsCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_isodata", (Tcl_CmdProc *)XC_IsodataCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_isosurf", (Tcl_CmdProc *)XC_IsosurfCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_isoplane", (Tcl_CmdProc *)XC_IsoplaneCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_isoexpand", (Tcl_CmdProc *)XC_IsoexpandCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_setGLparam", (Tcl_CmdProc *)XC_SetGLparamCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_getGLparam", (Tcl_CmdProc *)XC_GetGLparamCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_bz", (Tcl_CmdProc *)XC_BzCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_supercell", (Tcl_CmdProc *)XC_SuperCellCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_isospacesel", (Tcl_CmdProc *)XC_IsoSpaceSelCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_wigner", (Tcl_CmdProc *)XC_WignerCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_isodatagrid", (Tcl_CmdProc *)XC_IsoDataGridCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_colorscheme", (Tcl_CmdProc *)XC_ColorschemeCmd, 
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_readXSF", (Tcl_CmdProc *)XC_ReadXSFCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_readbandXSF", (Tcl_CmdProc *)XC_ReadBandXSFCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_f3toi4", (Tcl_CmdProc *)XC_F3toI4Cmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_fractcoor", (Tcl_CmdProc *)XC_FractCoorCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_molsurf", (Tcl_CmdProc *)XC_MolSurfCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_molsurfreg", (Tcl_CmdProc *)XC_MolSurfRegCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "xc_molsurfconfig", (Tcl_CmdProc *)XC_MolSurfConfigCmd,
		    (ClientData)Tk_MainWindow(interp), 
		    (Tcl_CmdDeleteProc *)NULL);

  i=-1;
  while( custom_command[++i].proc ) 
    Tcl_CreateCommand(interp, 
		      custom_command[i].name,
		      custom_command[i].proc,
		      (ClientData) Tk_MainWindow(interp),
		      (Tcl_CmdDeleteProc *) NULL);

  /*
   * application-specific Togl sub-commands
   */
  i=-1;
  while( togl_cmd[++i].proc ) {
    Tcl_CreateCommand(interp,
		      togl_cmd[i].name,
		      togl_cmd[i].proc,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
  }
  
  /* define xcrys(platform) variable */
#if defined(MAC_OSX)
  Tcl_SetVar2(interp, "xcrys", "platform", "macosx", TCL_GLOBAL_ONLY);
#elif defined(CYGWIN)
  Tcl_SetVar2(interp, "xcrys", "platform", "cygwin", TCL_GLOBAL_ONLY);
#elif defined(WIN32)
  Tcl_SetVar2(interp, "xcrys", "platform", "windows", TCL_GLOBAL_ONLY);
#else
  Tcl_SetVar2(interp, "xcrys", "platform", "unix", TCL_GLOBAL_ONLY);
#endif

  /* store the photo image to be used for printing into xcrys(print.image) */
  Tcl_SetVar2(interp, "xcrys", "print.image", printImage, TCL_GLOBAL_ONLY);
  
  togl_exists = XC_FALSE;

  return TCL_OK;
}


/* XC_OpenStr --> inplementation of 'xc_openstr' custom Tcl command
 * ---------------------
 * Usage: xc_openstr <fileformat> <file> <toglName> ?DisplayMode?
 *
 *                 fileformat   --- format of file (xcr, xyz, pdb)
 *                 
 *                 toglName     --- name of togl where structure will be 
 *                                  displayed
 *
 *                 DisplayMode  --- one of WF, PL, BS1 and BS2
 *
 * NOTE: so far we can open structure and display only in 2D modes;
 *       and only than we can switch to 3D modes
 */
int 
XC_OpenStrCmd(ClientData clientData, Tcl_Interp *interp,
	      int argc, const char *argv[])
{
  struct Togl *togl;
  FILE *fp;
  int format;
  /* vars. for MESA */

  if (argc < 4 || argc > 5) {
    Tcl_SetResult(interp, "Usage: xc_openstr <fileformat> <file> <toglName> ?DisplayMode?", TCL_STATIC);
    return TCL_ERROR;
  }
  /* T.K.__correcture */
  if (argc == 4) argv[4] = "WF"; 

  /* determine format of file */
  if ( (strcmp(argv[1],"xcr") == 0) || (strcmp(argv[1],"xsf") == 0) ) format = FORMAT_XSF;
  else if ( strcmp(argv[1],"xyz") == 0 ) format = FORMAT_XYZ;
  else if ( strcmp(argv[1],"pdb") == 0 ) format = FORMAT_PDB;
  else {
    char rss[1024];
    snprintf(rss, sizeof(rss),"unknown format \"%s\", must be one of \"xcr\", \"xyz\", \"pdb\"", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* Check if 'file' is OK */
  fp = fopen(argv[2],"r");
  if (fp == NULL) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "Can't open file \"%s\"",argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* find togl associated with toglName */
  if ( Togl_GetToglFromName(interp, argv[3], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* Check if DisplayMode is OK */
  if ( argc == 5 && strcmp(argv[4],"WF") != 0 
       && strcmp(argv[4],"PL") != 0 && strcmp(argv[4],"BS1") != 0 
       && strcmp(argv[4],"BS2") != 0 ) 
    {
      char rss[1024];
      snprintf(rss, sizeof(rss), "Unknown DisplayMode \"%s\". Must be one of WF, PL, BS1, BS2",argv[4]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    } 

  /* is structure already opened ??? */
  if (VPf.stropened) {
    const char *arg[] = { "xc_closestr" };
    XC_CloseStrCmd(clientData, interp, 1, arg);
  }

  /* reset all variables needed */
  ReloadVars();

  /* command is OK 
   * LET'S READ A FILE
   */
  if (!ReadStructFile(fp, argv[2], format, XSF_OPEN)) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"Error reading file \"%s\"",argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    fclose(fp);
    return TCL_ERROR;
  }

  /* structure is now opened, so we must assign ... */
  VPf.stropened = 1;
  /* now we can close the file */
  fclose(fp);

  /******************************************************/
  /*        now initialize various vectors              */
  /******************************************************/
  LoadIdentity( vec.crdmajor ); /* coor-sist base vectors */
  vecMatToVec( vec.crdmajor, vec.crdvec );

  /******************************************************/
  /* now display the structure according to DisplayMode */
  /******************************************************/
  /* because displayMode is one of 2D */
  dimType = XC_2D;
  /* enable things for 2D */
  EnableOr2D_Or3D();

  /* now make xcViewPort with VPf.projmade = GL_FALSE */
  VPf.projmade = GL_FALSE;
  xcViewPort();
  printf("after xcViewPort, argv=%s %s %s %s\n",argv[1],argv[2],
	 argv[3],argv[4]);

  xcdebug("before xcMakeProjection2D");
  /* if displaymode is not specified, make WF */
  if ( strcmp(argv[4],"WF") == 0 || argc == 4) {
    xcdebug("before xcMakeProjection2D");
    xcMakeProjection2D("WF");
    xcAssignDisplayFunc("WF");
  }
  if ( strcmp(argv[4],"PL") == 0) {
    printf("Goint to xcPointLineCreate");
    xcMakeProjection2D("PL");
    xcAssignDisplayFunc("PL");
    xcMakePointList();
  }
  if ( strcmp(argv[4],"BS1") == 0 ) {
    xcMakeProjection2D("BS1");
    xcAssignDisplayFunc("BS1");
    /* xcMakeBallLists(); */
  }
  if ( strcmp(argv[4],"BS2") == 0 ) {
    xcMakeProjection2D("BS2");
    xcAssignDisplayFunc("BS2");
    /* xcMakeBallLists(); */
  }

  /* now make xcViewPort with VPf.projmade = GL_TRUE */
  VPf.projmade = GL_TRUE;
  xcViewPort();

  /* Display a Structure */
  Togl_PostRedisplay(togl); 

  /******************************************************/
  /*            now create event handler                */
  /******************************************************/
  /* event_mask = Button1MotionMask | 
     Button2MotionMask | 
     ExposureMask | 
     StructureNotifyMask | 
     ButtonReleaseMask |
     SubstructureNotifyMask |
     FocusChangeMask;
     Tk_CreateEventHandler(tk_mesawin, event_mask,  
     (Tk_EventProc *)XC_MesaEvent, (ClientData) tk_mesawin);
  */

  printf("=>end of XC_OpenStrCmd\n");
  fflush(stdout);

  /* assing the sInfo global array */
  Set_sInfoArray( interp );
  /* if we come so far we are succesfully */

  return TCL_OK;
}


/* XC_DisplayMode2DCmd --> inplementation of 'xc_displayMode2D' custom 
 * Tcl command
 * ---------------------
 * Usage: xc_displayMode2D <toglName> <displayMode> 
 */
int 
XC_DisplayMode2DCmd(ClientData clientData, Tcl_Interp *interp,
		    int argc, const char *argv[])
/* =================================================================
 *     procedure assign wanted display function to xcDisplayFunc
 *               indirectly via "xcAssignDisplayFunc"
 * =================================================================
 */
{ 
  struct Togl *togl;

  /* was xc_wireframe invoke correctly ??*/
  if (argc != 3) {
    Tcl_SetResult(interp, "Usage: xc_displayMode2D <toglName> <displayMode>", TCL_STATIC);
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

  /* is displayMode OK ??? */
  if ( strcmp(argv[2],"WF") != 0 &&
       strcmp(argv[2],"PL") != 0 &&
       strcmp(argv[2],"PB") != 0 &&
       strcmp(argv[2],"BS1") != 0 &&
       strcmp(argv[2],"BS2") != 0 &&
       strcmp(argv[2],"SF") != 0 ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"Unknown <displayMode> \"%s\", must be one of WF, PL, PB, BS1, BS2, SF", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* if structure is not opened, return silently */
  if (!VPf.stropened) return TCL_OK;

  /* if we are in selection mode, mode-changing is not allowed, return 
     silenlty */
  if (VPf.selection) return TCL_OK;

  /* glEnable things for 2D */
  dimType = XC_2D;
  EnableOr2D_Or3D();

  if ( strcmp(argv[2],"WF") == 0 ) {
    sprintf(displayMode2D, "WF");
    xcMakeProjection2D("WF");
    xcAssignDisplayFunc("WF");
  } 
  else if ( strcmp(argv[2],"PL") == 0 ) {
    sprintf(displayMode2D, "PL");
    xcMakeProjection2D("PL");
    xcAssignDisplayFunc("PL");
    xcMakePointList();
  } 
  else if ( strcmp(argv[2],"PB") == 0 ) {
    sprintf(displayMode2D, "PB");
    MakeArcPoints();
    xcMakeProjection2D("PB");
    xcAssignDisplayFunc("PB");
    xcMakePointList();
  } 
  else if ( strcmp(argv[2],"BS1") == 0 ) {
    /* for each atom of a kind make a displaylist for ball */
    sprintf(displayMode2D, "BS1");
    xcMakeProjection2D("BS1");
    xcAssignDisplayFunc("BS1");
    /* xcMakeBallLists(); */
  } 
  else if ( strcmp(argv[2],"BS2") == 0 ) {
    sprintf(displayMode2D, "BS2");
    MakeArcPoints();
    /* for each atom of a kind make a displaylist for ball */
    xcMakeProjection2D("BS2");
    xcAssignDisplayFunc("BS2");
    /* xcMakeBallLists(); */
  } 
  else if ( strcmp(argv[2],"SF") == 0 ) {
    sprintf(displayMode2D, "SF");
    /* for each atom of a kind make a displaylist for ball */
    xcMakeProjection2D("SF");
    xcAssignDisplayFunc("SF");
    /* xcMakeBallLists(); */
  } 
  else {
    char rss[1024];
    snprintf(rss, sizeof(rss),"unknown display mode \"%s\", must be one of WF, PL,PB, BS1, BS2, SF\n",argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  /* now perform Z-orientation */
  hpsort_index1(tmp_nobjects,zorient, iwksp);    
  /* update a display */
  Togl_PostRedisplay(togl); 

  return TCL_OK;
}


/* XC_RotateCb --> implementation of 'xc_rotate' custom Togl sub-command 
 * ---------------------------------------------------------------------  
 * 
 * Usage: xc_rotate <toglName> dir rotstep
 *
 *        dir     --  +x | -x | +y | -y | +z | -z | ++xy | ++xz | ++yz | xy | xz | yz | ab | ac | bc
 *                
 *        rotstep --  rotation step (in degrees)           
 *
 * NOTE: this command is invoked by ButtonWidget-press;
 *       mouse driven rotation is now implemented by XC_B1MotionCmd
 */
int 
XC_RotateCb(ClientData clientData, Tcl_Interp *interp,
	    int argc, const char *argv[])
{
  Togl *togl;
  double rotstep, rotstep1;

  /* this proc is meant only for .mesa widget */
  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( togl != mesa_togl) {
    Tcl_SetResult(interp,"xc_rotate is meant only for .mesa",TCL_STATIC);
    return TCL_ERROR;
  }

  /* if structure is not opened, do noting, just return silently */
  if (!VPf.stropened) return TCL_OK;

  if (argc < 4 || argc > 5) {
    Tcl_SetResult(interp, "Usage: xc_rotate <toglName> <dir> <rotstep> ?<rotstep1>?", TCL_STATIC);
    return TCL_ERROR;
  }
  /*printf("Begin of XC_RotateCmd\n",NULL);*/

  /* is rotstep valid enough */
  if (Tcl_GetDouble(interp, argv[3], &rotstep) == TCL_ERROR) {
    Tcl_AppendResult(interp,"invalid rotstep \"%s\" in %s xc_rotate command", 
		     argv[0], argv[3],(char *)NULL);
    return TCL_ERROR;
  }    
  /* from degrees to rad */
  rotstep /= RAD2DEG;

  if (strncmp(argv[2],"++",2) == 0 ) {
    if (Tcl_GetDouble(interp, argv[4], &rotstep1) == TCL_ERROR) {
      Tcl_AppendResult(interp,
		       "invalid rotstep \"%s\" in %s xc_rotate command", 
		       argv[0], argv[3],(char *)NULL);
      return TCL_ERROR;
    }
    rotstep1 /= RAD2DEG;
  }

  if ( dimType == XC_2D ) {
    /* this is now rotation only for 2D display structures */
    /* when Z-rotating we don't need to "hpsort_index1" */
    if (argv[2][1] == 'z') {
      if (argv[2][0] == '-') xcRotateZ(-rotstep);
      else if (argv[2][0] == '+') xcRotateZ( rotstep);
      else if (argv[2][0] == 'x') {
	VecRotTo_XZ();
	xcAssignDisplayFunc(displayMode2D);
      }
      else if (argv[2][0] == 'y') {
	VecRotTo_YZ();
	xcAssignDisplayFunc(displayMode2D);
      }	
      else {
	Tcl_SetResult(interp, "invalid rotation direction in xc_rotate command, must be one of +x, -x, +y, -y, +z, -z, xy, xz, yz", TCL_STATIC);
	return TCL_ERROR;
      }
    } else {
      /* X & Y rotation */
      /* if X || Y rotation, we must "hpsort_index1" */
      /* because Y-coord is in Xlib upside down and in GL it is OK,
	 I'm must change +X/-X rotation to made things OK */
      if ( strcmp(argv[2],"+x") == 0 )     xcRotateX(-rotstep);
      else if (strcmp(argv[2],"-x") == 0 ) xcRotateX( rotstep);
      else if (strcmp(argv[2],"+y") == 0 ) xcRotateY( rotstep);
      else if (strcmp(argv[2],"-y") == 0 ) xcRotateY(-rotstep);
      else if (strcmp(argv[2],"++xy") == 0 ) {
	xcRotateX(rotstep);
	xcRotateY(rotstep1);
      }
      else if (strcmp(argv[2],"++xz") == 0 ) {
	xcRotateX(rotstep);
	xcRotateZ(rotstep1);
      }
      else if (strcmp(argv[2],"++yz") == 0 ) {
	xcRotateY(rotstep);
	xcRotateZ(rotstep1);
      }    
      else if (strcmp(argv[2],"xy") == 0 ) {
	VecRotTo_XY();
	xcAssignDisplayFunc(displayMode2D);
      }
      else if (strcmp(argv[2], "ab" ) == 0 ) {
	VecRotTo_AB();
	xcAssignDisplayFunc(displayMode2D);
      }
      else if (strcmp(argv[2], "ac" ) == 0 ) {
	VecRotTo_AC();
	xcAssignDisplayFunc(displayMode2D);
      }
      else if (strcmp(argv[2], "bc" ) == 0 ) {
	VecRotTo_BC();
	xcAssignDisplayFunc(displayMode2D);
      }      
      else {
	Tcl_SetResult(interp, "invalid rotation direction in xc_rotate command, must be one of +x, -x, +y, -y, +z, -z, ++xy, ++xz, ++yz, xy, xz, yz, ab, ac, bc", TCL_STATIC);
	return TCL_ERROR;
      } /* T.K. */
      hpsort_index1(tmp_nobjects, zorient, iwksp);       
    }
  } 
  else if ( dimType == XC_3D ) 
    {
      double cosfi, sinfi;
      /* because Y-coord is in Xlib upside down and in GL it is OK,
	 I'm must change +X/-X rotation to made things OK */
      if ( strcmp(argv[2],"+x") == 0 )     {
	cosfi = cos(-rotstep);
	sinfi = sin(-rotstep);
	GetRotX(-rotstep);
	VecRotateX( cosfi, sinfi);
      }
      else if (strcmp(argv[2],"-x") == 0 ) {
	cosfi = cos( rotstep);
	sinfi = sin( rotstep);
	GetRotX( rotstep);
	VecRotateX( cosfi, sinfi);
      }
      else if (strcmp(argv[2],"+y") == 0 ) {
	cosfi = cos( rotstep);
	sinfi = sin( rotstep);
	GetRotY( rotstep);
	VecRotateY( cosfi, sinfi);
      }
      else if (strcmp(argv[2],"-y") == 0 ) {
	cosfi = cos(-rotstep);
	sinfi = sin(-rotstep);
	GetRotY(-rotstep);
	VecRotateY( cosfi, sinfi);
      }
      else if (strcmp(argv[2],"+z") == 0 ) {
	cosfi = cos( rotstep);
	sinfi = sin( rotstep);
	GetRotZ( rotstep);
	VecRotateZ( cosfi, sinfi);
      }
      else if (strcmp(argv[2],"-z") == 0 ) {
	cosfi = cos(-rotstep);
	sinfi = sin(-rotstep);
	GetRotZ(-rotstep);
	VecRotateZ( cosfi, sinfi);
      }
      else if (strcmp(argv[2],"++xy") == 0 ) {
	cosfi = cos( rotstep);
	sinfi = sin( rotstep);
	GetRotX( rotstep);
	VecRotateX( cosfi, sinfi);
	/*GetOldMat();			
	  GetMajorMat();*/

	cosfi = cos(rotstep1);
	sinfi = sin(rotstep1);
	GetRotY( rotstep1 );
	VecRotateY( cosfi, sinfi);
      }
      else if (strcmp(argv[2],"++xz") == 0 ) {
	cosfi = cos (rotstep);
	sinfi = sin (rotstep);
	GetRotX (rotstep);
	VecRotateX (cosfi, sinfi);

	cosfi = cos( rotstep1);
	sinfi = sin( rotstep1);
	GetRotZ( rotstep1 );
	VecRotateZ( cosfi, sinfi);
      }
      else if (strcmp(argv[2],"++yz") == 0 ) {
	cosfi = cos( rotstep);
	sinfi = sin( rotstep);
	GetRotY( rotstep);
	VecRotateY( cosfi, sinfi);

	cosfi = cos( rotstep1);
	sinfi = sin( rotstep1);
	GetRotZ( rotstep1 );
	VecRotateZ( cosfi, sinfi);	  
      }

      else if (strcmp(argv[2], "xy") == 0) VecRotTo_XY();
      else if (strcmp(argv[2], "xz") == 0) VecRotTo_XZ();
      else if (strcmp(argv[2], "yz") == 0) VecRotTo_YZ();
      else if (strcmp(argv[2], "ab") == 0 ) VecRotTo_AB();
      else if (strcmp(argv[2], "ac") == 0 ) VecRotTo_AC();
      else if (strcmp(argv[2], "bc") == 0 ) VecRotTo_BC();
      else {
	Tcl_SetResult(interp, "invalid rotation direction in xc_rotate command, must be one of +x, -x, +y, -y, +z, -z, ++xy, ++xz, ++yz, xy, xz, yz, ab, ac, bc", TCL_STATIC);
	return TCL_ERROR;
      }
      /* assign previous MajorMat to OldMat and than get new MajorMat */ 
      GetOldMat(); /* this func. assign temp. MajorMat to OldMat and assign
		      nullMat to MajorMat */
      GetMajorMat();
    }

  /* now Update display */
  Togl_PostRedisplay(togl);

  return TCL_OK;
}


int
xcRotateX(double fi)
{
  int i;
  double cosfi, sinfi;

  cosfi = cos(fi);
  sinfi = sin(fi);

  /* for X rotation order is (Y,Z,X) */
  for(i = 1; i <= tmp_nobjects; i++)
    {
      Rotate(&((coor + i)->y1), &((coor + i)->z1), &((coor + i)->x1), 
	     cosfi, sinfi);
      if ( (coor + i)->flag >= BOND ) {
	Rotate(&((coor + i)->y2), &((coor + i)->z2), &((coor + i)->x2), 
	       cosfi, sinfi);
	*(zorient + i) = 0.5 * ((coor + i)->z1 + (coor + i)->z2) + 
	  Z_OFFSET((coor + i)->flag);
      } else /* ATOM, SELATOM */
	*(zorient + i) = (coor + i)->z1 + Z_OFFSET((coor + i)->flag); 
    }
  /* rotate coor-sistem base vectors */
  VecRotateX(cosfi, sinfi);

  return XC_OK;
}


int
xcRotateY(double fi)
{
  int i;
  double cosfi, sinfi;

  cosfi = cos(fi);
  sinfi = sin(fi);

  /* for Y rotation order is (Z,X,Y) */
  for(i = 1; i <= tmp_nobjects; i++)
    {
      Rotate(&((coor + i)->z1), &((coor + i)->x1), &((coor + i)->y1), 
	     cosfi, sinfi);
      if ( (coor + i)->flag >= BOND ) {
	Rotate(&((coor + i)->z2), &((coor + i)->x2), &((coor + i)->y2), 
	       cosfi, sinfi);
	*(zorient + i) = 0.5 * ((coor + i)->z1 + (coor + i)->z2) + 
	  Z_OFFSET((coor + i)->flag);
      } else /* ATOM */
	*(zorient + i) = (coor + i)->z1 + Z_OFFSET((coor + i)->flag);
    }
  /* rotate coor-sistem base vectors */
  VecRotateY(cosfi, sinfi);

  return XC_OK;
}


int
xcRotateZ(double fi)
{
  int i;
  double cosfi, sinfi;

  cosfi = cos(fi);
  sinfi = sin(fi);

  for(i = 1; i <= tmp_nobjects; i++)
    {
      Rotate(&((coor + i)->x1), &((coor + i)->y1), &((coor + i)->z1), 
	     cosfi, sinfi);
      if ( (coor + i)->flag >= BOND ) {
	Rotate(&((coor + i)->x2), &((coor + i)->y2), &((coor + i)->z2), 
	       cosfi, sinfi);
	*(zorient + i) = 0.5 * ((coor + i)->z1 + (coor + i)->z2) + 
	  Z_OFFSET((coor + i)->flag);
      } else /* ATOM */
	*(zorient + i) = (coor + i)->z1 + Z_OFFSET((coor + i)->flag);
    }
  /* rotate coor-sistem base vectors */
  VecRotateZ(cosfi, sinfi);

  return XC_OK;
}


int
xcRotateXY(double fiX, double fiY)
{
  int i;
  double cosfiX, sinfiX;
  double cosfiY, sinfiY;

  cosfiX = cos(fiX);
  sinfiX = sin(fiX);
  cosfiY = cos(fiY);
  sinfiY = sin(fiY);

  printf("nobject=%d, tmp_nobject=%d\n",nobjects, tmp_nobjects);
  /* for X rotation order is (Y,Z,X) */
  for(i = 1; i <= tmp_nobjects; i++)
    {
      /* rotate X */
      Rotate(&((coor + i)->y1), &((coor + i)->z1), &((coor + i)->x1), 
	     cosfiX, sinfiX);
      /* rotate Y */
      Rotate(&((coor + i)->z1), &((coor + i)->x1), &((coor + i)->y1), 
	     cosfiY, sinfiY);
      if ( (coor + i)->flag >= BOND ) {
	Rotate(&((coor + i)->y2), &((coor + i)->z2), &((coor + i)->x2), 
	       cosfiX, sinfiX);
	Rotate(&((coor + i)->z2), &((coor + i)->x2), &((coor + i)->y2), 
	       cosfiY, sinfiY);
	*(zorient + i) = 0.5 * ((coor + i)->z1 + (coor + i)->z2) +
	  Z_OFFSET((coor + i)->flag);
      } else /* ATOM */
	*(zorient + i) = (coor + i)->z1 + Z_OFFSET((coor + i)->flag); 
    }
  /* rotate coor-sistem base vectors */
  VecRotateXY(cosfiX, sinfiX, cosfiY, sinfiY);
  return XC_OK;
}


/* XC_TranslateCb --> implementation of 'xc_translate' custom Togl subcommand 
 * ---------------------------------------------------------------------------
 * 
 * Usage: xc_translate <toglName> <dir> <translstep>
 *
 *        dir     --  +x | -x | +y | -y | +z | -z
 *                
 *        translstep --  translation step (in % of canvas dim.) 
 *
 * NOTE: this command is invoked by ButtonWidget-press;
 *       mouse drivven translation is now implemented by XC_MesaEvent
 *
 */
int 
XC_TranslateCb(ClientData clientData, Tcl_Interp *interp,
	       int argc, const char *argv[])
{
  Togl *togl;
  double translstep;
  double translstep1;

  /* this proc is meant only for .mesa widget */
  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( togl != mesa_togl) {
    Tcl_SetResult(interp, "xc_translate is meant only for .mesa", TCL_STATIC);
    return TCL_ERROR;
  }

  /* if structure is not opened, do noting, just return silently */
  if (!VPf.stropened) return TCL_OK;

  if (argc != 4) {
    Tcl_SetResult(interp, "Usage: xc_translate <toglName> <dir> <translstep> ", TCL_STATIC);
    return TCL_ERROR;
  }

  /* is translstep valid enough */
  if ( Tcl_GetDouble(interp, argv[3], &translstep) == TCL_ERROR ) {
    Tcl_AppendResult(interp,
		     "invalid translstep \"%s\" in %s xc_translate command",
		     argv[0], argv[2], (char *)NULL);
    return TCL_ERROR;
  }

  if ( togl == mesa_togl ) {
    /* from % to screen units */
    translstep1 = translstep * (double) VPf.canvassize;

    if (strcmp(argv[2],"+x") == 0) 
      tr.xtransl += translstep1;
    else if (strcmp(argv[2],"-x") == 0) 
      tr.xtransl -= translstep1;
    else if (strcmp(argv[2],"+y") == 0) 
      tr.ytransl += translstep1;
    else if (strcmp(argv[2],"-y") == 0) 
      tr.ytransl -= translstep1;
    else if ( strcmp(argv[2],"+z") == 0) 
      tr.zoom *= 1.0 + translstep;
    else if ( strcmp(argv[2],"-z") == 0)
      tr.zoom /= (1.0 + translstep);    
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss),"invalid translation direction \"%s\" in xc_translate command, must be one of +x, -x, +y, -y, +z, -z",argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    /* now make a ViewPort transformation */
    xcViewPort();

    /* now reset lights */
    LoadLights();

    /* now update the display */
    Togl_PostRedisplay(togl);
  } else {
    /* do nothing for the time being */
    /*
      NEW_WIN_CONTEXT *wc = FindWinContextByTogl (togl);
      wc->tr.zoom *= (1.0 + zoom);
      crySetProjection (wc, togl);      
    */
  }

  return TCL_OK; 
}


/* XC_CloseStrCmd --> inplementation of 'xc_closestr' custom Tcl command
 * ---------------------
 * Usage: xc_closestr <toglName>
 *                 
 * Effect: all pointers that belong to displayed structure are freed
 */
int 
XC_CloseStrCmd(ClientData clientData, Tcl_Interp *interp,
	       int argc, const char *argv[])
{
  struct Togl *togl;

  if (argc != 2) {
    Tcl_SetResult(interp, "Usage: xc_closestr <toglName>", TCL_STATIC);
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

  /* structure can be closed only ones */
  if (VPf.stropened == 0) {
    Tcl_SetResult(interp, "structure has already been closed or not opened yet!!!", TCL_STATIC);
    return TCL_ERROR;
  }
  VPf.stropened = 0;

  /* make a black polygon -> so structure will disapear */
  xcClearScreen(togl);

  /* Togl_PostRedisplay(togl); */
  /* delete a MesaEventHandler */
  /* Tk_DeleteEventHandler(tk_mesawin, event_mask,  
     (Tk_EventProc *)XC_MesaEvent, (ClientData) tk_mesawin);*/

  /* maybe there are some 3D lists; destroy them */
  /*xcMaybeDestroyLists();*/

  /* T.K.__CORECTURE: XCDISPLAY = (VOID *) 0 */
  /* xcDisplayFunc( (void *)NULL ); */
  /* xcDisplay = (void *) NULL; */

  xcDisplayFunc(xcDummyDisplay);

  /* destroy the glXContext; free the colormap associated with it */
  /* xcDestroyMesaContext(); */

  /* clean all stuff like DATAGRID, SUPERCELL, BZ, WIGNER-SEITZ,
     ISOSURFACE, etc.
  */
  CloseGridList();

  /* now FREE all variables */
  FreeAllVariables();

  ReloadVars();  /* reset all variables that need to be reseted */ 

  return TCL_OK;
}


extern AtomBond *coor3D;
void
FreeAllVariables(void)
{
  xcFree((void*) sqn);
  xcFree((void*) iwksp);
  xcFree((void*) nat);
  xcFree((void*) xat);
  xcFree((void*) yat);
  xcFree((void*) zat);
  xcFree((void*) fv);
  xcFree((void*) bondend);
  xcFree((void*) xbond);
  xcFree((void*) ybond);
  xcFree((void*) zbond);
  xcFree((void*) xbond2);
  xcFree((void*) ybond2);
  xcFree((void*) zbond2);
  xcFree((void*) zbmid);
  xcFree((void*) natbond);
  xcFree((void*) sqnbond);
  xcFree((void*) coor);
  xcFree((void*) zorient);
  xcFree((void*) coor3D);
  /* following variables has been malloced only if nframes > 0 */
  if ( nframes > 0 ) {
    xcFree((void*) xframe);
    xcFree((void*) yframe);
    xcFree((void*) zframe);
    xcFree((void*) xframe2);
    xcFree((void*) yframe2);
    xcFree((void*) zframe2);
    xcFree((void*) frametype);
  }
  /*for(i=0; i<N_LIST_PRIMITIVES; i++)
    if (primitive_Lists[i]) free((FREE_ARG) primitive_Lists[i] );*/
}


/* XC_DisplayMode3DCmd --> inplementation of 'xc_displayMode3D' custom 
 * Tcl command
 * ---------------------
 * Usage: xc_displayMode3D <toglName> over| boolean <displayMode>
 */
int 
XC_DisplayMode3DCmd(ClientData clientData, Tcl_Interp *interp,
		    int argc, const char *argv[])
/* =================================================================
 *     procedure assign wanted display function to xcDisplayFunc
 * =================================================================
 */
{ 
  struct Togl *togl;

  if (argc != 4) {
    Tcl_SetResult(interp, "Usage: xc_displayMode3D <toglName> over|boolean <displayMode>", TCL_STATIC);
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

  /* if structure is not opened, return silently */
  if (!VPf.stropened) return TCL_OK;

  /* if we are in selection mode, mode-changing is not allowed, return 
     silenlty */
  if (VPf.selection) return TCL_OK;

  /* argv[1] must be either 'over', either 'boolean' */
  if ( strcmp(argv[2],"over") != 0 &&
       strcmp(argv[2],"boolean") != 0 ) {
    Tcl_SetResult(interp, "Usage: xc_displayMode3D over|boolean <displayMode>", TCL_STATIC);
    return TCL_ERROR;  
  }

  /* --- Is <displayMode> OK ??? --- */
  /* S  -- sticks
   * B  -- balls
   * BS -- balls&sticks
   * SF -- spacefill 
   * PB -- pipe&balls
   */
  if ( strcmp(argv[3],"S") != 0 &&
       strcmp(argv[3],"B") != 0 &&
       strcmp(argv[3],"SF") != 0 &&
       strcmp(argv[3],"PB") != 0 &&
       strcmp(argv[3],"BS") != 0 ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"Unknown <displayMode> \"%s\", must be one of S, B, SF, BS", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* First destroy BallLists, if they exist */
  /*xcMaybeDestroyBallLists();*/

  /* if previous displayMode was one of 2D, than we must prepare 
   * for 3D "AtomBond3D coor3D"
   */
  printf("XC_DisplayMode3DCmd:: --before From2Dto3D()\ndimType=%d",dimType);
  fflush(stdout);

  /* glEnable things for 3D */  
  dimType = XC_3D;
  EnableOr2D_Or3D();

  if ( strcmp(argv[2],"boolean") == 0 ) {
    if ( strcmp(argv[3],"S") == 0 ) {
      is.stickmode = !is.stickmode;
      is.spacefillmode = GL_FALSE;
      xcMakeProjection3D("sticks");
      /* xcMakeStick3DLists(); */
    }
    else if ( strcmp(argv[3],"B") == 0 ) {
      is.ballmode = !is.ballmode;
      is.spacefillmode = GL_FALSE;
      xcMakeProjection3D("balls");
      /* xcMakeBall3DLists(); */
    }  
    else if ( strcmp(argv[3],"BS") == 0 ) {
      /* actually this is more over mode */
      is.stickmode = GL_TRUE;
      is.ballmode = GL_TRUE;
      is.spacefillmode = GL_FALSE;
      xcMakeProjection3D("balls");
      /* xcMakeStick3DLists();
	 xcMakeBall3DLists(); */
    }
    else if ( strcmp(argv[3],"SF") == 0 ) {
      is.spacefillmode = !is.spacefillmode;
      is.stickmode = GL_FALSE;
      is.ballmode = GL_FALSE;
      xcMakeProjection3D("space");
      /* xcMakeSpaceFill3DLists(); */
    }
  }
  else if ( strcmp(argv[2],"over") == 0 ) {
    if ( strcmp(argv[3],"S") == 0 ) {
      is.stickmode = GL_TRUE;
      is.ballmode = GL_FALSE;
      is.spacefillmode = GL_FALSE;
      xcMakeProjection3D("sticks");
      /* xcMakeStick3DLists(); */
    }
    else if ( strcmp(argv[3],"B") == 0 ) {
      is.ballmode = GL_TRUE;
      is.stickmode = GL_FALSE;
      is.spacefillmode = GL_FALSE;
      xcMakeProjection3D("balls");
      /* xcMakeBall3DLists(); */
    }  
    else if ( strcmp(argv[3],"PB") == 0 ) {
      is.pipemode  = GL_TRUE;
      is.stickmode = GL_TRUE;
      is.ballmode  = GL_TRUE;
      is.spacefillmode = GL_FALSE;
      xcMakeProjection3D("balls");
      /* xcMakeStick3DLists(); */
      /* xcMakeBall3DLists();  */
    }
    if ( strcmp(argv[3],"BS") == 0 ) {
      is.pipemode  = GL_FALSE;
      is.stickmode = GL_TRUE;
      is.ballmode  = GL_TRUE;
      is.spacefillmode = GL_FALSE;
      xcMakeProjection3D("balls");
      /* xcMakeStick3DLists(); */
      /* xcMakeBall3DLists();  */
    }
    if ( strcmp(argv[3],"SF") == 0 ) {
      is.spacefillmode = GL_TRUE;
      is.stickmode = GL_FALSE;
      is.ballmode = GL_FALSE;
      xcMakeProjection3D("space");
      /* xcMakeSpaceFill3DLists(); */
    }
  }

  /* LOAD LIGHTS & MATERIAL PROPERTYES & FRONTFACE & BLENDFUNC */
  LoadLights();
  LoadStructMaterial(); 
  LoadBlendfunc_And_Frontface();

  xcDisplayFunc(xcDisplay3D);

  Togl_PostRedisplay(togl);

  return TCL_OK;
}


/* XC_DrawStyle3DCmd --> inplementation of 'xc_drawStyle3D' custom 
 * Tcl command
 * ---------------------
 * Usage: xc_drawStyle3D <toglName> <drawStyle3D>
 */
int 
XC_DrawStyle3DCmd(ClientData clientData, Tcl_Interp *interp,
		  int argc, const char *argv[])
{
  struct Togl *togl;
  /* =================================================================
   *     procedure assign wanted drawStyle and redraw in that style
   * =================================================================
   */

  if (argc != 3) {
    Tcl_SetResult(interp, "Usage: xc_drawStyle3D <toglName> <drawStyle>", TCL_STATIC);
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

  /* if structure is not opened, return silently */
  if (!VPf.stropened) return TCL_OK;

  /* check if we are in 3D mode, if not return silently */
  if (dimType != XC_3D) return TCL_OK;

  /*EV*/
  if ( strcmp(argv[2],"solid") == 0 ) {
    is.solid    = GL_TRUE;
  }
  else if ( strcmp(argv[2],"wire") == 0 ) {
    is.solid    = GL_FALSE;
  }
  else if ( strcmp(argv[2],"anaglyph") == 0 ) {
    if (is.anaglyph == GL_TRUE) {
      is.anaglyph = GL_FALSE;
    } else {
      is.anaglyph = GL_TRUE;
    }
  }
  else if ( strcmp(argv[2],"stereo") == 0 ) {
    /* toogle the stereo flag */
    if (is.stereo == GL_TRUE) {
      /*fprintf(stderr,"Stereo = off\n");*/
      is.stereo = GL_FALSE;
      glDrawBuffer(GL_BACK);
    } else {
      /*fprintf(stderr,"Stereo = whats to be on\n");*/
      GLboolean can_do_stereo;
      glGetBooleanv(GL_STEREO, &can_do_stereo);
      if (can_do_stereo == GL_TRUE) is.stereo = GL_TRUE;
      else { 
	is.stereo = GL_FALSE;
	/*fprintf(stderr,"Stereo = off\n");*/
      }
    }
  } else {      
    char rss[1024];
    snprintf(rss, sizeof(rss),
	     "Unknown <drawStyle> \"%s\", must be one of \"solid\", \"wire\", \"anaglyph\", or \"stereo-\"", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  Togl_PostRedisplay(togl);

  return TCL_OK;
}


/* XC_ShadeModel3DCmd --> implementation of 'xc_shadeModel3D' custom 
 * Tcl command
 * ---------------------
 * Usage: xc_shadeModel3D <toglName> <shadeModel>
 */
int 
XC_ShadeModel3DCmd(ClientData clientData, Tcl_Interp *interp,
		   int argc, const char *argv[])
{
  struct Togl *togl;
  /* =================================================================
   *     procedure assign wanted drawStyle and redraw in that style
   * =================================================================
   */

  if (argc != 3) {
    Tcl_SetResult(interp, "Usage: xc_shadeModel3D <toglName> <shadeModel>", TCL_STATIC);
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

  /* if structure is not opened, return silently */
  if (!VPf.stropened) return TCL_OK;

  /* check if we are in 3D mode, if not return silently */
  if (dimType != XC_3D) return TCL_OK;

  if ( strcmp(argv[2],"smooth") == 0 ) {
    is.smooth = GL_TRUE;
    glShadeModel(GL_SMOOTH);
  }
  else if ( strcmp(argv[2],"flat") == 0 ) {
    is.smooth = GL_FALSE;
    glShadeModel(GL_FLAT);
  }
  else {       
    char rss[1024];
    snprintf(rss, sizeof(rss),
	     "Unknown <shadeModel> \"%s\", must be one \"smooth\", \"flat\"", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  Togl_PostRedisplay(togl);

  return TCL_OK;
}


/* XC_PointSizeCmd --> implementation of 'xc_pointSize' custom 
 * Tcl command
 * ---------------------
 * Usage: xc_pointSize <toglName> <pointSize>
 */
int 
XC_PointSizeCmd(ClientData clientData, Tcl_Interp *interp,
		int argc, const char *argv[])
{
  struct Togl *togl;
  int pointsize;

  if (argc != 3) {
    Tcl_SetResult(interp, "Usage: xc_pointSize <toglName> <pointSize>", TCL_STATIC);
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

  if ( Tcl_GetInt(interp, argv[1], &pointsize) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" in xc_poinsize command", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  printf("PointSize:: %d\n",pointsize);
  glPointSize( (GLfloat) pointsize );

  /* if structure is not opened, return silently */
  if (!VPf.stropened) return TCL_OK;

  Togl_PostRedisplay(togl);

  return TCL_OK;
}


/* XC_ResetVarCmd --> implementation of xc_resetvar custom Tcl command
 *
 * Usage: xc_resetvar <toglName> <var>
 *
 *                     <var> -- number to identificate a variable;
 *                              look in the struct.h -- defines with
 *                              R_ prefix; 
 * Side Effects: wanted variable is reseted to default value, and
 *               ReDisplay is made if and only if structure is opened
 */
int 
XC_ResetVarCmd(ClientData clientData, Tcl_Interp *interp,
	       int argc, const char *argv[])
{
  struct Togl *togl;
  int var;

  if (argc != 3) {
    Tcl_SetResult(interp, "Usage: xc_resetvar <toglName> <var>", TCL_STATIC);
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

  if ( Tcl_GetInt(interp, argv[2], &var) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" in xc_resetvar command", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  ResetVar(togl, var);

  return TCL_OK;
}


/* XC_OldAtmColCmd --> load old atomic colors
 *
 * Usage: xc_oldatmcol
 */
int 
XC_OldAtmColCmd(ClientData clientData, Tcl_Interp *interp,
                int argc, const char *argv[])
{
  LoadOldAtomicColors();
  return TCL_OK;
}

/* XC_NewValueCmd --> implementation of xc_newvalue custom Tcl command
 *
 * Usage: xc_newvalue <toglName> <var> ?value1? ?value2? ?value3? ?value4? 
 *
 *                     <var> -- number to identificate a variable;
 *                              look in the struct.h -- defines with
 *                              L_ prefix; 
 * Side Effects: wanted variable is set to new value, and
 *               ReDisplay is made if and only if structure is opened
 */
int 
XC_NewValueCmd(ClientData clientData, Tcl_Interp *interp,
	       int argc, const char *argv[])
{
  struct Togl *togl;
  int var, i;
  double value[4] = { 0.0, 0.0, 0.0, 0.0 };

  if (argc < 3 || argc > 7) {
    Tcl_SetResult(interp, "Usage: xc_newvalue <toglName> <var> ?value1? ?value2? ?value3? ?value4?", TCL_STATIC);
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

  if ( Tcl_GetInt(interp, argv[2], &var) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer for <var>, but got \"%s\" in xc_newvalue command", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  for (i=3; i<argc; i++)  
    if ( Tcl_GetDouble(interp, argv[i], &value[i - 3]) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted number, but got \"%s\" while executing %s %s %s", argv[i], argv[0], argv[1], argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }  


  if ( var == SET_ATOMLABEL_DO_DISPLAY ) {
    /* usage: xc_newvalue togl SET_ATOMLABEL_DO_DISPLAY atomID 0|1}*/
    int atomID = (int) value[0];
    int toggle = (int) value[1];

    if ( argc != 5 ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "wrong number of arguments %d, must be xc_newvalue togl SET_ATOMLABEL_DO_DISPLAY atomID 0|1", argc);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    if (atomID < 1 || atomID > natoms) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "atomID out of range, should be between [0,%d]", natoms);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    if (toggle) atomLabel[atomID].do_display = XC_YES;
    else atomLabel[atomID].do_display = XC_NO;    

    Togl_PostRedisplay(togl);
  } 

  else if ( var == SET_GLOBALATOMLABEL_DO_DISPLAY ) {
    /* usage: xc_newvalue togl SET_GLOBALATOMLABEL_DO_DISPLAY 0|1 */
    int toggle = (int) value[0];

    if ( argc != 4 ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "wrong number of arguments %d, must be xc_newvalue togl SET_GLOBALATOMLABEL_DO_DISPLAY 0|1", argc);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    if (toggle) globalAtomLabel.do_display = XC_YES;
    else globalAtomLabel.do_display = XC_NO;    

    Togl_PostRedisplay(togl);
  } 

  else if ( var == SET_DO_NOT_DISPLAY_ATOMLABEL ) {
    /* usage: xc_newvalue togl SET_DO_NOT_DISPLAYATOM_LABEL atomID 0|1 */
    int atomID = (int) value[0];
    int toggle = (int) value[1];

    if ( argc != 5 ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "wrong number of arguments %d, must be xc_newvalue togl SET_DO_NOT_DISPLAY_LABEL atomID 0|1", argc);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    if (atomID < 1 || atomID > natoms) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "atomID out of range, should be between [0,%d]", natoms);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    if (toggle) do_not_display_atomlabel[atomID] = XC_YES;
    else do_not_display_atomlabel[atomID] = XC_NO;    

    Togl_PostRedisplay(togl);
  }
  else {
    LoadNewValue(togl, var, value[0], value[1], value[2], value[3]);
  }

  return TCL_OK;
}


/* XC_GetDefault --> implementation of xc_getdefault custom Tcl command
 *
 * Usage: xc_getdefault <var> ?value?
 *
 *                     <var> -- number to identificate a variable;
 *                              look in the struct.h -- defines with
 *                              D_ prefix; 
 * Side Effects: default value of "<var> variable" is returned
 */
int 
XC_GetDefaultCmd(ClientData clientData, Tcl_Interp *interp,
		 int argc, const char *argv[])
{
  int var, value = 0;
  char *result = (char *) Tcl_Alloc( sizeof(char) * 64); /* maximum lenght of result string is 64 characters */
  if (argc < 2) {
    Tcl_SetResult(interp, "Usage: xc_getdefault <var> ?value?", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( Tcl_GetInt(interp, argv[1], &var) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" in xc_getdefault command", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc == 3 )
    if ( Tcl_GetInt(interp, argv[2], &value) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" in xc_getdefault command", argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

  /* ============================================================ */
  if ( var == D_ATCOL_ONE ) 
    /* case D_ATCOL */
    sprintf(result,"%f %f %f",DefAtCol[value][0],DefAtCol[value][1],
	    DefAtCol[value][2]);
  else if ( var == D_ATRAD_ONE )
    /* case D_ATRAD_ONE */
    sprintf(result,"%f", DEF_ATRADF * rcov[value]);
  else if ( var == D_RCOV_ONE)
    /* case D_RCOV_ONE */
    sprintf(result,"%f", rcov[value]);    
  else if ( var == D_FRAMECOL )
    /* case D_FRAMECOL */
    sprintf(result,"%f %f %f", DefFrameCol[0], DefFrameCol[1], DefFrameCol[2]);
  else  if ( var == D_UNIBONDCOLOR ) 
    /* case D_UNIBONDCOLOR */
    sprintf(result,"%f %f %f",
	    DefUnibondCol[0], DefUnibondCol[1], DefUnibondCol[2]);

  else if ( var == GET_FOG_COLOR ) {
    sprintf(result, "%f %f %f %f", 
	    def_fog_color[0], def_fog_color[1], def_fog_color[2], def_fog_color[3]);
  }

  else if ( var == D_FORCE_COLOR ) {
    sprintf(result, "%f %f %f %f", 
	    DefAtCol[0][0],  DefAtCol[0][1],  DefAtCol[0][2], 1.0);    
  }

  else if ( var == D_XYZ_AXIS_COLOR ) {
    sprintf(result, "%f %f %f %f",
	    def_xyz_axis_color[0],
	    def_xyz_axis_color[1],
	    def_xyz_axis_color[2],
	    def_xyz_axis_color[3]);
  }

  else if ( var == D_XYZ_XYPLANE_COLOR ) {
    sprintf(result, "%f %f %f %f",
	    def_xyz_xyplane_color[0],
	    def_xyz_xyplane_color[1],
	    def_xyz_xyplane_color[2],
	    def_xyz_xyplane_color[3]);
  }

  else sprintf(result, "%f", GetDefault(var));

  Tcl_SetResult(interp, result, TCL_DYNAMIC);

  return TCL_OK;
}


/* XC_GetValue --> implementation of xc_getvalue custom Tcl command
 *
 * Usage: xc_getvalue <var> ?value?
 *
 *                     <var> -- number to identificate a variable;
 *                              look in the struct.h -- defines with
 *                              D_ prefix; 
 * Side Effects: temporary value of "<var> variable" is returned
 */
int 
XC_GetValueCmd(ClientData clientData, Tcl_Interp *interp,
	       int argc, const char *argv[])
{
  int var, value = 0;
  char *result = (char *) Tcl_Alloc( sizeof(char) * 256); /* maximum lenght of result  string is 256 characters */
  if (argc < 2) {
    Tcl_SetResult(interp, "Usage: xc_getvalue <var> ?value?", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( Tcl_GetInt(interp, argv[1], &var) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" in xc_getvalue command", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc == 3 )
    if ( Tcl_GetInt(interp, argv[2], &value) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" in xc_getvalue command", argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

  /* ======================== */
  if ( var == D_ATCOL_ONE ) {
    /* case D_ATCOL */
    sprintf(result,"%f %f %f",atcol[value][0],atcol[value][1],
	    atcol[value][2]);
    printf("GETVALUE> color = %s\n",result);
    fflush(stdout);
  } 
  /* case D_ATRAD_ONE */
  else if ( var == D_ATRAD_ONE ) { 
    sprintf(result, "%f", atrad[value]);
  } 
  /* case D_RCOV_ONE */
  else if ( var == D_RCOV_ONE ) { 
    sprintf(result, "%f", rcov[value]);
  } 
  /* case D_FRAMECOL */
  else if ( var == D_FRAMECOL ) {  
    sprintf(result,"%f %f %f", framecol[0], framecol[1], framecol[2]);
  }
  /* case D_MAXSTRUCTSIZE */
  else if ( var == D_MAXSTRUCTSIZE ) {
    sprintf(result,"%f %f %f", max.x, max.y, max.z);
  }
  else if ( var == D_UNIBONDCOLOR ) {
    /* case D_UNIBONDCOLOR */
    sprintf(result,"%f %f %f", unibondCol[0], unibondCol[1], unibondCol[2]);
  }
  else if ( var == D_BACKGROUND ) {
    /* case D_UNIBONDCOLOR */
    sprintf(result,"%f %f %f %f", bg[0], bg[1], bg[2], bg[3]);
  }
  else if ( var == GET_NATOMS ) {
    if ( ! VPf.stropened) {
      sprintf(result, "-1");
    } else {
      sprintf(result, "%d", natoms);
    }
  }
  else if ( var == GET_NAT ) {
    if ( ! VPf.stropened) {
      sprintf(result, "-1");
    } else {
      if ( value < 1 || value > natoms ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "atom index %d too large, must be between [1,%d]",
		 value, natoms);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    sprintf(result, "%d", *(nat + value));  
  }
  else if ( var == GET_ATOMLABEL_LABEL ) {
    if ( value < 1 || value > natoms ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "atom index %d too large, must be between [1,%d]",
	       value, natoms);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( atomLabel[value].label ) {
      sprintf(result, "%s", atomLabel[value].label);
    } else {
      result[0] = '\0';
    }
  }
  else if ( var == GET_ATOMLABEL_BRIGHTCOLOR ) {
    if ( value < 1 || value > natoms ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "atom index %d too large, must be between [1,%d]",
	       value, natoms);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( atomLabel[value].base ) {
      sprintf(result, "%f %f %f",
	      atomLabel[value].bright_color[0], 
	      atomLabel[value].bright_color[1], 
	      atomLabel[value].bright_color[2]);
    } else {
      sprintf(result, "1.0 1.0 1.0");
    }
  }
  else if ( var == GET_ATOMLABEL_DARKCOLOR ) {
    if ( value < 1 || value > natoms ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "atom index %d too large, must be between [1,%d]",
	       value, natoms);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( atomLabel[value].base ) {
      sprintf(result, "%f %f %f",
	      atomLabel[value].dark_color[0], 
	      atomLabel[value].dark_color[1], 
	      atomLabel[value].dark_color[2]);
    } else {
      sprintf(result, "0.0 0.0 0.0");
    }
  }

  else if ( var == GET_ATOMLABEL_DO_DISPLAY ) {
    if (value < 1 || value > natoms) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "atom index %d out of range, should be between [0,%d]", 
	       value, natoms);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    sprintf(result, "%d", atomLabel[value].do_display);    
  } 

  else if ( var == GET_ATOMLABEL_ALL_ID ) {
    /* return all the ID's of all atoms than have cutom-labels */
    int i;
    for (i=1; i<=natoms; i++) {
      if (atomLabel[i].base) {
	sprintf(result, "%d ", i);
	Tcl_AppendElement(interp, result);
      }
    }
    Tcl_Free(result);
    return TCL_OK;
  }

  else if ( var == GET_GLOBALATOMLABEL_DO_DISPLAY ) {
    sprintf(result, "%d", globalAtomLabel.do_display);
  }

  else if ( var == GET_FOG_COLOR ) {
    sprintf(result, "%f %f %f %f", 
	    fog.color[0], fog.color[1], fog.color[2], fog.color[3]);
  }

  else if ( var == D_FORCE_COLOR ) {
    sprintf(result, "%f %f %f %f", 
	    FV.color[0],  FV.color[1],  FV.color[2], 1.0);    
  }

  else if ( var == D_XYZ_AXIS_COLOR ) {
    sprintf(result, "%f %f %f %f",
	    xyz.axis_color[0],
	    xyz.axis_color[1],
	    xyz.axis_color[2],
	    xyz.axis_color[3]);
  }

  else if ( var == D_XYZ_XYPLANE_COLOR ) {
    sprintf(result, "%f %f %f %f",
	    xyz.xyplane_color[0],
	    xyz.xyplane_color[1],
	    xyz.xyplane_color[2],
	    xyz.xyplane_color[3]);
  }

  else {
    /* ELSE get the result from GetValue(var) function */
    sprintf(result, "%.15f", GetValue(var));
  }
  Tcl_SetResult(interp, result, TCL_DYNAMIC);

  return TCL_OK;
}


/* XC_UpdateStrCmd --> inplementation of 'xc_updatestr' custom Tcl command
 * ---------------------
 * Usage: xc_updatestr <file> <toglName> ?format?
 *
 *                 toglName --- name of togl where structure will be 
 *                              displayed
 *                 format   --- format of the file ( xsf | pdb ); default = xsf
 *
 * Side effects: when changing "number of cells drawn" for periodic structure
 *               we want situation to be as before, just number of cells is 
 *               supposed to change, so we must preserve a lot of stuff
 */
int 
XC_UpdateStrCmd(ClientData clientData, Tcl_Interp *interp,
		int argc, const char *argv[])
{
  struct Togl *togl;
  FILE *fp;
  int format = FORMAT_XSF;

  fprintf(stderr, "*** XC_UpdateStrCmd; argc == %d\n", argc);
  fflush(stderr);
  
  if (argc < 3 || argc > 4) {
    Tcl_SetResult(interp, "Usage: xc_updatestr <file> <toglName> ?format?", TCL_STATIC);
    return TCL_ERROR;
  }
  
  /* --- if structure is not opened just return silently --- */
  if (!VPf.stropened) return TCL_OK;

  /* Check if 'file' is OK */
  fp = fopen(argv[1],"r");
  if (fp == NULL) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "Can't open file \"%s\"",argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* find togl associated with toglName */
  if ( Togl_GetToglFromName(interp, argv[2], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( argc == 4 ) {
    /* determine file format */
    if ( (strcmp(argv[3],"xcr") == 0) || (strcmp(argv[3],"xsf") == 0) ) format = FORMAT_XSF;
    else if ( strcmp(argv[3],"xyz") == 0 ) format = FORMAT_XYZ;
    else if ( strcmp(argv[3],"pdb") == 0 ) format = FORMAT_PDB;
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss),"unknown format \"%s\", must be one of \"xsf\", \"xyz\", \"pdb\"", argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  /* free neccesarry variables */
  FreeAllVariables();

  /* reset all variables needed */
  natoms = 0;
  nbonds = 0;
  nframes = 0;
  tmp_nobjects = 0;
  nobjects = 0;

  /* command is OK 
   * LET'S READ A FILE
   */  
  if (!ReadStructFile(fp, argv[1], format, XSF_UPDATE)) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"Error reading file \"%s\"",argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    fclose(fp);
    VPf.stropened = 0;
    xcClearScreen(togl);
    /*xcMaybeDestroyLists();*/
    xcDisplayFunc(xcDummyDisplay);
    ReloadVars();
    Togl_PostRedisplay(togl);  
    return TCL_ERROR;
  }

  /* now we can close the file */
  fclose(fp);

  /* Update DisplayFunc */
  UpdateDispFunc();

  /* update a projection */
  UpdateProjection();

  xcViewPort();

  /* assing the sInfo global array */
  Set_sInfoArray( interp );

  /* Display a Structure */
  Togl_PostRedisplay(togl);

  /* if we come so far we were succesfully */
  return TCL_OK;
}


/* 
   Revision: Thu Nov 12 08:58:59 CET 1998

   this routine was used in the times of mesaWin. At present we use Togl and
   the routine is maintained just for backward compatibility; it was adapted
   for togl
*/
/*****************************************************************************
 * sometimes due to some strange reconfiguration of <mesaWin> the displayed  *
 * structure vanishes; this command is provided to recover from that kind    *
 * of event                                                                  *
 *****************************************************************************/
int 
XC_SwapBufferCmd(ClientData clientData, Tcl_Interp *interp,
		 int argc, const char *argv[])
{
  struct Togl *togl;
  if (argc != 2) {
    Tcl_SetResult(interp, "Usage: xc_swapbuffer <toglName>", TCL_STATIC);
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

  if ( VPf.stropened )
    Togl_SwapBuffers( togl );

  return TCL_OK;
}


/* 
   Revision: Thu Nov 12 08:59:56 CET 1998

   this routine was used in the times of mesaWin. At present we use Togl and
   the routine is maintained just for backward compatibility; it was adapted
   for togl
*/
/*****************************************************************************
 * sometimes we just want to run the Display Function                        *
 *****************************************************************************/
int 
XC_DisplayCmd(ClientData clientData, Tcl_Interp *interp,
	      int argc, const char *argv[])
{
  struct Togl *togl; 

  if (argc != 2) {
    Tcl_SetResult(interp, "Usage: xc_display <toglName>", TCL_STATIC);
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

  if ( VPf.stropened ) Togl_PostRedisplay(togl);

  return TCL_OK;
}



/*
 *
 * Usage: xc_rotationmatrix get
 *     or xc_rotationmatrix set v1 v2 v3 v4 v5 v6 ....
 *
 * the appropriate vector is double vec.crdvec[16] !!!
 */
int XC_RotationMatrixCmd(ClientData clientData, Tcl_Interp *interp,
			 int argc, const char *argv[])
{
  int i;

  if ( strcmp(argv[1], "get") == 0 ) 
    {
      char *string = calloc(30, sizeof(char));

      for (i=0; i<16; i++) {
	sprintf(string, "%22.15e ",vec.crdvec[i]);
	Tcl_AppendResult(interp, string, (char*)NULL);
      }

      free( (FREE_ARG) string);
    }
  else if ( strcmp(argv[1], "set") == 0 ) 
    {
      double number;
      if ( argc != 18 ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "wrong # args: should be \"xc_rotationmatrix set v1 ... v16\"");
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      for (i=0; i<16; i++) 
	{
	  if ( Tcl_GetDouble(interp, argv[i+2], &number) == TCL_ERROR ) {
	    char rss[1024];
	    snprintf(rss, sizeof(rss), "wanted double but got %s", argv[i+2]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	  vec.crdvec[i] = number;
	}
      vecVecToMat ( vec.crdvec, vec.crdmajor );      

      if ( dimType == XC_2D ) {
	/* recalculate the 2D-coor accoring to new rotation matrix */
	/* this call will do the job */
	xcAssignDisplayFunc(displayMode2D);
      }
    }
  else 
    {
      char rss[1024];
      snprintf(rss, sizeof(rss), "unknown mode %s, must be get or set", argv[1]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

  return TCL_OK;
}


/*
 *
 * Usage: xc_translparam get
 *     or xc_translparam set xtransl ytransl zoom
 */
int XC_TranslParamCmd(ClientData clientData, Tcl_Interp *interp,
		      int argc, const char *argv[])
{
  int i;

  if ( strcmp(argv[1], "get") == 0 ) 
    {
      char *string = calloc(30, sizeof(char));           
      sprintf(string, "%22.15e ", tr.xtransl); Tcl_AppendResult(interp, string, (char*)NULL);
      sprintf(string, "%22.15e ", tr.ytransl); Tcl_AppendResult(interp, string, (char*)NULL);
      sprintf(string, "%22.15e ", tr.zoom);    Tcl_AppendResult(interp, string, (char*)NULL);
      free( (FREE_ARG) string);
    }
  else if ( strcmp(argv[1], "set") == 0 ) 
    {
      double number[3];      
      if ( argc != 5 ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "wrong # args: should be \"translparam set xtransl ytransl zoom\"");
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      for (i=0; i<3; i++) 
	{
	  if ( Tcl_GetDouble(interp, argv[i+2], number + i) == TCL_ERROR ) {
	    char rss[1024];
	    snprintf(rss, sizeof(rss), "wanted double but got %s", argv[i+2]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	}
      tr.xtransl = number[0];
      tr.ytransl = number[1];
      tr.zoom    = number[2];
      UpdateProjection();
      xcViewPort();
    }
  else 
    {
      char rss[1024];
      snprintf(rss, sizeof(rss), "unknown mode %s, must be get or set", argv[1]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

  return TCL_OK;
}


/*
  Checks if the stereo is supported by a hardware

  Note: probably this should be called after the togl creation

  Usage: xc_stereo
*/
int 
XC_StereoCmd(ClientData clientData,Tcl_Interp *interp, int argc, const char *argv[])
{
  char *result = (char *) Tcl_Alloc( sizeof(char) * 128);
  GLboolean stereo;

#if defined(TOGL_X11)
  Display *display;
  int screen;
  XVisualInfo *info;
  int list1[] = { GLX_RGBA, GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1,
		  GLX_BLUE_SIZE, 1, GLX_STEREO, (int) None, };

  display = XOpenDisplay(NULL);
  screen = DefaultScreen(display);
  info = glXChooseVisual(display, screen, list1);
  if (info) {
    stereo = GL_TRUE;
    XFree(info);
  }
  XCloseDisplay(display);
#else
  glGetBooleanv(GL_STEREO, &stereo);
#endif

  if (stereo == GL_TRUE) 
    {
      fprintf(stderr, "\n*** the hardware supports the stereo ***\n\n");
      sprintf(result, "%d", 1);
    } 
  else 
    {
      fprintf(stderr, "\n*** the hardware does not support the stereo ***\n\n");
      sprintf(result, "%d", 0);    
    }

  Tcl_SetResult(interp, result, TCL_DYNAMIC); 
  return TCL_OK;
}

