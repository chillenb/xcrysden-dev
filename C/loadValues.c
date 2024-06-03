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
 * Source: $XCRYSDEN_TOPDIR/C/loadValues.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <togl.h>
#include <stdio.h>
#include <math.h>
#include "struct.h" 
#include "atoms.h"
#include "vector.h"
#include "xcfunc.h"
#include "memory.h"

extern struct Togl *mesa_togl;
extern Options3D   is;
extern OrthoProj   ort;
extern StructSize  ss;
extern ForceVector FV;
extern RenderVectors *forceVectors;
extern H_Bond      hbonds;

XYZ_Attrib xyz;
GLfloat def_xyz_axis_color[4]    = { 1, 1, 1, 1 };
GLfloat def_xyz_xyplane_color[4] = { 0, 0, 1, 1 };

GLfloat def_fog_colormode   = XC_FOG_CUSTOMCOLOR;
GLfloat def_fog_color[4]    = {0.0, 0.0, 0.0, 1.0};
GLint   def_fog_mode        = XC_FOG_LINEAR;
GLfloat def_fog_density     = 0.1;
GLfloat def_fog_ort_start_f = 0.8;
GLfloat def_fog_ort_end_f   = 1.8;
GLfloat def_fog_persp_f1    = 0.9;
GLfloat def_fog_persp_f2    = 1.0;
XCfog fog;

int   def_antialias_degree = 1;
float def_antialias_offset = 0.39;
XCantialias antialias;

extern GLfloat DefAtCol_Ambient_by_Diffuse;
GLfloat AtCol_Ambient_by_Diffuse;

/* FUNCTION PROTOTYPES */
void LoadDefaultValues(void);
static void DefValues(void);
void ResetVar(struct Togl *togl, int var);
void ReloadVars(void);
void LoadNewValue(struct Togl *togl, int var, double value1, double value2, 
		  double value3, double value4);
double GetDefault(int var);
double GetValue(int var);
void LoadOldAtomicColors(void);

/* --- xcDisplayFunc.c --- */
/* extern void (*xcDisplay)(void); */
void xcMakePointList(void);
/*  extern void xcMakeFrame3DLists(void);  */
/*  extern void xcMakeBallLabel3DList(void); */
/*  extern void xcMakeSpaceLabel3DList(void); */
extern void xcTurnXYZOff(struct Togl *togl);
extern void xcClearScreen(struct Togl *togl);
extern void xcChangeBackground(struct Togl *togl, GLclampf bckg[4]);

/* --- xcDisplayFunc2.c --- */
extern void xcTurnLabelsOn(void);
extern void xcTurnFramesOn(void);
extern void xcTurnFramesOff(void);

/* --- remakestr.c --- */
extern void ReDisplay(struct Togl *togl, GLuint remake3D);

/* --- xcColorScheme.c --- */
extern void LoadAtmCol(const int flag);


/* --- cryWinContext.c --- */
extern NEW_WIN_CONTEXT *FindWinContextByTogl(struct Togl *togl);

extern float xcLinf(float value);

void
LoadDefaultValues(void)
{
  ReloadVars();
  DefValues();
}


static void
DefValues(void)
{
  int i, j;
  
  /* this must be changed first */
  for(i=1; i <= MAXNAT; i++)
    rcov[i] = DEF_RCOVF * rcovdef[i];      

  /* VPf.scrf        = DEF_SCRF; */
  VPf.ballf         = DEF_BALLF;
  VPf.rodf          = DEF_RODF;
  VPf.rcovf         = DEF_RCOVF;
  VPf.WFlinewidth   = DEF_LineWidthWF;
  VPf.WF3Dlinewidth = DEF_LineWidthWF3D;
  VPf.OUTlinewidth  = DEF_LineWidthOUT;
  VPf.PLlinewidth   = DEF_LineWidthPL;
  VPf.PLradius      = DEF_RadiusPL;
  rrod              = DEF_ATRADF * DEF_BALLF * DEF_RODF * rcovdef[1];
  VPf.atradf        = DEF_ATRADF;
  VPf.tessFactor    = DEF_TESSELLATION;
  VPf.xyzOn         = 0;
  VPf.labelsOn      = 0;
  VPf.framesOn      = 0;
  VPf.framewidth    = DEF_FRAMEWIDTH;
  VPf.framef        = DEF_FRAMERODF;
  VPf.selection     = GL_FALSE;
  VPf.atomadd       = GL_FALSE;
  VPf.projmade      = GL_FALSE;
  VPf.wignerseitz   = GL_FALSE;
  VPf.supercell     = GL_FALSE;
  for(i=ISOOBJ_BASE; i<MAX_ISOOBJECTS; i++) {
    VPf.isosurf[i]     = GL_FALSE;
    VPf.colorplane[i]  = GL_FALSE;
    VPf.isoline[i]     = GL_FALSE;
  }
  VPf.isospacesel2D = GL_FALSE;
  VPf.isospacesel3D = GL_FALSE;
  VPf.force         = GL_FALSE;
  VPf.Hbond         = DEF_HBOND;

  VPf.unibond          = DEF_UNIBOND;
  VPf.perspective      = DEF_PERSPECTIVE;  
  VPf.perspective_size = DEF_PERSPECTIVEFRONT;
  VPf.perspective_fovy = DEF_PERSPECTIVEFOVY;
  VPf.perspective_far  = DEF_PERSPECTIVEBACK;
  
  VPf.fog       = DEF_FOG;
  VPf.antialias = DEF_ANTIALIAS;
  
  /* new isosurfaces */
  VPf.nsurface = 0;

  rframe          = DEF_FRAMERODF * rcovdef[1];
  /* default for frames fill be LineFrame */
  is.lineframe    = GL_TRUE;

  /* ALL THAT'S HERE IS DEFAULT */
  for(i=0; i <= MAXNAT; i++)
    {
      atrad[i] = DEF_ATRADF * rcovdef[i];
      rball[i] = DEF_BALLF * atrad[i];
      for(j=0; j < 3; j++)
	atcol[i][j] = DefAtCol[i][j];
      LoadAtmCol(0);
    }
  FindMaxRad();

  for(j=0; j < 3; j++)
	framecol[j] = DefFrameCol[j];

  /* to prevent VP.VPfactor being Nan at begining in xcvieport.c */
  ort.size = 1.0;

  for(i=0; i<4; i++) {
    bg[i] = DefBg[i];
    unibondCol[i] = DefUnibondCol[i];
  }

  /* DEFAULTS for forceVectors*/
  FV.lengthfactor = XC_FORCE_LENGTHFACTOR;
  FV.threshold    = XC_FORCE_TRESHHOLD;
  FV.ScaleFunc    = xcLinf; 
  FV.rod_thickf   = VECTOR_THICKF;
  FV.arr_thickf   = VECTOR_ARRTHICKF;
  FV.arr_lenf     = VECTOR_ARROWSIZE;
  COPY_V(3, FV.color, DefAtCol[0]);   
  FV.color[3] = 1.0;

  /* DEFAULTS for FOG (i.e. depth cuing) */

  COPY_V(4, fog.color, def_fog_color);
  fog.mode        = def_fog_mode;
  fog.density     = def_fog_density;
  fog.colormode   = def_fog_colormode;
  fog.ort_start_f = def_fog_ort_start_f;
  fog.ort_end_f   = def_fog_ort_end_f;
  fog.persp_f1    = def_fog_persp_f1;
  fog.persp_f2    = def_fog_persp_f2;

  antialias.degree = def_antialias_degree;
  antialias.offset = def_antialias_offset;

  AtCol_Ambient_by_Diffuse = DefAtCol_Ambient_by_Diffuse;

  current_file_format = FORMAT_NONE;

  /* H-bonds */
  if ( !hbonds.H_like_list ) {
    hbonds.H_like_list = (int*) xcCalloc((size_t)2, sizeof(int));
  }
  if ( !hbonds.O_like_list ) {
    hbonds.O_like_list = (int*) xcCalloc((size_t)6, sizeof(int));
  }
  hbonds.H_like_list[0] = 1;  /* Hydrogen    */
  hbonds.H_like_list[1] = 0;  /* end-of-list */
  /**/
  hbonds.O_like_list[0] = 7;  /* Nitrogen    */
  hbonds.O_like_list[1] = 8;  /* Oxygen      */
  hbonds.O_like_list[2] = 9;  /* Fluorine    */
  hbonds.O_like_list[3] = 16; /* Sulfur     */
  hbonds.O_like_list[4] = 17; /* Chlorine    */
  hbonds.O_like_list[5] = 0;  /* end-of-list */

  COPY_V(4, xyz.axis_color,    def_xyz_axis_color);
  COPY_V(4, xyz.xyplane_color, def_xyz_xyplane_color);

  /* turn off stereo mode and anaglyphs by default */
  is.stereo   = GL_FALSE;
  is.anaglyph = GL_FALSE;
}


/*****************************************************************************
 * this is defined in "struct.h"
 * ----------------------------
 * define R_ATRAD           0
 * define R_RBALL           1
 * define R_RROD            2
 * define R_ATCOL           3
 * define R_WFLINEWIDTH     4
 * define R_PLLINEWIDTH     5
 * define R_PLRADIUS        6
 * define R_SCRF            7
 * define R_ALL             8
 * define R_FRAMECOL        9
 * define R_FRAMELINEWIDTH 10
 * define R_FRAMERODF      11
 *------------------------------
 */

void
ResetVar(struct Togl *togl, int var)
{
  int i,j;
  GLuint remake3D = REMAKE3D_NONE;
    
  switch(var)
    {
    case R_ATRAD:
      for(i=0; i <= MAXNAT; i++) {
	atrad[i] = DEF_ATRADF * rcovdef[i];
	rball[i] = VPf.ballf * DEF_ATRADF * rcovdef[i];
      }
      FindMaxRad();
      rrod = VPf.ballf * VPf.rodf * DEF_ATRADF * rcovdef[1];
      VPf.atradf = DEF_ATRADF;
      remake3D = REMAKE3D_BALL | REMAKE3D_SPACE | REMAKE3D_BOND |
	REMAKE3D_BALL_LABEL | REMAKE3D_SPACE_LABEL;
      break;

    case R_RCOV:
      for(i=1; i <= MAXNAT; i++)
	rcov[i] = DEF_RCOVF * rcovdef[i];      
      FindMaxRad();
      rrod = VPf.ballf * VPf.rodf * DEF_ATRADF * rcovdef[1];
      VPf.atradf = DEF_ATRADF;
      remake3D = REMAKE3D_BALL | REMAKE3D_SPACE | REMAKE3D_BOND |
	REMAKE3D_BALL_LABEL | REMAKE3D_SPACE_LABEL;
      break;

    case R_TESSELLATION:
      VPf.tessFactor    = DEF_TESSELLATION;
      remake3D = REMAKE3D_BALL | REMAKE3D_SPACE | REMAKE3D_BOND;
      break;

    case R_UNIBOND:
      VPf.unibond = DEF_UNIBOND;
      remake3D = REMAKE3D_BOND;
      break;

    case R_BACKGROUND:
      for (i=0; i<4; i++)
	bg[i] = DefBg[i];
      break;

    case R_PERSPECTIVE:
      VPf.perspective = DEF_PERSPECTIVE;
      break;

    case R_PERSPECTIVEBACK:
      VPf.perspective_far = DEF_PERSPECTIVEBACK;
      break;

    case R_PERSPECTIVEFOVY:
      VPf.perspective_fovy = DEF_PERSPECTIVEFOVY;
      break;

    case R_PERSPECTIVEFRONT:
      VPf.perspective_size = DEF_PERSPECTIVEFRONT;
      break;

    case R_FOG:
      if ( togl == mesa_togl ) {
	VPf.fog = DEF_FOG;
      } else {
	NEW_WIN_CONTEXT *wc;
	wc = FindWinContextByTogl( togl );
	wc->VPf.fog = DEF_FOG;
      }
      break;

    case R_AMBIENT_BY_DIFFUSE:
      AtCol_Ambient_by_Diffuse = DefAtCol_Ambient_by_Diffuse;
      break;
      
    case R_ANTIALIAS:
      if ( togl == mesa_togl ) {
	VPf.antialias = DEF_ANTIALIAS;
      } else {
	NEW_WIN_CONTEXT *wc;
	wc = FindWinContextByTogl( togl );
	wc->VPf.antialias = DEF_ANTIALIAS;
      }
      break;

    case R_RBALL:
      for(i=0; i <= MAXNAT; i++)
	rball[i] = DEF_BALLF * DEF_ATRADF * rcovdef[i];
      VPf.ballf = DEF_BALLF;
      rrod = VPf.ballf * VPf.rodf * DEF_ATRADF * rcovdef[1];
      remake3D = REMAKE3D_BALL | REMAKE3D_BOND | REMAKE3D_BALL_LABEL;
      break;
    
    case R_RROD:
      rrod = DEF_BALLF * DEF_ATRADF * DEF_RODF * rcovdef[1];
      VPf.rodf = DEF_RODF;
      remake3D = REMAKE3D_BOND;
      break;
    
    case R_ATCOL:
      for(i = 0; i <= MAXNAT; i++)
	for(j = 0; j < 3; j++)
	  atcol[i][j] = DefAtCol[i][j];
      LoadAtmCol(0);
      remake3D = REMAKE3D_BALL | REMAKE3D_SPACE | REMAKE3D_BOND;
      break;
    
    case R_WFLINEWIDTH:
      VPf.WFlinewidth = DEF_LineWidthWF;
      break;
    
    case R_WF3DLINEWIDTH:
      VPf.WF3Dlinewidth = DEF_LineWidthWF3D;
      break;
    
    case R_OUTLINEWIDTH:
      VPf.OUTlinewidth = DEF_LineWidthOUT;
      break;
    
    case R_PLLINEWIDTH:
      VPf.PLlinewidth = DEF_LineWidthPL;
      break;
    
    case R_PLRADIUS:
      VPf.PLradius = DEF_RadiusPL;
      break;
    
    case R_SCRF:
      /* this is no more used */
      /* 
	 VPf.scrf = DEF_SCRF;	 
	 if (VPf.stropened) {
	 xcViewPort();
	 Togl_PostRedisplay(togl);
	 }
      */
      return;
    
    case R_ALL:
      DefValues();
      remake3D = REMAKE3D_ALL;
      break;
      
    case R_FRAMECOL:
      for(j = 0; j < 3; j++)
	  framecol[j] = DefFrameCol[j];
      break;

    case R_UNIBONDCOLOR:
      for(j = 0; j < 3; j++)
	unibondCol[j] = DefUnibondCol[j];
      remake3D = REMAKE3D_BOND;
      break;

    case R_FRAMELINEWIDTH:
      VPf.framewidth = DEF_FRAMEWIDTH;
      break;

    case R_FRAMERODF:
      rframe = DEF_FRAMERODF * rcovdef[1];
      remake3D = REMAKE3D_FRAME;
      break;

    case R_CURRENTFILEFORMAT:
      current_file_format = FORMAT_NONE;
      break;
      
    case R_FORCE_RODTHICKF:
      FV.rod_thickf   = VECTOR_THICKF;
      BuildForceVectors(&FV);
      UpdateProjection();    
      break;

    case R_FORCE_ARRTHICKF:
      FV.arr_thickf   = VECTOR_ARRTHICKF;
      BuildForceVectors(&FV);
      UpdateProjection();    
      break;

    case R_FORCE_ARRLENF:
      FV.arr_lenf     = VECTOR_ARROWSIZE;
      BuildForceVectors(&FV);
      UpdateProjection();    
      break;

    case R_FORCE_COLOR:
      COPY_V(3, FV.color, DefAtCol[0]);
      FV.color[3] = 1.0;
      break;

    case R_XYZ_AXIS_COLOR:
      COPY_V(4, xyz.axis_color, def_xyz_axis_color);
      break;

    case R_XYZ_XYPLANE_COLOR:
      COPY_V(4, xyz.xyplane_color, def_xyz_xyplane_color);
      break;
    }

  if (VPf.stropened) {
    /* I must update a Display */
    ReDisplay(togl, remake3D);
  }
}


void
ReloadVars(void)
{
  /* T.K */
  VPf.xyzOn = 0;
  VPf.labelsOn = 0;
  VPf.framesOn = 0;
  VPf.selection = 0;
  natoms = 0;
  nbonds = 0;
  nframes = 0;
  nobjects = 0;
  tmp_nobjects = 0;
  tr.xtransl = 0.0;
  tr.ytransl = 0.0;
  tr.zoom = DEF_ZOOM;
  tr.b1motion = 0;
  tr.b2motion = 0;
  tr.shiftB1motion = 0;
}


/*****************************************************************************
 * this is defined in "struct.h"
 * ----------------------------
 * define L_SPACE_COV       0 atrad = rcov (effect to all elements - SPACEFILL)
 * define L_SPACE_VDW       1 atrad = rvdw (effect to all elements - SPACEFILL)
 * define L_BALL_COV        2 atrad = rcov (effect to all elements - BALLS)
 * define L_BALL_VDW        3 atrad = rvdw (effect to all elements - BALLS) 
 * define L_ATRAD_ONE       4 change just one atmoic radius
 * define L_ATRAD_SCALE     5 change all atomic radius for some scale factor
 * define L_BALLF           6 change ball factor
 * define L_RODF            7 change rod factor
 * define L_ATCOL_ONE       8 change color of one atom
 * define L_WFLINEWIDTH     9
 * define L_PLLINEWIDTH     10
 * define L_PLRADIUS        11
 * define L_SCRF            12
 * define L_COV_SCALE       13
 * define L_XYZ_ON          14  for coordinate sistem  
 * define L_LABELS_ON       15  for atomic labels 
 * define L_FRAME_ON        16
 * define L_FRAMECOL        17 width of line-frame 
 * define L_FRAMELINEWIDTH  18
 * define L_FRAMERODF       19 width of rod-frame
 * define L_LINEFRAME       20 display frames as rods or as lines when in 3D
 *                             mode
 * define L_TR_XTRANSL      21
 * define L_TR_YTRANSL      22
 * define L_BACKGROUND      23 
 *------------------------------
 */

void
LoadNewValue(struct Togl *togl, int var, double value1, double value2, 
	     double value3, double value4)
{
  int i;
  GLuint remake3D = REMAKE3D_NONE;

  switch(var)
    {
    case L_SPACE_COV:
      for(i=1; i <= MAXNAT; i++) 
	atrad[i] = VPf.atradf * rcovdef[i]; 
      FindMaxRad();
      remake3D = REMAKE3D_SPACE | REMAKE3D_SPACE_LABEL;
      break;

    case L_BALL_COV:
      printf("L_BALL_COV\n");
      fflush(stdout);
      for(i=0; i <= MAXNAT; i++)
	rball[i] = VPf.ballf *  VPf.atradf * rcovdef[i];
      rrod = VPf.rodf * rball[1];
      remake3D = REMAKE3D_BALL | REMAKE3D_BOND | REMAKE3D_BALL_LABEL;
      break;
    
    case L_SPACE_VDW:
      printf("L_SPACE_VDW\n");
      fflush(stdout);
      for(i=1; i <= MAXNAT; i++) 
	atrad[i] = VPf.atradf * rvdw[i];
      FindMaxRad();
      remake3D = REMAKE3D_SPACE | REMAKE3D_SPACE_LABEL;
      break;

    case L_BALL_VDW:
      printf("L_BALL_VDW\n");
      fflush(stdout);
      for(i=0; i <= MAXNAT; i++) 
	rball[i] = VPf.ballf * VPf.atradf * rvdw[i];
      rrod = VPf.rodf * rball[1];
      remake3D = REMAKE3D_BALL | REMAKE3D_BOND | REMAKE3D_BALL_LABEL;
      break;
    
      /* t.k.; note: should make a rcov_custom variable ... */
    case L_RCOV_ONE:
      rcov[(int) value1] = value2;      
      remake3D = REMAKE3D_BOND;
      if (VPf.stropened) 
	/* it's illegal to ReMakeStr if structure is not opened */
	ReMakeStr();
      break;
      /* /t.k. */

    case L_ATRAD_ONE:
      atrad[(int) value1] = value2; 
      FindMaxRad();
      rball[(int) value1] = atrad[(int) value1] * VPf.ballf;
      remake3D = REMAKE3D_SPACE | REMAKE3D_BALL | REMAKE3D_BALL_LABEL |
	REMAKE3D_SPACE_LABEL;
      if (VPf.stropened) 
	/* it's illegal to ReMakeStr if structure is not opened */
	ReMakeStr();
      break;
    
    case L_ATRAD_SCALE:
      printf("L_ATRAD_SCALE:: value1 = %f\n",value1);
      for(i=0; i <= MAXNAT; i++) {
	atrad[i]  /= VPf.atradf; 
	rball[i]  /= VPf.atradf;
	atrad[i]  *= value1;
	rball[i]  *= value1;
      }
      FindMaxRad();
      rrod /= VPf.atradf;
      rrod *= value1;
      VPf.atradf = value1;	
      remake3D = REMAKE3D_SPACE | REMAKE3D_BALL | REMAKE3D_BOND |
	REMAKE3D_BALL_LABEL | REMAKE3D_SPACE_LABEL;      
      break;

    case L_TESSELLATION:
      VPf.tessFactor = MAX( value1, 2.0 );
      remake3D = REMAKE3D_SPACE | REMAKE3D_BALL | REMAKE3D_BOND;
      break;

    case L_UNIBOND:
      VPf.unibond = (GLboolean) value1;
      remake3D = REMAKE3D_BOND;
      break;
	
    case L_UNIBONDCOLOR:
      unibondCol[0] = value1;
      unibondCol[1] = value2;
      unibondCol[2] = value3;
      remake3D = REMAKE3D_BOND;
      break;      

    case L_FOG:
      if ( togl == mesa_togl ) {
	VPf.fog = (GLboolean) value1;
      } else {
	NEW_WIN_CONTEXT *wc;
	wc = FindWinContextByTogl( togl );
	wc->VPf.fog = (GLboolean) value1;
      }
      break;

    case L_AMBIENT_BY_DIFFUSE:
      AtCol_Ambient_by_Diffuse = (GLfloat) value1;
      break;

    case L_ANTIALIAS:
      if ( togl == mesa_togl ) {
	VPf.antialias = (GLboolean) value1;
      } else {
	NEW_WIN_CONTEXT *wc;
	wc = FindWinContextByTogl( togl );
	wc->VPf.antialias = (GLboolean) value1;
      }
      break;

    case L_PERSPECTIVE:
      VPf.perspective = (GLboolean) value1;
      break;

    case L_PERSPECTIVEBACK:
      VPf.perspective_far = value1;
      break;
      
    case L_PERSPECTIVEFOVY:
      VPf.perspective_fovy = value1;
      break;
      
    case L_PERSPECTIVEFRONT:
      VPf.perspective_size = value1;
      break;
      
    case L_BALLF:
      VPf.ballf = value1;
      for(i=0; i <= MAXNAT; i++) 
	rball[i] = VPf.ballf * atrad[i];
      rrod = VPf.rodf * rball[1];
      remake3D = REMAKE3D_BALL | REMAKE3D_BOND | REMAKE3D_BALL_LABEL;
      break;
    
    case L_RODF:
      VPf.rodf = value1;
      rrod = VPf.atradf * VPf.ballf * VPf.rodf * atrad[1];      
      remake3D = REMAKE3D_BOND;
      break;

    case L_ATCOL_ONE:
      atcol[(int) value1][0] = (float) value2;
      atcol[(int) value1][1] = (float) value3;
      atcol[(int) value1][2] = (float) value4;	
      LoadAtmCol(0);
      remake3D = REMAKE3D_SPACE | REMAKE3D_BALL | REMAKE3D_BOND;
      printf("L_ATCOL_ONE\n");
      printf("color = %f %f %f\n",atcol[(int) value1][0],
	     atcol[(int) value1][1],
	     atcol[(int) value1][2]);
      fflush(stdout);      
      break;

    case L_WFLINEWIDTH:
      VPf.WFlinewidth = (int) value1;
      break;

    case L_WF3DLINEWIDTH:
      VPf.WF3Dlinewidth = (int) value1;
      break;

    case L_OUTLINEWIDTH:
      VPf.OUTlinewidth = (int) value1;
      break;

    case L_PLLINEWIDTH:
      VPf.PLlinewidth = (int) value1;
      break;

    case L_PLRADIUS:
      VPf.PLradius = (int) value1;
      xcMakePointList();
      break;

    case L_SCRF:
      /* this is no more used */
      /* VPf.scrf = value1; */
      break;
      
    case L_COV_SCALE:
      VPf.rcovf = value1;
      printf("L_COV_SCALE:: VPf.rcovf = %f\n",VPf.rcovf);
      /* this one is complicated !!!!
       * if this is changed, than chemical connectivity criteria is changed
       * this means that we must reassing bonds; 
       *
       * Side Effects:
       * ------------
       * Bond criteria is changed -> that reflect in new connectivity, but
       * atrad, ballr, etc. are not changed; L_COV_SCALE is just "change
       * chemical connecivity criteia" !!!!
       */
      for(i=1; i <= MAXNAT; i++) 
	rcov[i] = VPf.rcovf * rcovdef[i];
      remake3D = REMAKE3D_BOND;
      if (VPf.stropened) 
	/* it's illegal to ReMakeStr if structure is not opened */
	ReMakeStr();
      break;

    case L_XYZ_ON:
      printf("L_XYZ_ON> value1=%f\n",value1);
      fflush(stdout);
      if ( value1 > 0.0 ) {
	VPf.xyzOn = 1;
	if ( VPf.stropened ) Togl_PostRedisplay(togl);
      }
      else xcTurnXYZOff(togl);
      return;

    case L_LABELS_ON:
      if ( value1 > 0.0 ) {
	VPf.labelsOn = 1;
	if ( VPf.stropened ) {
	  /* maybe LabelList is not yet made */
	  /* xcMakeBallLabel3DList();  if LabelList is not yet made, than this 
				     * fn. will make it, else will do nothing 
				     */
/*  	  xcMakeSpaceLabel3DList(); */
	  xcTurnLabelsOn();
	}
      }
      else xcTurnLabelsOff();      
      break;

    case L_FRAME_ON:
      if ( value1 > 0.0 ) {
	VPf.framesOn = 1;
	if ( VPf.stropened ) {
	  /* maybe FrameList is not yet made */
	  /* xcMakeFrame3DLists();  if FrameList is not made yet, than this 
				    fn. will make it, else will do nothing */
	  xcTurnFramesOn();
	}
      }
      else xcTurnFramesOff();      
      break;

    case L_FRAMECOL:
      framecol[0] = value1;
      framecol[1] = value2;
      framecol[2] = value3;
      break;
   
    case L_FRAMELINEWIDTH:
      VPf.framewidth = (int) value1;
      break;

    case L_FRAMERODF:
      VPf.framef = value1;
      rframe = VPf.framef * rcovdef[1];
      remake3D = REMAKE3D_FRAME;
      break;
      
    case L_LINEFRAME:
      /* if value1 is 1 than is.lineframe is set to 1 */
      if (value1 > 0.0 ) is.lineframe = 1;
      else is.lineframe = 0;
      /* now create xcMakeFrame3DLists, if it's not created yet */
      /*xcMakeFrame3DLists();*/ /* this is also done in L_FRAMEON !!!
			         * Here is just in case 
			         */
      break;
      
    case L_TR_XTRANSL:
      tr.xtransl = value1;
      xcViewPort();
      Togl_PostRedisplay(togl);
      return;

    case L_TR_YTRANSL:
      tr.ytransl = value1;
      xcViewPort();
      Togl_PostRedisplay(togl);
      return;

    case L_BACKGROUND:
      if ( togl == mesa_togl ) {
	bg[0] = value1;
	bg[1] = value2;
	bg[2] = value3;
	bg[3] = value4;
	xcChangeBackground( togl, bg );
      } else {
	NEW_WIN_CONTEXT *wc;
	/* not .mesa */
	wc = FindWinContextByTogl( togl );
	wc->bg[0] = value1;
	wc->bg[1] = value2;
	wc->bg[2] = value3;
	wc->bg[3] = value4;	
	xcChangeBackground( togl, wc->bg );
      }
      return;    

    case L_CURRENTFILEFORMAT:
      current_file_format = (int) value1;
      break;

    case SET_FOG_COLORMODE:
      fog.colormode = (int) value1;
      break;

    case SET_FOG_COLOR:
      fog.color[0] = (GLfloat) value1;
      fog.color[1] = (GLfloat) value2;
      fog.color[2] = (GLfloat) value3;
      fog.color[3] = (GLfloat) value4;
    case SET_FOG_MODE:
      fog.mode = value1;
      break;

    case SET_FOG_DENSITY:
      fog.density = (GLfloat) value1;
      break;

    case SET_FOG_ORT_START_F:
      fog.ort_start_f = (GLfloat) value1;
      break;

    case SET_FOG_ORT_END_F:
      fog.ort_end_f = (GLfloat) value1;
      break;

    case SET_FOG_PERSP_F1:
      fog.persp_f1 = (GLfloat) value1;
      break;

    case SET_FOG_PERSP_F2:
      fog.persp_f2 = (GLfloat) value1;
      break;

    case SET_ANTIALIAS_DEGREE:
      antialias.degree = (int) value1;
      break;

    case SET_ANTIALIAS_OFFSET:
      antialias.offset = (float) value1;
      break;

    case L_FORCE_RODTHICKF:
      FV.rod_thickf   = (float) value1;
      BuildForceVectors(&FV);
      UpdateProjection();    
      break;

    case L_FORCE_ARRTHICKF:
      FV.arr_thickf   = (float) value1;
      BuildForceVectors(&FV);
      UpdateProjection();    
      break;

    case L_FORCE_ARRLENF:
      FV.arr_lenf     = (float) value1;
      BuildForceVectors(&FV);
      UpdateProjection();    
      break;

    case L_FORCE_COLOR:
      FV.color[0] = (float) value1;
      FV.color[1] = (float) value2;
      FV.color[2] = (float) value3;
      FV.color[3] = (float) value4;
      break;

    case R_XYZ_AXIS_COLOR:
      xyz.axis_color[0] = (float) value1;
      xyz.axis_color[1] = (float) value2;
      xyz.axis_color[2] = (float) value3;
      xyz.axis_color[3] = (float) value4;
      break;

    case R_XYZ_XYPLANE_COLOR:
      xyz.xyplane_color[0] = (float) value1;
      xyz.xyplane_color[1] = (float) value2;
      xyz.xyplane_color[2] = (float) value3;
      xyz.xyplane_color[3] = (float) value4;
      break;
    }

  if (VPf.stropened)
    /* now reDisplay */
    ReDisplay(togl, remake3D);
}


/*****************************************************************************
 * D_ for "get default values"
 * ----------------------------
 * define D_SCRF            0
 * define D_BALLF           1
 * define D_RODF            2
 * define D_COVF            3
 * define D_ALL             4
 * define D_WFLINEWIDTH     5
 * define D_PLLINEWIDTH     6
 * define D_PLRADIUS        7	
 * define D_ATCOL_ONE       8    [this is done in XC_GetDefaultCmd]
 * define D_ATRAD_SCALE     9
 * define D_ATRAD_ONE      10
 * define D_FRAMECOL       11    [this is done in XC_GetDefaultCmd]
 * define D_FRAMELINEWIDTH 12
 * define D_FRAMERODF      13
 *------------------------------
 */
double
GetDefault(int var)
{
  switch(var)
    {
    case D_BALLF:
      return DEF_BALLF;

    case D_RODF:
      return DEF_RODF;

    case D_COVF:
      return DEF_RCOVF;

    case D_WFLINEWIDTH:
      return DEF_LineWidthWF;
    
    case D_WF3DLINEWIDTH:      return DEF_LineWidthWF3D;
    
    case D_OUTLINEWIDTH:
      return DEF_LineWidthOUT;
    
    case D_PLLINEWIDTH:
      return DEF_LineWidthPL;
    
    case D_PLRADIUS:
      return DEF_RadiusPL;
    
    case D_SCRF:
      /* this is no more used */
      /* return DEF_SCRF; */
      return 1.0;

    case D_ATRAD_SCALE:
      return DEF_ATRADF;
      
    case D_TESSELLATION:
      return DEF_TESSELLATION;

    case D_UNIBOND:
      return DEF_UNIBOND;

    case D_FOG:
      return DEF_FOG;

    case D_AMBIENT_BY_DIFFUSE:
      return DefAtCol_Ambient_by_Diffuse;
      
    case D_ANTIALIAS:
      return DEF_ANTIALIAS;

    case D_PERSPECTIVE:
      return DEF_PERSPECTIVE;

    case D_PERSPECTIVEBACK:
      return DEF_PERSPECTIVEBACK;
      
    case D_PERSPECTIVEFOVY:
      return DEF_PERSPECTIVEFOVY;
      
    case D_PERSPECTIVEFRONT:
      return DEF_PERSPECTIVEFRONT;
      
    case D_FRAMELINEWIDTH:
      return DEF_FRAMEWIDTH;

    case D_FRAMERODF:
      return DEF_FRAMERODF;

    case D_CURRENTFILEFORMAT:
      current_file_format = FORMAT_NONE;
      break;

    case GET_FOG_COLORMODE:
      return def_fog_colormode;

    case GET_FOG_MODE:
      return def_fog_mode;

    case GET_FOG_DENSITY:
      return def_fog_density;

    case GET_FOG_ORT_START_F:
      return def_fog_ort_start_f;

    case GET_FOG_ORT_END_F:
      return def_fog_ort_end_f;

    case GET_FOG_PERSP_F1:
      return def_fog_persp_f1;

    case GET_FOG_PERSP_F2:
      return def_fog_persp_f2;

    case GET_ANTIALIAS_DEGREE:
      return def_antialias_degree;

    case GET_ANTIALIAS_OFFSET:
      return def_antialias_offset;

    case D_FORCE_RODTHICKF:
      return VECTOR_THICKF;

    case D_FORCE_ARRTHICKF:
      return VECTOR_ARRTHICKF;

    case D_FORCE_ARRLENF:
      return VECTOR_ARROWSIZE;
    } 
  return 0.0;
}


/*****************************************************************************/
/* GetValue use the same constants as GetDefault !!!!!!!!!!
 */
double 
GetValue(int var) 
{ 
  switch(var) 
    { 
    case D_BALLF:  
      return VPf.ballf;

    case D_RODF: 
      return VPf.rodf;

    case D_COVF: 
      return VPf.rcovf;

    case D_WFLINEWIDTH: 
      return VPf.WFlinewidth;
    
    case D_WF3DLINEWIDTH: 
      return VPf.WF3Dlinewidth;
    
    case D_OUTLINEWIDTH: 
      return VPf.OUTlinewidth;
    
    case D_PLLINEWIDTH: 
      return VPf.PLlinewidth;
    
    case D_PLRADIUS: 
      return VPf.PLradius;
    
    case D_SCRF:
      /* this is no more used */
      /* return VPf.scrf; */
      return 1.0;
     
    case D_ATRAD_SCALE:      
      return VPf.atradf;

    case D_TESSELLATION:
      return VPf.tessFactor;
      
    case D_UNIBOND:
      return VPf.unibond;
      
    case D_FOG:
      return VPf.fog;
      
    case D_AMBIENT_BY_DIFFUSE:
      return AtCol_Ambient_by_Diffuse;

    case D_ANTIALIAS:
      return VPf.antialias;      

    case D_PERSPECTIVE:
      return VPf.perspective;
      
    case D_PERSPECTIVEBACK:
      return VPf.perspective_far;
      
    case D_PERSPECTIVEFOVY:
      return VPf.perspective_fovy;
      
    case D_PERSPECTIVEFRONT:
      return VPf.perspective_size;
      
    case D_FRAMELINEWIDTH:
      return VPf.framewidth;

    case D_FRAMERODF:
      return VPf.framef;

    case D_CURRENTFILEFORMAT:
      return current_file_format;

    case GET_SS_MINX:
      return ss.minX;

    case GET_SS_MINY:
      return ss.minY;

    case GET_SS_MINZ:
      return ss.minZ;

    case GET_SS_MAXX:
      return ss.maxX;
      
    case GET_SS_MAXY:
      return ss.maxY;
      
    case GET_SS_MAXZ:
      return ss.maxZ;

    case GET_AT_MINX:
      return min.x;

    case GET_AT_MINY:
      return min.y;

    case GET_AT_MINZ:
      return min.z;

    case GET_AT_MAXX:
      return max.x;
      
    case GET_AT_MAXY:
      return max.y;
      
    case GET_AT_MAXZ:
      return max.z;

    case GET_FOG_COLORMODE:
      return fog.colormode;

    case GET_FOG_MODE:
      return fog.mode;

    case GET_FOG_DENSITY:
      return fog.density;

    case GET_FOG_ORT_START_F:
      return fog.ort_start_f;

    case GET_FOG_ORT_END_F:
      return fog.ort_end_f;

    case GET_FOG_PERSP_F1:
      return fog.persp_f1;

    case GET_FOG_PERSP_F2:
      return fog.persp_f2;

    case GET_ANTIALIAS_DEGREE:
      return antialias.degree;

    case GET_ANTIALIAS_OFFSET:
      return antialias.offset;

    case GET_ALAT:
      if (xcr.lconvvec) {
	return NORM_V3 (vec.conv[0]);
      } else if (xcr.lprimvec) {
	return NORM_V3 (vec.prim[0]);
      } 
      return 1.0;
	
    case D_FORCE_RODTHICKF:
      return FV.rod_thickf;
      
    case D_FORCE_ARRTHICKF:
      return FV.arr_thickf;

    case D_FORCE_ARRLENF:
      return FV.arr_lenf;
    } 
  return 0.0;
} 


void
LoadOldAtomicColors(void)
{
  int i,j;  
  for(i=0; i <= MAXNAT; i++)
    for(j=0; j < 3; j++)
      atcol[i][j] = DefAtColOld[i][j];    
  LoadAtmCol(0);
}
