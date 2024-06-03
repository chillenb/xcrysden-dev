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
* Source: $XCRYSDEN_TOPDIR/C/xcFont.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

TODO: 

4. make a routine for mapping TkFont -> XLFD and finding correct 
normal/bold/medium/oblique/italic/roman       
*/

#ifndef WIN32
#  define TOGL_X11
#  define X11
#else
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  undef WIN32_LEAN_AND_MEAN
#  include <winnt.h>
#endif

#include <togl.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <tk.h>
#ifdef WIN32
#  include <tkPlatDecls.h>
#endif

#ifndef WIN32
#  include <X11/Xutil.h>
#endif

#include "struct.h"
#include "xcGLparam.h"
#include "memory.h"
#include "xcfunc.h"

#define XC_FONT_TKYES_XNO 999

extern GLuint         fontOffset;
extern RasterFontSize rf;
AtomicLabel globalAtomLabel = { 
  0,  0,  0,
  {1.0, 1.0, 1.0},
  {0.0, 0.0, 0.0},
  XC_YES, (char*)NULL, (Tk_Font)NULL
};
AtomicLabel *atomLabel = (AtomicLabel*)NULL;
short *do_not_display_atomlabel = (short*)NULL; /* apply to all kind
						   of labels default
						   and custom  */

struct xcToglFont {
  char              *font;
  Tk_Font            tkfont;
  GLuint             base;
  int                height;
  int                width;
  struct xcToglFont *prev;
};
static struct xcToglFont *xcFont, *xcFontPtr = (struct xcToglFont*) NULL;

typedef struct {
  char    *font;
  char    *string;
  Tk_Font tkfont;
  int     height;
  int     width;
} FontType;


/* xcFont.c */
int XC_SetFontCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
int XC_SetAtomLabelCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
int XC_ClearAtomLabelCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
int XC_QueryFontCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
void xcFont_PrintString(const char *s);
void xcTkFontFreeAll(void);
GLuint Togl_LoadBitmapFontOld(const Togl *togl, const char *fontname);
void Togl_UnloadBitmapFontOld(const Togl *togl, GLuint fontbase);

static void _clearAtomLabel(int atomID);
static int _setFontColor(Tcl_Interp *interp, int argc, int bright_ind, const char *bright_argv, int dark_ind, const char *dark_argv, AtomicLabel *alabel);
static void xcFont_addNew(struct xcToglFont *f, struct Togl *togl, FontType this);
static struct xcToglFont *xcFont_find(char *font);
static int _assignFont(int font_found, struct Togl *togl, const char *fontString, FontType thisFont);


/* ------------------------------------------------------------------------
 *
 * xc_setfont $togl XLFD_fontname|tkfontname ?bright-fontcolor? ?dark-fontcolor?
 *
 * Format of bright-fontcolor & dark-fontcolor is RGB, where the
 * components are within [0-1].
 * ------------------------------------------------------------------------ */

int
XC_SetFontCb(ClientData clientData, Tcl_Interp *interp,
	     int argc, const char *argv[])
{
  Togl *togl;
  FontType          thisFont;
  int               result, font_found = 0;

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if (argc < 3 || argc > 5) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "wrong number of arguments, should be xc_setfont togl togl_font ?bright-fontcolor? ?dark-fontcolor?");
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* assign the font color first */

  if ( _setFontColor(interp, argc, 3, argv[3], 4, argv[4], &globalAtomLabel) == TCL_ERROR ) {
    return TCL_ERROR;
  }

  /* now parse the font-name */

  if (strlen(argv[2]) > 0) {
    /*
      int i;
      for (i=0; i < NFONTS; i++) {
      if ( strcmp (argv[2], fontType[i].string) == 0 ) {
      thisFont.font   = fontType[i].font;
      thisFont.height = fontType[i].height;
      thisFont.width  = fontType[i].width;
      font_found    = 1;
      break;
      }
      }
    */

    result = _assignFont(font_found, togl, argv[2], thisFont);
    /*fprintf(stderr,"_assignFont = %d\n", result);*/
    if ( result == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "couldn't load font %s", argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    } else if ( result == XC_FONT_TKYES_XNO ) {
      return TCL_OK;
    } else {
      /* xcFont should exists by now */
      if ( !xcFont->base ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "couldn't load font %s", argv[2]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }      

    /* the font setting was successful, update the globalAtomLabel */
    globalAtomLabel.base   = xcFont->base;
    globalAtomLabel.tkfont = xcFont->tkfont;
    globalAtomLabel.height = xcFont->height;
    globalAtomLabel.width  = xcFont->width;
    globalAtomLabel.do_display = XC_YES;

    /*    
	  fprintf(stderr, "FONT-INFO:: (xcFontPtr=%d)\n", xcFontPtr);
	  fprintf(stderr, "  base=%d\n  tkfont=%d\n  h=%d, w=%d\n",
	  xcFont->base, &xcFont->tkfont, xcFont->height, xcFont->width);    
    */
  }

  Togl_PostRedisplay (togl);
  return TCL_OK;
}


/* ------------------------------------------------------------------------

   xc_setatomlabel  togl  atomID  labelString  ?XLFD_fontname?  ?bright-fontcolor? ?dark-fontcolor?

   dark-fontcolor   == color for ballstick+spacefill for lighting-Off display modes
   bright-fontcolor == color for the rest of display modes

   ------------------------------------------------------------------------ */

int
XC_SetAtomLabelCb(ClientData clientData, Tcl_Interp *interp,
		  int argc, const char *argv[])
{
  Togl *togl;
  FontType   thisFont;
  int        atomID;

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if (argc < 4 || argc > 7) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "wrong number of arguments, must be xc_setatomlabel  togl  atomID  labelString  ?font?  ?fontColor1?  ?fontColor2?");
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &atomID) == TCL_ERROR) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "wanted integer but got %s", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( atomID < 1 || atomID > natoms ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "atomID %d out of range, should be within [1,%d]", atomID, natoms);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* copy the atom label */
  atomLabel[atomID].label = (char*) xcRealloc ( atomLabel[atomID].label, sizeof(char)*((size_t) strlen(argv[3])+1) );
  strcpy(atomLabel[atomID].label, argv[3]);

  if ( argc >= 5 && strlen(argv[4]) > 0 ) {
    int result;

    result = _assignFont(0, togl, argv[4], thisFont);
    /*fprintf(stderr,"_assignFont = %d\n", result);*/
    if ( result == TCL_ERROR || !xcFont->base) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "couldn't load font %s", argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    } else if ( result == XC_FONT_TKYES_XNO ) {
      return TCL_OK;
    }

    /* now assign the font */

    /* fprintf(stderr,"1-base:: %d\n", xcFont->base); */

    atomLabel[atomID].base   = xcFont->base;
    atomLabel[atomID].tkfont = xcFont->tkfont;
    atomLabel[atomID].height = xcFont->height;
    atomLabel[atomID].width  = xcFont->width;    
  } else if (!xcFontPtr) {
    /* default font */

    /* fprintf(stderr,"2-base:: %d\n", fontOffset); */

    atomLabel[atomID].base   = fontOffset;
    atomLabel[atomID].tkfont = (Tk_Font)NULL;
    atomLabel[atomID].height = rf.height;
    atomLabel[atomID].width  = rf.wid; 
  } else {
    /* assign the last font loaded global font */

    /* fprintf(stderr,"3-base:: %d\n", xcFontPtr->base); */

    atomLabel[atomID].base   = xcFontPtr->base;
    atomLabel[atomID].tkfont = xcFontPtr->tkfont;
    atomLabel[atomID].height = xcFontPtr->height;
    atomLabel[atomID].width  = xcFontPtr->width;    
  }
  atomLabel[atomID].do_display = XC_YES;

  /* assign the font color */

  if ( _setFontColor(interp, argc, 5, argv[5], 6, argv[6], atomLabel + atomID) == TCL_ERROR ) {
    return TCL_ERROR;
  }

  Togl_PostRedisplay (togl);
  return TCL_OK;
}



/* ------------------------------------------------------------------------

   xc_clearatomlabel  togl atomID  ?atomID? ...

   or

   xc_clearatomlabel togl all

   ------------------------------------------------------------------------ */

int
XC_ClearAtomLabelCb(ClientData clientData, Tcl_Interp *interp,
		    int argc, const char *argv[])
{
  Togl *togl;
  int i;

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if (argc < 3 ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "Usage:\n   xc_clearatomlabels  togl  atomID  ?atomID? ...\n or \n   togl xc_cleanatomlabels all");
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( strcmp(argv[2], "all") == 0 ) {
    /* the "xc_cleanatomlabels all form" */
    for (i=1; i<=natoms; i++) _clearAtomLabel(i);
  } 

  else {
    int atomID;  
    for (i=2; i<argc; i++) {
      if (Tcl_GetInt(interp, argv[i], &atomID) == TCL_ERROR) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "wanted integer but got %s", argv[i]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      if ( atomID < 1 || atomID > natoms ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "atom index %d out of range, should be between [0,%d]", 
		 atomID, natoms);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      _clearAtomLabel(atomID);
    }
  }

  Togl_PostRedisplay (togl);
  return TCL_OK;
}

static void _clearAtomLabel(int atomID) {  
  atomLabel[atomID].base       = 0;
  atomLabel[atomID].do_display = XC_NO;
  xcFree(atomLabel[atomID].label);
  atomLabel[atomID].label      = (char*)NULL; 
}




/* ------------------------------------------------------------------------
 *
 * $togl xc_queryfont togl XLFD_fontname
 *
 * Format of bright-fontcolor & dark-fontcolor is RGB, where the
 * components are within [0-1].
 * ------------------------------------------------------------------------ */

int
XC_QueryFontCb(ClientData clientData, Tcl_Interp *interp,
	       int argc, const char *argv[])
{
  Togl *togl;
  FontType   thisFont;
  int        result, font_found = 0;

  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if (argc != 3) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "wrong number of arguments, should be xc_queryfont togl XLFD_fontname");
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  result = _assignFont(font_found, togl, argv[2], thisFont);

  /* fprintf(stderr,"_assignFont = %d\n", result); */

  if ( result == XC_FONT_TKYES_XNO ) {
    return TCL_OK;
  } else if ( result == TCL_ERROR ) {
    char *result = Tcl_Alloc (sizeof(char) * 3);      
    sprintf(result, "-1");
    Tcl_SetResult(interp, result, TCL_DYNAMIC);    
  } else {
    /* xcFont should exists by now */
    if ( xcFont->base ) {
      char *result = Tcl_Alloc (sizeof(char) * 3);      
      sprintf(result, "+1");
      Tcl_SetResult(interp, result, TCL_DYNAMIC);    
    } else {
      char *result = Tcl_Alloc (sizeof(char) * 3);      
      sprintf(result, "-1");
      Tcl_SetResult(interp, result, TCL_DYNAMIC);    
    }
  }    
  return TCL_OK;
}



void xcFont_PrintString (const char *s) {
  glCallLists( strlen(s), GL_UNSIGNED_BYTE, s );
}



void xcTkFontFreeAll(void) {
  struct xcToglFont *f = xcFontPtr;
  while (f) {
    if (f->tkfont) Tk_FreeFont (f->tkfont);
    f = f->prev;
  }
}



static int _setFontColor(Tcl_Interp *interp, int argc, 
			 int bright_ind, const char *bright_argv, 
			 int dark_ind, const char *dark_argv, AtomicLabel *alabel) {
  GetGlParam  bright, dark;

  if ( argc > bright_ind && strlen(bright_argv) > 0 ) {
    if ( ! xcSplitList( XC_GET_RGB, interp, &bright_argv, &bright ) ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "error parsing bright font color, should be {red green blue}, but got %s", bright_argv);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    COPY_V(3, alabel->bright_color, bright.vec);
  }

  if ( argc > dark_ind && strlen(dark_argv) > 0 ) {
    if ( ! xcSplitList( XC_GET_RGB, interp, &dark_argv, &dark ) ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "error parsing dark font color, should be {red green blue}, but got %s", dark_argv);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    COPY_V(3, alabel->dark_color, dark.vec);
  }

  return TCL_OK;
}


static void xcFont_addNew(struct xcToglFont *f, struct Togl *togl, FontType this) {
  f->font   = this.font;
  f->base   = Togl_LoadBitmapFontOld (togl, this.font);
  /*fprintf(stderr, "f->base==%d\n", f->base);*/
  f->tkfont = this.tkfont;
  f->height = this.height;
  f->width  = this.width;
  f->prev   = xcFontPtr;

  xcFontPtr = f;
}


static struct xcToglFont *xcFont_find(char *font) {
  struct xcToglFont *f = xcFontPtr;
  while (f) {
    /*fprintf(stderr, "f: %s %s\n", font, f->font);*/
    if (strcmp(font,f->font) == 0) return f;
    f = f->prev;
  }
  return (struct xcToglFont*) NULL;
}



#ifdef WIN32

/*
 * The following structure represents Windows' implementation of a font.
 */

typedef struct {
  Tk_Font font;		/* Stuff used by generic font package.  Must
			 * be first in structure. */
  HFONT hFont;		/* Windows information about font. */
  HWND hwnd;			/* Toplevel window of application that owns
				 * this font, used for getting HDC. */
  int widths[256];		/* Widths of first 256 chars in this font. */
} WinFont;

#include "togl_struct.h"
#endif /* WIN32 */


static int 
_assignFont(int font_found, struct Togl *togl, const char *fontString, FontType thisFont) 
{
  struct xcToglFont *_font;

#ifndef WIN32
  /* X11 */
  XFontStruct *fontinfo;
#else
  /* WIN32 */
  WinFont    *winfont;
  TEXTMETRIC tm;
#endif

  if (!font_found) {
    Tcl_Interp  *interp = Togl_Interp(togl);
    /* 
       use the Tk-wrapper to prevent segmentation-fault of XLoadQueryFont 
       when font-name is invalid 
    */
    Tk_Font tkfont;

    tkfont = Tk_GetFont(interp, Togl_TkWin(togl), fontString);
    /*fprintf(stderr,"*** _assignFont: tkfont == %d\n", tkfont);*/
    if (!tkfont) {
      return TCL_ERROR;
    }
    /*_font = Tk_NameOfFont (tkfont);*/
#ifndef WIN32
    /* X11 */
    fontinfo = (XFontStruct *) XLoadQueryFont( Tk_Display(Togl_TkWin(togl)), fontString );
    if (!fontinfo) {
      char *result = Tcl_Alloc (sizeof(char) * 3);      
      sprintf(result, "-1");
      Tk_FreeFont (tkfont);
      Tcl_SetResult(interp, result, TCL_DYNAMIC);
      return XC_FONT_TKYES_XNO;
    }    

    thisFont.height = fontinfo->max_bounds.ascent - fontinfo->max_bounds.descent;
    thisFont.width  = fontinfo->max_bounds.rbearing - fontinfo->max_bounds.lbearing;
#else
    /* WIN32 */

    /*HFONT      oldFont;*/
    /*oldFont = SelectObject(togl->tglGLHdc, winfont->hFont);*/

    winfont = (WinFont*) tkfont;
    if (!winfont) {
      return 0;
    }

    SelectObject(togl->tglGLHdc, winfont->hFont);
    GetTextMetrics(togl->tglGLHdc, &tm);
    /*fontString      = Tk_NameOfFont(tkfont);
      fprintf(stderr,"FONT-STRING: %s\n", fontString);*/
    thisFont.height = tm.tmHeight - tm.tmInternalLeading;
    thisFont.width  = tm.tmAveCharWidth;
#endif

    thisFont.font = (char*) xcMalloc(sizeof(char) * (strlen(fontString)+1));
    strcpy(thisFont.font, fontString);
    thisFont.tkfont = tkfont;
  }

  /* now find if the font was already loaded */

  _font = xcFont_find(thisFont.font);
  /*fprintf(stderr,"Font12: %d , %s\n", _font, thisFont.font); fflush(stderr);*/

  if ( !_font ) {
    /* the font was not yet loaded */
    xcFont = (struct xcToglFont *) xcMalloc ( sizeof (struct xcToglFont) );
    xcFont_addNew (xcFont, togl, thisFont);
  } else {
    /* the font was already loaded */
    xcFont = _font;
  }
  return TCL_OK;
}




/* ------------------------------------------------------------------------ */




/*
  The Togl_LoadBitmapFont() of Togl2.0 is bad as it seg-faults. Hence
  we will use the routine from Togl-1.7 which works and rename it to
  Togl_LoadBitmapFontOld. But to this end we need the definition of
  the struct Togl.
*/
#include "GL/glx.h"

struct Togl_PackageGlobals
{
    Tk_OptionTable optionTable; /* Used to parse options */
    Togl   *toglHead;           /* Head of linked list of all Togl widgets */
    int     nextContextTag;     /* Used to assign similar context tags */
};
typedef struct Togl_PackageGlobals Togl_PackageGlobals;

struct Togl
{
    Togl   *Next;               /* next in linked list */

#if defined(TOGL_WGL)
    HDC     tglGLHdc;           /* Device context of device that OpenGL calls
                                 * will be drawn on */
    HGLRC   tglGLHglrc;         /* OpenGL rendering context to be made current */
    int     CiColormapSize;     /* (Maximum) size of colormap in color index
                                 * mode */
#  ifdef STEREO_I_H
    StereoI *pStereoI;
#  endif
#elif defined(TOGL_X11)
    GLXContext GlCtx;           /* Normal planes GLX context */
#elif defined(TOGL_AGL)
    AGLContext aglCtx;
#endif
    int     contextTag;         /* all contexts with same tag share display
                                 * lists */

    XVisualInfo *VisInfo;       /* Visual info of the current */

    Display *display;           /* X's token for the window's display. */
    Tk_Window TkWin;            /* Tk window structure */
    Tcl_Interp *Interp;         /* Tcl interpreter */
    Tcl_Command widgetCmd;      /* Token for togl's widget command */
    Togl_PackageGlobals *tpg;   /* Used to access globals */
#ifndef NO_TK_CURSOR
    Tk_Cursor Cursor;           /* The widget's cursor */
#endif
    int     Width, Height;      /* Dimensions of window */
    int     SetGrid;            /* positive is grid size for window manager */
    int     TimerInterval;      /* Time interval for timer in milliseconds */
    Tcl_TimerToken timerHandler;        /* Token for togl's timer handler */
    Bool    RgbaFlag;           /* configuration flags (ala GLX parameters) */
    int     RgbaRed;
    int     RgbaGreen;
    int     RgbaBlue;
    Bool    DoubleFlag;
    Bool    DepthFlag;
    int     DepthSize;
    Bool    AccumFlag;
    int     AccumRed;
    int     AccumGreen;
    int     AccumBlue;
    int     AccumAlpha;
    Bool    AlphaFlag;
    int     AlphaSize;
    Bool    StencilFlag;
    int     StencilSize;
    Bool    PrivateCmapFlag;
    Bool    OverlayFlag;
    int     Stereo;
    double  EyeSeparation;
    double  Convergence;
    int     AuxNumber;
    Bool    Indirect;
    int     PixelFormat;
    int     SwapInterval;
    const char *ShareList;      /* name (ident) of Togl to share dlists with */
    const char *ShareContext;   /* name (ident) to share OpenGL context with */

    const char *Ident;          /* User's identification string */
    ClientData Client_Data;     /* Pointer to user data */

    Bool    UpdatePending;      /* Should normal planes be redrawn? */

    Tcl_Obj *CreateProc;        /* Callback when widget is realized */
    Tcl_Obj *DisplayProc;       /* Callback when widget is redrawn */
    Tcl_Obj *ReshapeProc;       /* Callback when window size changes */
    Tcl_Obj *DestroyProc;       /* Callback when widget is destroyed */
    Tcl_Obj *TimerProc;         /* Callback when widget is idle */

    /* Overlay stuff */
#if defined(TOGL_X11)
    GLXContext OverlayCtx;      /* Overlay planes OpenGL context */
#elif defined(TOGL_WGL)
    HGLRC   tglGLOverlayHglrc;
#endif

    Window  OverlayWindow;      /* The overlay window, or 0 */
    Tcl_Obj *OverlayDisplayProc;        /* Overlay redraw proc */
    Bool    OverlayUpdatePending;       /* Should overlay be redrawn? */
    Colormap OverlayCmap;       /* colormap for overlay is created */
    int     OverlayTransparentPixel;    /* transparent pixel */
    Bool    OverlayIsMapped;

    GLfloat *EpsRedMap;         /* Index2RGB Maps for Color index modes */
    GLfloat *EpsGreenMap;
    GLfloat *EpsBlueMap;
    GLint   EpsMapSize;         /* = Number of indices in our Togl */
    int     currentStereoBuffer;
#ifdef HAVE_AUTOSTEREO
    int     as_initialized;     /* for autostereo package */
    ASHandle ash;               /* for autostereo package */
#endif
    int     badWindow;          /* true when Togl_CreateWindow fails */
};

#  if defined(TOGL_WGL)
#    include "tkWinInt.h"
#    include "tkFont.h"

/* 
 * The following structure represents Windows' implementation of a font.
 */

typedef struct WinFont
{
    TkFont  font;               /* Stuff used by generic font package.  Must be
                                 * first in structure. */
    HFONT   hFont;              /* Windows information about font. */
    HWND    hwnd;               /* Toplevel window of application that owns
                                 * this font, used for getting HDC. */
    int     widths[256];        /* Widths of first 256 chars in this font. */
} WinFont;
#  endif /* TOGL_WGL */


#  define MAX_FONTS 1000
static GLuint ListBase[MAX_FONTS];
static GLuint ListCount[MAX_FONTS];

/* 
 * "Standard" fonts which can be specified to Togl_LoadBitmapFontOld()
 */
#undef TOGL_BITMAP_8_BY_13
#undef TOGL_BITMAP_9_BY_15
#undef TOGL_BITMAP_TIMES_ROMAN_10
#undef TOGL_BITMAP_TIMES_ROMAN_24
#undef TOGL_BITMAP_HELVETICA_10
#undef TOGL_BITMAP_HELVETICA_12
#undef TOGL_BITMAP_HELVETICA_18
#define TOGL_BITMAP_8_BY_13		((char *) 1)
#define TOGL_BITMAP_9_BY_15		((char *) 2)
#define TOGL_BITMAP_TIMES_ROMAN_10	((char *) 3)
#define TOGL_BITMAP_TIMES_ROMAN_24	((char *) 4)
#define TOGL_BITMAP_HELVETICA_10	((char *) 5)
#define TOGL_BITMAP_HELVETICA_12	((char *) 6)
#define TOGL_BITMAP_HELVETICA_18	((char *) 7)
#define DEFAULT_FONTNAME	        "fixed"

/* 
 * Load the named bitmap font as a sequence of bitmaps in a display list.
 * fontname may be one of the predefined fonts like TOGL_BITMAP_8_BY_13
 * or an X font name, or a Windows font name, etc.
 */
GLuint
Togl_LoadBitmapFontOld(const Togl *togl, const char *fontname)
{
    static Bool FirstTime = True;

#  if defined(TOGL_X11)
    XFontStruct *fontinfo;
#  elif defined(TOGL_WGL)
    WinFont *winfont;
    HFONT   oldFont;
    TEXTMETRIC tm;
#  endif
    /* TOGL_X11 */
    int     first, last, count;
    GLuint  fontbase;
    const char *name;

    /* Initialize the ListBase and ListCount arrays */
    if (FirstTime) {
        int     i;

        for (i = 0; i < MAX_FONTS; i++) {
            ListBase[i] = ListCount[i] = 0;
        }
        FirstTime = False;
    }

    /* 
     * This method of selecting X fonts according to a TOGL_ font name
     * is a kludge.  To be fixed when I find time...
     */
    if (fontname == TOGL_BITMAP_8_BY_13) {
        name = "8x13";
    } else if (fontname == TOGL_BITMAP_9_BY_15) {
        name = "9x15";
    } else if (fontname == TOGL_BITMAP_TIMES_ROMAN_10) {
        name = "-adobe-times-medium-r-normal--10-100-75-75-p-54-iso8859-1";
    } else if (fontname == TOGL_BITMAP_TIMES_ROMAN_24) {
        name = "-adobe-times-medium-r-normal--24-240-75-75-p-124-iso8859-1";
    } else if (fontname == TOGL_BITMAP_HELVETICA_10) {
        name = "-adobe-helvetica-medium-r-normal--10-100-75-75-p-57-iso8859-1";
    } else if (fontname == TOGL_BITMAP_HELVETICA_12) {
        name = "-adobe-helvetica-medium-r-normal--12-120-75-75-p-67-iso8859-1";
    } else if (fontname == TOGL_BITMAP_HELVETICA_18) {
        name = "-adobe-helvetica-medium-r-normal--18-180-75-75-p-98-iso8859-1";
    } else if (!fontname) {
        name = DEFAULT_FONTNAME;
    } else {
        name = (const char *) fontname;
    }

    assert(name);

#  if defined(TOGL_X11)
    fontinfo = (XFontStruct *) XLoadQueryFont(Tk_Display(togl->TkWin), name);
    if (!fontinfo) {
        return 0;
    }
    first = fontinfo->min_char_or_byte2;
    last = fontinfo->max_char_or_byte2;
#  elif defined(TOGL_WGL)
    winfont = (WinFont *) Tk_GetFont(togl->Interp, togl->TkWin, name);
    if (!winfont) {
        return 0;
    }
    oldFont = SelectObject(togl->tglGLHdc, winfont->hFont);
    GetTextMetrics(togl->tglGLHdc, &tm);
    first = tm.tmFirstChar;
    last = tm.tmLastChar;
#  elif defined(TOGL_AGL_CLASSIC) || defined(TOGL_AGL)
    first = 10;                 /* don't know how to determine font range on
                                 * Mac... */
    last = 127;
#  endif
    /* TOGL_X11 */

    count = last - first + 1;
    fontbase = glGenLists((GLuint) (last + 1));
    if (fontbase == 0) {
#  ifdef TOGL_WGL
        SelectObject(togl->tglGLHdc, oldFont);
        Tk_FreeFont((Tk_Font) winfont);
#  endif
        /* TOGL_WGL */
        return 0;
    }
#  if defined(TOGL_WGL)
    wglUseFontBitmaps(togl->tglGLHdc, first, count, (int) fontbase + first);
    SelectObject(togl->tglGLHdc, oldFont);
    Tk_FreeFont((Tk_Font) winfont);
#  elif defined(TOGL_X11)
    glXUseXFont(fontinfo->fid, first, count, (int) fontbase + first);
#  elif defined(TOGL_AGL_CLASSIC) || defined(TOGL_AGL)
    aglUseFont(togl->aglCtx, 1, 0, 14,  /* for now, only app font, regular
                                         * 14-point */
            10, 118, fontbase + first);
#  endif

    /* Record the list base and number of display lists for
     * Togl_UnloadBitmapFont(). */
    {
        int     i;

        for (i = 0; i < MAX_FONTS; i++) {
            if (ListBase[i] == 0) {
                ListBase[i] = fontbase;
                ListCount[i] = last + 1;
                break;
            }
        }
    }

    return fontbase;
}



/* 
 * Release the display lists which were generated by Togl_LoadBitmapFontOld().
 */
void
Togl_UnloadBitmapFontOld(const Togl *togl, GLuint fontbase)
{
    int     i;

    (void) togl;
    for (i = 0; i < MAX_FONTS; i++) {
        if (ListBase[i] == fontbase) {
            glDeleteLists(ListBase[i], ListCount[i]);
            ListBase[i] = ListCount[i] = 0;
            return;
        }
    }
}


