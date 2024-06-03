struct Togl
{
    Togl   *Next;               /* next in linked list */

#if defined(TOGL_WGL)
    HDC     tglGLHdc;           /* Device context of device that OpenGL calls
                                 * will be drawn on */
    HGLRC   tglGLHglrc;         /* OpenGL rendering context to be made current */
    int     CiColormapSize;     /* (Maximum) size of colormap in color index
                                 * mode */
#elif defined(TOGL_X11)
    GLXContext GlCtx;           /* Normal planes GLX context */
#elif defined(TOGL_AGL_CLASSIC) || defined(TOGL_AGL)
    AGLContext aglCtx;
#endif                          /* TOGL_WGL */

    Display *display;           /* X's token for the window's display. */
    Tk_Window TkWin;            /* Tk window structure */
    Tcl_Interp *Interp;         /* Tcl interpreter */
    Tcl_Command widgetCmd;      /* Token for togl's widget command */
#ifndef NO_TK_CURSOR
    Tk_Cursor Cursor;           /* The widget's cursor */
#endif
    int     Width, Height;      /* Dimensions of window */
    int     SetGrid;            /* positive is grid size for window manager */
    int     TimerInterval;      /* Time interval for timer in milliseconds */
#if (TCL_MAJOR_VERSION * 100 + TCL_MINOR_VERSION) >= 705
    Tcl_TimerToken timerHandler;        /* Token for togl's timer handler */
#else
    Tk_TimerToken timerHandler; /* Token for togl's timer handler */
#endif
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
    Bool    StereoFlag;
#ifdef __sgi
    Bool    OldStereoFlag;
#endif
    int     AuxNumber;
    Bool    Indirect;
    int     PixelFormat;
    const char *ShareList;      /* name (ident) of Togl to share dlists with */
    const char *ShareContext;   /* name (ident) to share OpenGL context with */

    const char *Ident;          /* User's identification string */
    ClientData Client_Data;     /* Pointer to user data */

    Bool    UpdatePending;      /* Should normal planes be redrawn? */

    Togl_Callback *CreateProc;  /* Callback when widget is created */
    Togl_Callback *DisplayProc; /* Callback when widget is rendered */
    Togl_Callback *ReshapeProc; /* Callback when window size changes */
    Togl_Callback *DestroyProc; /* Callback when widget is destroyed */
    Togl_Callback *TimerProc;   /* Callback when widget is idle */

    /* Overlay stuff */
#if defined(TOGL_X11)
    GLXContext OverlayCtx;      /* Overlay planes OpenGL context */
#elif defined(TOGL_WGL)
    HGLRC   tglGLOverlayHglrc;
#endif                          /* TOGL_X11 */

    Window  OverlayWindow;      /* The overlay window, or 0 */
    Togl_Callback *OverlayDisplayProc;  /* Overlay redraw proc */
    Bool    OverlayUpdatePending;       /* Should overlay be redrawn? */
    Colormap OverlayCmap;       /* colormap for overlay is created */
    int     OverlayTransparentPixel;    /* transparent pixel */
    Bool    OverlayIsMapped;

    /* for DumpToEpsFile: Added by Miguel A. de Riera Pasenau 10.01.1997 */
    XVisualInfo *VisInfo;       /* Visual info of the current */
    /* context needed for DumpToEpsFile */
    GLfloat *EpsRedMap;         /* Index2RGB Maps for Color index modes */
    GLfloat *EpsGreenMap;
    GLfloat *EpsBlueMap;
    GLint   EpsMapSize;         /* = Number of indices in our Togl */
};

