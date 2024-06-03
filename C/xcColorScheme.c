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
* Source: $XCRYSDEN_TOPDIR/C/xcColorScheme.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <togl.h>
#include <tk.h>
#include <stdio.h>
#include <string.h>
#include "struct.h"
#include "xcGLparam.h"
#include "isosurf.h"
#include "xcfunc.h"

#define COLSH_DEFAULT 0
#define COLSH_NN      1
#define COLSH_SLAB    2

#define COLSH_MONOCHROME   0
#define COLSH_RGB          1
#define COLSH_RAINBOW      2
#define COLSH_GEOGRAPHIC   3
#define COLSH_N_COLOR_BASE 4

#define COLSH_OVERRIDE  0
#define COLSH_COMBINED  1

#define COLSH_XP  0
#define COLSH_XM  1
#define COLSH_YP  2
#define COLSH_YM  3
#define COLSH_ZP  4
#define COLSH_ZM  5

struct COLOR_SCHEME {
  int whatscheme;
  int colorscheme;
  int colortype;
  GetGlParam nn_coor;
  int slab_dir;
  GetGlParam slabrange;
  float alpha; /* alpha factor for combined color */
  float dist_r; /* distance radius */
} colSh;

static BASECOLOR bcol[COLSH_N_COLOR_BASE];

extern float monocol[2][4];
extern float rainbow[6][4];
extern float rgbcol[3][4];
extern float geographic[11][4];

/* --- static function prototypes --- */
static void Colorise(void);
static void AssignColorInit(void);
static void AssignColor(GLfloat col[3], float x, int colorscheme, 
			int colortype, float oldcol[3], float alpha);

/*****************************************************************************
xc_colorscheme structure toglName (-nn {x y z}|-slab +x|-x|+y|-y|+z|-z
-slabrange {min max}
-r radius
-default)
-colorscheme monochrome|rgb|
rainbow|geographic
-colortype override|combined
-alpha <alpha>
*****************************************************************************/
int 
XC_ColorschemeCmd(ClientData clientData,Tcl_Interp *interp,
		  int argc, const char *argv[])
{
  static int first_time = 1;
  register int i;
  struct Togl *togl;

  /* initialization */
  if (first_time) {
    AssignColorInit();
    first_time = 0;
    colSh.slabrange.vec[0] = 0.0;
    colSh.slabrange.vec[1] = 1.0;
  }
  /* this is default */
  colSh.whatscheme = COLSH_DEFAULT;

  /* so far just "xc_colorscheme structure" is supported */
  if ( strncmp(argv[1], "struct", 6) != 0 ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),
	     "Usage: xc_colorscheme structure toglName (-nn {x y z}|-slab +x|-x|+y|-y|+z|-z | -slabrange {min max} | -default) -colorscheme monochrome|rgb|rainbow|geographic -colortype override|combined -alpha <alpha>");
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

  for(i=3; i<argc; i+=2) {
    if ( strcmp(argv[i], "-nn") == 0 ) {
      if ( !xcSplitList( XC_GET_XYZ, interp, argv + i + 1, &colSh.nn_coor ) )
	return TCL_ERROR;
      colSh.whatscheme = COLSH_NN;
    }
    else if ( strcmp(argv[i], "-slab") == 0 ) {
      colSh.whatscheme = COLSH_SLAB;     
      if ( strcmp(argv[i+1], "+x") == 0 )
	colSh.slab_dir = COLSH_XP;
      else if ( strcmp(argv[i+1], "-x") == 0 )
	colSh.slab_dir = COLSH_XM;
      else if ( strcmp(argv[i+1], "+y") == 0 )
	colSh.slab_dir = COLSH_YP;
      else if ( strcmp(argv[i+1], "-y") == 0 )
	colSh.slab_dir = COLSH_YM;
      else if ( strcmp(argv[i+1], "+z") == 0 )
	colSh.slab_dir = COLSH_ZP;
      else if ( strcmp(argv[i+1], "-z") == 0 )
	colSh.slab_dir = COLSH_ZM;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss), "unknown flag %s for -slab option, must be one of +x, -x, +y, -y, +z, -z", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else if ( strcmp(argv[i], "-slabrange") == 0 ) {
      if ( !xcSplitList( XC_GET_FLOAT2, interp, argv + i + 1, &colSh.slabrange ) )
	return TCL_ERROR;
    }

    else if ( strcmp(argv[i], "-r") == 0 ) {
      double r;
      if ( Tcl_GetDouble(interp, argv[i+1], &r) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "wanted double, but got %s", argv[i]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      } else {
	colSh.dist_r = (float) r;
      }
    }
    else if ( strcmp(argv[i], "-default" ) == 0 ) {
      colSh.whatscheme = COLSH_DEFAULT;
      break;
    }
    else if ( strcmp(argv[i], "-colorscheme" ) == 0 ) {
      if ( strncmp(argv[i+1], "monochrome", 4) == 0 )
	colSh.colorscheme = COLSH_MONOCHROME;
      else if ( strcmp(argv[i+1], "rgb") == 0 )
	colSh.colorscheme = COLSH_RGB;
      else if ( strcmp(argv[i+1], "rainbow") == 0 )
	colSh.colorscheme = COLSH_RAINBOW;
      else if ( strcmp(argv[i+1], "geographic") == 0 )
	colSh.colorscheme = COLSH_GEOGRAPHIC;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss), "unknown flag %s for -colorscheme option, must be one of monochrome, rgb, rainbow, geographic", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else if ( strcmp(argv[i], "-colortype") == 0 ) {
      if ( strcmp(argv[i+1], "override") == 0 ) 
	colSh.colortype = COLSH_OVERRIDE;
      else if ( strcmp(argv[i+1], "combined") == 0 ) 
	colSh.colortype = COLSH_COMBINED;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss), "unknown flag %s for -colortype option, must be one of override, combined", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else if ( strcmp(argv[i], "-alpha") == 0 ) {
      /* this is used as an alpha factor for combined color */
      double alpha;
      if ( Tcl_GetDouble(interp, argv[i+1], &alpha) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "wanted double, but got %s", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      colSh.alpha = (float) alpha;
    }
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss), "unknown option %s for xc_colorscheme command, must be one of -nn, -slab, -default, -colorscheme, -colortype -alpha", argv[i]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  Colorise();
  /* here insert all updated stuff */
  if (VPf.stropened) 
    ReDisplay(togl, REMAKE3D_SPACE | REMAKE3D_BALL | REMAKE3D_BOND);

  return TCL_OK;
}


void
LoadAtmCol(const int flag) 
{
  register int i, j;
  if (VPf.stropened || flag) {
    /* atcol[0] == atcol[SELCOL]; thatwhy from i=0 */
    atm.col[0][0] = atcol[0][0];
    atm.col[0][1] = atcol[0][1];
    atm.col[0][2] = atcol[0][2];
    for(i=1; i<=natoms; i++)
      for(j=0; j<3; j++)
	atm.col[i][j] = atcol[ nat[i] ][j];
  }
}


static void 
Colorise(void) 
{
  register int i;
  float d, dis = 0.0;

  if ( colSh.whatscheme == COLSH_DEFAULT )
    LoadAtmCol(0);
  else if ( colSh.whatscheme == COLSH_NN ) {
    /* first find the maximum distance */
    for(i=1; i<=natoms; i++) {
      d = dist3f( colSh.nn_coor.vec[0] - xat[i] - mx,
		  colSh.nn_coor.vec[1] - yat[i] - my,
		  colSh.nn_coor.vec[2] - zat[i] - mz );
      if ( dis < d ) dis = d;
    }

    for(i=1; i<=natoms; i++) {
      d = dist3f( colSh.nn_coor.vec[0] - xat[i] - mx,
		  colSh.nn_coor.vec[1] - yat[i] - my,
		  colSh.nn_coor.vec[2] - zat[i] - mz );      
      AssignColor(atm.col[i], d/dis, colSh.colorscheme, 
		  colSh.colortype, atcol[ nat[i] ], colSh.alpha);
    }
  }
  else if ( colSh.whatscheme == COLSH_SLAB ) {
    for(i=1; i<=natoms; i++) {
      switch (colSh.slab_dir)
	{
	case COLSH_XP:
	  AssignColor(atm.col[i], (max.x + xat[i]) / (2.0 * max.x), 
		      colSh.colorscheme, colSh.colortype, 
		      atcol[ nat[i] ], colSh.alpha);
	  break;
	case COLSH_XM:
	  AssignColor(atm.col[i], (max.x - xat[i]) / (2.0 * max.x), 
		      colSh.colorscheme, colSh.colortype, 
		      atcol[ nat[i] ], colSh.alpha);
	  break;	  
	case COLSH_YP:
	  AssignColor(atm.col[i], (max.y + yat[i]) / (2.0 * max.y), 
		      colSh.colorscheme, colSh.colortype, 
		      atcol[ nat[i] ], colSh.alpha);
	  break;
	case COLSH_YM:
	  AssignColor(atm.col[i], (max.y - yat[i]) / (2.0 * max.y), 
		      colSh.colorscheme, colSh.colortype, 
		      atcol[ nat[i] ], colSh.alpha);
	  break;
	case COLSH_ZP:
	  AssignColor(atm.col[i], (max.z + zat[i]) / (2.0 * max.z), 
		      colSh.colorscheme, colSh.colortype, 
		      atcol[ nat[i] ], colSh.alpha);
	  break;
	case COLSH_ZM:
	  AssignColor(atm.col[i], (max.z - zat[i]) / (2.0 * max.z), 
		      colSh.colorscheme, colSh.colortype, 
		      atcol[ nat[i] ], colSh.alpha);
	  break;
	}
    }
  }
}


static void
AssignColor(GLfloat col[3], float x, int colorscheme, 
	    int colortype, float oldcol[3], float alpha)
{
  float fac, f;
  register int i, j;

  f = x;

  if ( colSh.whatscheme == COLSH_SLAB ) {
    if ( x < colSh.slabrange.vec[0] ) {
      f = 0.0;
    } else if ( x > colSh.slabrange.vec[1] ) {
      f = 1.0;
    } else {
      f = (x - colSh.slabrange.vec[0]) / (colSh.slabrange.vec[1]-colSh.slabrange.vec[0]);
    }
  } else if ( colSh.whatscheme == COLSH_NN ) {
    if ( x < colSh.dist_r ) {
      f = x / colSh.dist_r;
    } else {
      f = 1.0;
    }
  }

  /*fprintf(stderr,"COLOR__SCHEME: x=%10.5f, factor=%10.5f\n", x, f);*/

  if (colorscheme > COLSH_GEOGRAPHIC) colorscheme = COLSH_GEOGRAPHIC;
  else if (colorscheme < COLSH_MONOCHROME)  colorscheme = COLSH_MONOCHROME;

  for(i=0; i<bcol[colorscheme].ncol-1; i++)
    /* remember baseintv goes from 1.0 to 0.0 */
    if ( f <= bcol[colorscheme].baseintv[i] && 
	 f >= bcol[colorscheme].baseintv[i+1] ) {
      fac = (bcol[colorscheme].baseintv[i] - f) / 
	(bcol[colorscheme].baseintv[i] - bcol[colorscheme].baseintv[i+1]);
      for(j=0; j<3; j++)
	col[j] = (1.0 - fac) * bcol[colorscheme].basecol[i][j] +
	  fac * bcol[colorscheme].basecol[i+1][j];
      break;
    }  

  /* if colortype is COLSH_COMBINED above obtained color should be mixed 
     with current atcol */
  if (colortype == COLSH_COMBINED) {
    /* correct combined equation if needed */
    col[0] = alpha * col[0] + (1.0 - alpha)*oldcol[0];
    col[1] = alpha * col[1] + (1.0 - alpha)*oldcol[1];
    col[2] = alpha * col[2] + (1.0 - alpha)*oldcol[2];
  }
}


static void 
AssignColorInit(void)
{
  register int i, j;

  /* MONOCHROME */
  bcol[COLSH_MONOCHROME].ncol = 2;
  bcol[COLSH_MONOCHROME].baseintv[0] = 1.0;
  bcol[COLSH_MONOCHROME].baseintv[1] = 0.0;
  for(i=0; i<2; i++)
    for(j=0; j<4; j++)
      bcol[COLSH_MONOCHROME].basecol[i][j] = monocol[i][j];

  /* RAINBOW */
  bcol[COLSH_RAINBOW].ncol = 6;
  for (i=0;i<6;i++) 
    bcol[COLSH_RAINBOW].baseintv[i] = 1.0 - (float) i / 5.0;
  for(i=0; i<6; i++)
    for(j=0; j<4; j++)
      bcol[COLSH_RAINBOW].basecol[i][j] = rainbow[i][j];

  /* RGB */
  bcol[COLSH_RGB].ncol = 3;
  for (i=0;i<3;i++)
    bcol[COLSH_RGB].baseintv[i] = 1.0 - (float) i / 2.0;
  for(i=0; i<3; i++)
    for(j=0; j<4; j++)
      bcol[COLSH_RGB].basecol[i][j] = rgbcol[i][j];

  /* GEOGRAPHIC */
  bcol[COLSH_GEOGRAPHIC].ncol = 11;
  for (i=0;i<11;i++)
    bcol[COLSH_GEOGRAPHIC].baseintv[i] = 1.0 - (float) i / 10.0;
  for(i=0; i<11; i++)
    for(j=0; j<4; j++)
      bcol[COLSH_GEOGRAPHIC].basecol[i][j] = geographic[i][j];
}
