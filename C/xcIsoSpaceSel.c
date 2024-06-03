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
* Source: $XCRYSDEN_TOPDIR/C/xcIsoSpaceSel.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <togl.h>
#include <tk.h>
#include <GL/gl.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xcGLparam.h"
#include "struct.h"
#include "isosurf.h"
#include "primitives.h"
#include "vector.h"
#include "xcfunc.h"

extern Options3D is;
static float poi[4][3], poi_nml[3];
static CellCage cage;
static double cage_orig[3];

extern int normalizepvfv(float *vec);

static int Paralle2D( GetComOption points, GetComOption mar, 
		      boolean rectangular );
static int Center2D( GetComOption center, GetComOption nml, 
		     float size, float rot );
static void GetFractCoor(double ivec[][4], double frcoor[4][3]);

/*****************************************************************************
 * IN THIS FILE THE FOLLOWING COMMANDS ARE IMPLEMENTED:: 
 * ---------------------------------------------------
 * 
 *     xc_isospacesel
 *****************************************************************************/


/******************************************************************************
 * XC_IsoSpaceSelCmd --> inplementation of 'xc_isospacesel' custom Tcl command 
 * --------------------- 
 * Usage:  xc_isospacesel <toglName> init
 *                              or
 *                                  paralle2D -points {p1 p2 p3} \; px={x y z}
 *                                            -margins {ab cd ad bc} \
 *                                            -rectangular 0|1
 *                              or
 *                                  center2D  -center {x y z} \
 *                                            -normal {x y z} \
 *                                            -size   size \
 *                                            -rotate rot
 *                              or
 *                                  cell3D -ctype prim|conv \
 *                                         -margins {a b c a* b* c*}
 *                              or
 *                                  paralle3D -points {p1 p2 p3} \; px={x y z}
 *                                            -margins {ab cd ad bc down up} \
 *                                            -rectangular 0|1
 *
 *                              or  center3D  -center {x y z} \
 *                                            -normal {x y z} \
 *                                            -margins {down up} \
 *                                            -size   size \
 *                                            -rotate rot
 *
 *                                  fractcoor -vtype prim|conv
 *                              or 
 *                                  clear
 */
int
XC_IsoSpaceSelCmd(ClientData clientData, Tcl_Interp *interp,
		  int argc, const char *argv[])
{ 
  struct Togl *togl;
  register int i, j;
  double size = 1.0, rot = 0.0, (*cvec)[4];
  static float rotate = 0.0;
  char *result = Tcl_Alloc( sizeof(char) * 512); /* maximum lenght of result 
						    string is 512 characters */

  boolean  rectangular  = 0;
  boolean  paralle2D_is = 0;
  boolean  center2D_is  = 0;
  boolean  cell3D_is    = 0;
  GetComOption points   = { 0, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0} };
  GetComOption margins  = { 0, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0} };
  GetComOption center   = { 0, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0} };
  GetComOption normal   = { 0, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0} };

  /* t.k */
  sprintf(result," ");

  /* HERE */
  if ( argc < 3) {
    Tcl_SetResult(interp, "Usage:   xc_isospacesel <toglName> init\n or\n paralle2D -points {x1 y1 z1  x2 y2 z2  x3 y3 z3} -margins {ab cd ad bc} -rectangular 0|1\n or\n center2D  -center {x y z} -normal {x y z} -size size -rotate rot\n or\n clear", TCL_STATIC);
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

  if ( strcmp(argv[2], "init") == 0 ) {
    rotate = 0.0;
    return TCL_OK;
  }
  else if ( strcmp(argv[2], "paralle2D") == 0 ) {
    paralle2D_is = 1;
    for (i=3; i<argc; i += 2) {
      if ( strcmp(argv[i], "-points") == 0 ) {
	/* 
	   syntax: -points {x1 y1 z1   x2 y2 z2   x3 y3 z3   x4 y4 z4}
	*/
	if ( !xcSplitList( XC_GET_FLOAT12, interp, &argv[i+1], &points ) ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "parse error; when parsing -points {%s}; should be -points {x1 y1 z1   x2 y2 z2   x3 y3 z3   x4 y4 z4}", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }
      else if ( strcmp(argv[i], "-margins") == 0 ) {
	/*
	  syntax: -margins {ab cd ad bc}
	*/
	if ( !xcSplitList( XC_GET_FLOAT4, interp, &argv[i+1], &margins ) ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "parse error; when parsing -margins {%s}; should be -margins {ab cd ad bc}", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;	
	} 
      }
      else if ( strncmp(argv[i], "-rect", 5) == 0 ) {
	/*
	  syntax: -rectangular 0|1
	*/
	if ( Tcl_GetInt( interp, argv[i+1], &rectangular ) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "wanted integer, but got %s", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;		  
	} 
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss), "unknown options %s for xc_isospacesel <mesaWin> paralle2D, should be one of -points, -margins or -rectangular", argv[i]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
  }
  else if ( strcmp(argv[2], "center2D") == 0 ) {
    center2D_is = 1;
    for (i=3; i<argc; i += 2) {    
      if ( strcmp(argv[i], "-center") == 0 ) {
	/*
	  syntax: -center {x y z}
	*/
	if ( !xcSplitList( XC_GET_XYZ, interp, &argv[i+1], &center ) ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "parse error; when parsing -center {%s}; should be -center {x y z}", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}	
      }
      else if ( strcmp(argv[i], "-normal") == 0 ) {
	/*
	  syntax: -normal {x y z}
	*/
	if ( !xcSplitList( XC_GET_XYZ, interp, &argv[i+1], &normal ) ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "parse error; when parsing -normal {%s}; should be -normal {x y z}", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}		
      }
      else if ( strcmp(argv[i], "-size") == 0 ) {
	if ( Tcl_GetDouble( interp, argv[i+1], &size ) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "wanted double, but got %s", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;		  
	}
      }
      else if ( strncmp(argv[i], "-rot", 4) == 0 ) {
	if ( Tcl_GetDouble( interp, argv[i+1], &rot ) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "wanted double, but got %s", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss), "unknown option %s for xc_isospacesel <mesaWin> center2D, should be one of -center, -normal, -size or -rotate", argv[i]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
  }
  else if ( strcmp(argv[2], "cell3D") == 0 ) {
    cell3D_is = 1;
    for (i=3; i<argc; i += 2) {  
      if ( strcmp(argv[i], "-ctype") == 0 ) {
	if ( strcmp(argv[i+1], "prim") == 0 ) {
	  cvec = vec.prim;
	}
	else if ( strcmp(argv[i+1], "conv") == 0 ) {
	  cvec = vec.conv;
	} else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), 
		   "unknown cell type %s, must be prim or conv", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }
      else if ( strcmp(argv[i], "-margins") == 0 ) {
	if ( !xcSplitList( XC_GET_FLOAT6, interp, &argv[i+1], &margins ) ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "parse error; when parsing -margins {%s}; should be -margins {a b c a* b* c*}", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}	
      }
    }
  }
  else if ( strcmp(argv[2], "fractcoor") == 0 ) {
    double frcoor[4][3];
    /* use vec.recprim or vec.recconv */    
    /* so far argv[3] should be -vtype */
    if ( argc != 5 ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "Usage: xc_isospacesel <toglName> fractcoor -vtype prim|conv");
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( strcmp(argv[4], "prim") == 0 ) {
      /* according to primitive vectors */
      GetFractCoor(vec.recprim, frcoor);
    } else if ( strcmp(argv[4], "conv") == 0 ) {
      /* according to conventional vectors */
      GetFractCoor(vec.recconv, frcoor);
    } else {
      char rss[1024];
      snprintf(rss, sizeof(rss), "unknown vector type %s, must be prim or conv", argv[4]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    /* get the results */
    sprintf(result, "%f %f %f   %f %f %f   %f %f %f   %f %f %f",
	    frcoor[0][0], frcoor[0][1], frcoor[0][2],  
	    frcoor[1][0], frcoor[1][1], frcoor[1][2],
	    frcoor[2][0], frcoor[2][1], frcoor[2][2],  
	    frcoor[3][0], frcoor[3][1], frcoor[3][2]);
    Tcl_SetResult(interp, result, TCL_DYNAMIC);
    return TCL_OK;
  }
  else if ( strcmp(argv[2], "clear") == 0 ) {
    rotate = 0.0;
    VPf.isospacesel2D = GL_FALSE;
    VPf.isospacesel3D = GL_FALSE;
  }
  else {
    char rss[1024];
    snprintf(rss, sizeof(rss), "unknown type \"%s\" of xc_isospacesel <mesaWin> command, should be one of init, paralle2D, center2D or clear", argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }


  if (paralle2D_is) {
    /* 
       paralle2D mode
    */
    if (!Paralle2D( points, margins, rectangular )) {
      Tcl_SetResult(interp, "Usage: xc_iso isoplane <type> <colorBase> <scalefunc> <what> ?<islide>?", TCL_STATIC);
      return TCL_ERROR;
    }
    VPf.isospacesel2D = GL_TRUE;
    VPf.isospacesel3D = GL_FALSE;
  } 
  else if (center2D_is) {
    /*
      center2D mode
    */
    rotate += rot;
    if (!Center2D( center, normal, size, rotate )) {
      Tcl_SetResult(interp, "zero vector was submitted for normal, when executing xc_isospacesel <mesaWin> center2D -normal {x y z}", TCL_STATIC);	
      return TCL_ERROR;
    }
    VPf.isospacesel2D = GL_TRUE;
    VPf.isospacesel3D = GL_FALSE;
  } else if (cell3D_is) {
    /*
      cell3D mode
    */
    double mat[4][4];

    for (i=0; i<4; i++)
      for (j=0; j<4; j++)
	mat[i][j] = 0.0;

    VPf.isospacesel3D = GL_TRUE;
    VPf.isospacesel2D = GL_FALSE;
    for (i=0; i<xcr.dim; i++) {
      cage_orig[i] = -margins.vec[0]*cvec[0][i] - margins.vec[1]*cvec[1][i] -
	margins.vec[2]*cvec[2][i];
      for (j=0; j<xcr.dim; j++) 
	mat[j][i] =  (1.0 + margins.vec[j] + margins.vec[j+3]) * cvec[j][i];
    }
    if ( xcr.dim < 3 ) {
      cage_orig[0] = -margins.vec[0]*cvec[0][0] - margins.vec[1]*cvec[1][0];
      cage_orig[1] = -margins.vec[0]*cvec[0][1] - margins.vec[1]*cvec[1][1];
      cage_orig[2] = -margins.vec[2];
      mat[2][2]    =  margins.vec[2] + margins.vec[5];
    } 
    if ( xcr.dim < 2 ) {
      cage_orig[0] = -margins.vec[0]*cvec[0][0];
      cage_orig[1] = -margins.vec[1];
      cage_orig[2] = -margins.vec[2];
      mat[1][1]    =  margins.vec[1] + margins.vec[4];
    } 
    if ( xcr.dim < 1 ) {
      for(i=0; i<3; i++)
	cage_orig[i] = -margins.vec[i];
      mat[0][0] =  margins.vec[0] + margins.vec[3];
    }

    SetUnitCellCage( mat, &cage );
    MVf.isospaceselsize = 
      DetermineParapipedSized(mat[0], mat[1], mat[2], cage_orig);

    /* we must return back the origin & three vectors */
    sprintf(result, "%f %f %f   %f %f %f   %f %f %f   %f %f %f   ",
	    cage_orig[0], cage_orig[1], cage_orig[2],  
	    mat[0][0], mat[0][1], mat[0][2],
	    mat[1][0], mat[1][1], mat[1][2],
	    mat[2][0], mat[2][1], mat[2][2]);
  }
  /* 
     determine the parallelogram size
  */
  if ( paralle2D_is || center2D_is ) {
    float o[3], v0[3], v1[3];
    /* 0--3
       |  |
       1--2    vec0 = p[0] - p[1], vec[1] = p[2] - p[1] */
    o[0] = -poi[1][0] + mx;
    o[1] = -poi[1][1] + my;
    o[2] = -poi[1][2] + mz;
    for (i=0; i<3; i++) {
      v0[i] = poi[0][i] - poi[1][i];
      v1[i] = poi[2][i] - poi[1][i];
    }
    MVf.isospaceselsize = DetermineParalleSize( v0, v1, o );
    /* we must return back the edges of parallelogram/parallelepiped */
    sprintf(result, "%f %f %f   %f %f %f   %f %f %f   %f %f %f   ",
	    poi[0][0], poi[0][1], poi[0][2],  poi[1][0], poi[1][1], poi[1][2],
	    poi[2][0], poi[2][1], poi[2][2],  poi[3][0], poi[3][1], poi[3][2]);
  }

  /* because isospacesel may be greater than structer, 
     take care of projection */
  if ( VPf.isospacesel2D || VPf.isospacesel3D ) {
    if (is.stickmode && !is.ballmode) xcMakeProjection3D("sticks");
    if (is.ballmode) xcMakeProjection3D("balls");
    if (is.spacefillmode) xcMakeProjection3D("space");
  }

  Togl_PostRedisplay(togl);
  Tcl_SetResult(interp, result, TCL_DYNAMIC);    
  return TCL_OK;
}


static int
Paralle2D( GetComOption points, GetComOption mar, boolean rectangular )
{
  /*
    points.vec[0..2]  point A
    points.vec[3..5]  point B
    points.vec[6..8]  point C
    points.vec[9..11] point D

    mar.vec[0]    margin AB
    mar.vec[1]    margin CD
    mar.vec[2]    margin AD
    mar.vec[3]    margin BC
  */
  register int i, j;
  int ii[2];
  float det, k;
  float *p[4], bc[3], ba[3], nml[3], bc_nml[3], d[2][3];

  p[0] = points.vec;
  p[1] = points.vec + 3;
  p[2] = points.vec + 6;
  p[3] = points.vec + 9;
  for (i=0; i<3; i++) {    
    bc[i] = points.vec[6+i] - points.vec[3+i];
    ba[i] = points.vec[i]   - points.vec[3+i];
  }

  VecProduct3f( bc, ba, nml );
  if ( distfv( nml ) < 1.0e-6 ) {
    /* linear dependence */
    return XC_ERROR;
  }
  VecProduct3f( bc, nml, bc_nml );

  if ( bc_nml[0]*ba[0] + bc_nml[1]*ba[1] + bc_nml[2]*ba[2] < 0 )
    RevertVectorfv( bc_nml );
  /* normalize bc_nml */
  normalizepvfv( bc_nml );

  /* was "-rectangular 1" option specified ???? */
  if ( rectangular ) {
    /* check if it is already rectangular ???? */
    int aa;
    /* 
       if AC < BD then move A->A' && C->C'
       else            move B->B' && D->D'
    */
    if ( dist3f( p[0][0] - p[2][0], p[0][1] - p[2][1], p[0][2] - p[2][2] ) <
	 dist3f( p[1][0] - p[3][0], p[1][1] - p[3][1], p[1][2] - p[3][2] ) ) {
      /* move A -> A' */
      aa = 1;
    } else {
      /* move B -> B' */
      aa = 0;
    }

    det = 0.0;
    for (j=0; j<3; j++) {
      if (j==1) {
	ii[0] = 0;
	ii[1] = 1;
      } else if (j==2) {
	ii[0] = 0;
	ii[1] = 2;
      } else {
	ii[0] = 1;
	ii[1] = 2;
      }
      for (i=0; i<2; i++) {
	d[i][0] = bc[ii[i]];
	d[i][1] = -bc_nml[ii[i]];
	if ( aa ) {
	  d[i][2] = p[1][ii[i]] - p[0][ii[i]];
	} else {
	  d[i][2] = p[0][ii[i]] - p[1][ii[i]];
	}
      }
      det = CramerRule2x2Detf( d, 0 );
      if ( fabs(det) > 1.0e-6 ) break;
    }
    if ( fabs(det) < 1.0e-6 ) {
      /* linear dependence */
      return XC_ERROR;
    }
    k = CramerRule2x2Detf( d, 1 ) / det;
    for (i=0; i<3; i++) {
      if ( aa ) {
	p[0][i] += k * bc[i];
	p[2][i] -= k * bc[i]; /* C point is moved just in opposite direction */
      } else {
	p[1][i] += k * bc[i];
	p[3][i] -= k * bc[i];
      }
    }
  }

  /* take margins into account */
  /* margins are specified in fractional units */
  if ( rectangular ) {
    /* assing new, rectangular bc & ba vectors */
    for (i=0; i<3; i++) {
      bc[i] = p[2][i] - p[1][i];
      ba[i] = p[0][i] - p[1][i];
    }
  }
  for (i=0; i<3; i++) {
    p[0][i] += -mar.vec[0] * bc[i] + mar.vec[2] * ba[i];
    p[1][i] += -mar.vec[0] * bc[i] - mar.vec[3] * ba[i];
    p[2][i] +=  mar.vec[1] * bc[i] - mar.vec[3] * ba[i];
    p[3][i] +=  mar.vec[1] * bc[i] + mar.vec[2] * ba[i];
  }

  /* assign parallelogram's normal */
  for (i=0; i<3; i++) {
    poi_nml[i] = nml[i];
    for (j=0; j<4; j++)
      poi[j][i] = p[j][i];
  }
  return XC_OK;
}


static int
Center2D( GetComOption center, GetComOption nml, float size, float rot )
{
  register int i, j;
  float v[3], vnml[3], p[4][3];


  /* check if zero lenght normal was submited */
  if ( fabs(nml.vec[0]) < 1.0e-3 &&
       fabs(nml.vec[1]) < 1.0e-3 &&
       fabs(nml.vec[2]) < 1.0e-3 ) return XC_ERROR; /* zero-lenght normal */

  if ( (nml.vec[0] * nml.vec[1] * nml.vec[2]) < 1.0e-4 ) {
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    if ( fabs(nml.vec[2]) < 1.e-3 ) 
      v[2] = 1.0;
    if ( fabs(nml.vec[1]) < 1.e-3 )
      v[1] = 1.0;
    if ( fabs(nml.vec[0]) < 1.e-3 )
      v[0] = 1.0;
  } else {
    v[0] = 1.0;
    v[1] = 1.0;
    v[2] = - (nml.vec[0]*v[0] + nml.vec[1]*v[1]) / nml.vec[2];
  }

  VecProduct3f( v, nml.vec, vnml );
  normalizepvfv( v );
  normalizepvfv( vnml );
  if ( fabs(v[0]*vnml[0]+v[1]*vnml[1]+v[2]*vnml[2]) > 1.e-6 )
    xcError("VecProduct3f is false!!!");

  /* let the size be half of square edge, not square diagonal */
  size *= 2.0 / sqrt(3.0);
  for (i=0; i<3; i++) {
    v[i]    *= size;
    vnml[i] *= size;
    p[0][i] = center.vec[i] + vnml[i];
    p[1][i] = center.vec[i] + v[i];
    p[2][i] = center.vec[i] - vnml[i];
    p[3][i] = center.vec[i] - v[i];
    printf("point:: i=%d, %f %f %f %f\n", i, p[0][i], p[1][i], p[2][i], p[3][i]);
    fflush(stdout);
  }

  /* now rotate the points around normal */
  rot /= RAD2DEG;
  if ( rot != 0 )
    for (i=0; i<4; i++)
      xcRotate3fv( rot, nml.vec, center.vec, p[i] );

  /* assign poi & poin_nml */
  for (i=0; i<3; i++) {
    poi_nml[i] = nml.vec[i];
    for (j=0; j<4; j++)
      poi[j][i] = p[j][i];
  }

  return XC_OK;
}


void
IsoSpaceSel_Parallelogram(void)
{
  GLint     shade_model;
  GLboolean two_side, cull_face;
  GLboolean switch_shademodel = GL_FALSE;
  GLboolean switch_lightmodel = GL_FALSE;
  GLboolean switch_cullface   = GL_FALSE;

  /*glPushAttrib();*/

  /* TRANSPARENCY */
  LoadBlendfunc_And_Frontface();
  glEnable( GL_BLEND );
  glDepthMask( GL_FALSE );
  glBlendFunc( blend_cellcage.sfunc, blend_cellcage.dfunc );

  /* SHADEMODEL */
  glGetIntegerv( GL_SHADE_MODEL, &shade_model );
  if ( shade_model != GL_FLAT ) { 
    switch_shademodel = GL_TRUE;
    glShadeModel( GL_FLAT );
  } 

  /* LIGHTMODEL */
  glGetBooleanv( GL_LIGHT_MODEL_TWO_SIDE, &two_side );
  if ( !two_side ) {
    switch_lightmodel = GL_TRUE;
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );
  }

  /* CULL FACE */
  glGetBooleanv( GL_CULL_FACE, &cull_face );
  if ( cull_face ) {
    switch_cullface = GL_TRUE;
    glDisable( GL_CULL_FACE );
  }

  LoadCageOrVecMaterial( GLPAR_ISOPARALLE );

  glPushMatrix();
  glTranslated( -mx, -my, -mz );
  xcParallelogram( XCPRIM_SOLID_AND_BORDER, poi, poi_nml );
  glPopMatrix();

  /* get back the previous GL state */
  glDisable(GL_BLEND); 
  glDepthMask(GL_TRUE);  

  if ( switch_shademodel ) {
    glShadeModel( GL_SMOOTH );
  }

  if ( switch_lightmodel ) {
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE );
  }

  if ( switch_cullface ) {
    glEnable( GL_CULL_FACE );
  }

  LoadStructMaterial();
  /*glPopMatrix();*/
}


void
IsoSpaceSel_3D(void)
{
  GLint     shade_model;
  GLboolean two_side, cull_face;
  GLboolean switch_shademodel = GL_FALSE;
  GLboolean switch_lightmodel = GL_FALSE;
  GLboolean switch_cullface   = GL_FALSE;

  /* TRANSPARENCY */
  LoadBlendfunc_And_Frontface();
  glEnable( GL_BLEND );
  glDepthMask( GL_FALSE );
  glBlendFunc( blend_cellcage.sfunc, blend_cellcage.dfunc );

  /* SHADEMODEL */
  glGetIntegerv( GL_SHADE_MODEL, &shade_model );
  if ( shade_model != GL_FLAT ) { 
    switch_shademodel = GL_TRUE;
    glShadeModel( GL_FLAT );
  } 

  /* LIGHTMODEL */
  glGetBooleanv( GL_LIGHT_MODEL_TWO_SIDE, &two_side );
  if ( !two_side ) {
    switch_lightmodel = GL_TRUE;
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );
  }

  /* CULL FACE */
  glGetBooleanv( GL_CULL_FACE, &cull_face );
  if ( cull_face ) {
    switch_cullface = GL_TRUE;
    glDisable( GL_CULL_FACE );
  }

  LoadCageOrVecMaterial( GLPAR_ISOPARALLE );

  glPushMatrix();
  glTranslated( -mx + cage_orig[0], -my + cage_orig[1], -mz + cage_orig[2]);
  glLineWidth( 2.0 );
  xcWireCage( cage );
  xcSolidCage( cage );
  glPopMatrix();

  /* get back the previous GL state */
  glDisable(GL_BLEND); 
  glDepthMask(GL_TRUE);  

  if ( switch_shademodel ) {
    glShadeModel( GL_SMOOTH );
  }

  if ( switch_lightmodel ) {
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE );
  }

  if ( switch_cullface ) {
    glEnable( GL_CULL_FACE );
  }

  LoadStructMaterial();
}


static void
GetFractCoor(double ivec[][4], double frcoor[4][3])
{
  int i;

  for (i=0; i<4; i++) {
    frcoor[i][0] = ivec[0][0] * poi[i][0] + ivec[0][1] * poi[i][1] +
      ivec[0][2] * poi[i][2];

    frcoor[i][1] = ivec[1][0] * poi[i][0] + ivec[1][1] * poi[i][1] +
      ivec[1][2] * poi[i][2];

    frcoor[i][2] = ivec[2][0] * poi[i][0] + ivec[2][1] * poi[i][1] +
      ivec[2][2] * poi[i][2];
  }
}
