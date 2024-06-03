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
* Source: $XCRYSDEN_TOPDIR/C/xcGLparam.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <togl.h>
#include <tk.h>
#include <stdlib.h>
#include "struct.h"
#include "xcGLparam.h"
#include "isosurf.h"
#include "lighting.h"
#include "xcfunc.h"


extern ISO_ATTRIB isoDisp;
extern XCfog      fog;
extern GLfloat    def_fog_color[4];

/* t.k.: this is for debuging (for xcMakeprojection3D) */
extern Options3D is;
extern int togl_exists;
static void CopyGLParVec(int code, float *mat, GetGlParam par);


/*****************************************************************************/
/* XC_SetGLparamCmd --> inplementation of 'xc_setGLparam' custom Tcl command 
 * --------------------- 
 * Usage:  xc_setGLparam material \
 *			     -what structure|(isosurf pos|neg|one front|back) \
 *			     -shininess {v} \
 *			     -specular {r g b a} \
 *			     -ambient  {r g b a} \
 *			     -diffuse  {r g b a} \
 *                           -emission {r g b a}
 *	    or
 *	    
 *	    xc_setGLparam light \
 *			     -light {v} \
 *			     -ambient {r g b a} \
 *			     -diffuse {r g b a} \
 *			     -fract_position {x y z w} \
 *			     -specular {r g b a} \
 *			     -spot_dir {x,y,z} \
 *			     -spot_exp {v} \
 *			     -spot_cutoff {v} \
 *			     -const_atten {v} \
 *			     -lin_atten {v} \
 *			     -quad_atten {v}\n\
 *
 *	    or
 *	    
 *	    xc_setGLparam lightmodel \
 *			     -two_side 0|1 \
 *                           -two_side_iso 0|1 \
 *			     -ambient {r g b a} \
 *			     -local_viewer 0|1 \
 *			     -enable_light {v}|-disable_light {v}"
 *
 *          or
 *
 *          xc_setGLparam frontface \
 *                            -what (isosurf_one | isosurf_pos | \
 *                                   isosurf_neg | colorplane)
 *                            -frontface CW|CCW
 *
 *          or
 * 
 *          xc_setGLparam isonormal \
 *                            -what (isosurf_one | isosurf_pos | \
 *                                   isosurf_neg)
 *
 *          or
 *
 *          xc_setGLparam blendfunc \
 *                             -what isosurf
 *                             -sfunc <source_function>
 *                             -dfunc <destination_function>
 *          
 *          or
 *
 *          xc_setGLparam fog -color {r g b a}
 *         
 */ 
int 
XC_SetGLparamCmd(ClientData clientData, Tcl_Interp *interp,
		 int argc, const char *argv[])
{
  int i;

  struct GLParam {
    int    what;
    int    mat_what;
    int    mat_what_isosurf;
    int    mat_what_isoside;
    int    face_what;
    int    isonormal_what;
    GLenum face_frontface;
    int    blend_what;
    GLenum blend_sfunc;
    GLenum blend_dfunc;
  } glparam;

  GetGlParam specular      = { 0 };
  GetGlParam shininess     = { 0 };
  GetGlParam ambient       = { 0 };
  GetGlParam diffuse       = { 0 };
  GetGlParam emission      = { 0 };
  GetGlParam lightN        = { 0 };
  GetGlParam fract_position= { 0 };
  GetGlParam spot_dir      = { 0 };
  GetGlParam spot_exp      = { 0 };
  GetGlParam spot_cutoff   = { 0 };
  GetGlParam const_atten   = { 0 };
  GetGlParam lin_atten     = { 0 };
  GetGlParam quad_atten    = { 0 };
  GetGlParam two_side      = { 0 };
  GetGlParam two_side_iso  = { 0 };
  GetGlParam local_viewer  = { 0 };
  GetGlParam enable_light  = { 0 };
  GetGlParam disable_light = { 0 };
  GetGlParam fogcolor      = { 0 };

  boolean is_blend_sfunc = 0;
  boolean is_blend_dfunc = 0;

  if ( (argc % 2) != 0 ) {
    Tcl_SetResult(interp, "Usage: xc_setGLparam material \
-what structure|(isosurf pos|neg|one front|back) \
-shininess {v} -specular {r g b a} -ambient {r g b a} -diffuse {r g b a} -emission {r g b a}\n\
or\n\
xc_setGLparam light -light {v} \
-ambient {r g b a} -diffuse {r g b a} -fract_position {x y z w} \
-specular {r g b a} -spot_dir {x,y,z} -spot_exp {v} -spot_cutoff {v} \
-const_atten {v} -lin_atten {v} -quad_atten {v}\n\
or\n\
xc_setGLparam lightmodel -two_side 0|1 -two_side_iso 0|1 -ambient {r g b a} -local_viewer 0|1\
-enable_light {v}|-disable_light {v}\n\
or\n\
xc_setGLparam frontface -what (isosurf_one | isosurf_pos | isosurf_neg | colorplane) -frontface CW|CCW\n\
or\n\
xc_setGLparam isonormal -what (isosurf_one | isosurf_pos | isosurf_neg)\n\
or\n\
xc_setGLparam blendfunc -what isosurf -sfunc <source_function> -dfunc <destination_function>\n\
or\n\
xc_setGLparam fog -get fogcolor", TCL_STATIC);

    return TCL_ERROR;
  }

  /***************************************************/
  /* MATERIAL                                        */
  if ( strcmp(argv[1],"material") == 0 ) {
    glparam.what = GLPAR_MATERIAL;
    for (i=2; i<argc; i+=2) {
      if ( strcmp(argv[i],"-what") == 0 ) {
	if ( strcmp(argv[i+1],"structure") == 0 ) {
	  glparam.mat_what = GLPAR_MAT_WHAT_STRUCTURE;
	}
	else if ( strcmp(argv[i+1],"isosurf") == 0 ) {
	  i+=2;
	  glparam.mat_what = GLPAR_MAT_WHAT_ISOSURF;
	  if ( strcmp(argv[i],"pos") == 0 ) {
	    glparam.mat_what_isosurf = GLPAR_MAT_WHAT_ISOPOS;
	  }
	  else if ( strcmp(argv[i],"neg") == 0 ) {
	    glparam.mat_what_isosurf = GLPAR_MAT_WHAT_ISONEG;
	  }
	  else if ( strcmp(argv[i],"one") == 0 ) {
	    glparam.mat_what_isosurf = GLPAR_MAT_WHAT_ISOONE;
	  }
	  else {
	    char rss[1024];
	    snprintf(rss, sizeof(rss),"unknown value \"%s\", must be one of pos, neg, one, while executing %s %s %s %s ....", 
		     argv[i], argv[0], argv[1], argv[2], argv[3]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	  if ( strcmp(argv[i+1],"front") == 0 ) {
	    glparam.mat_what_isoside = GLPAR_MAT_WHAT_ISOFRONTSIDE;
	  }
	  else if ( strcmp(argv[i+1],"back") == 0 ) {
	    glparam.mat_what_isoside = GLPAR_MAT_WHAT_ISOBACKSIDE;
	  }
	  else {
	    char rss[1024];
	    snprintf(rss, sizeof(rss),"unknown value \"%s\", must be one of front, back, while executing %s %s %s %s ....", 
		     argv[i], argv[0], argv[1], argv[2], argv[3]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	}
	else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"unknown value \"%s\", must be one of structure, isosurf, while executing %s %s %s %s ...", 
		   argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }
      else if ( strcmp(argv[i],"-specular") == 0 ) {
	if ( !xcSplitList( XC_GET_RGBA, interp, argv + i + 1, &specular ) )
	  return TCL_ERROR;	    
      }
      else if ( strcmp(argv[i],"-shininess") == 0 ) {
	if ( !xcSplitList( XC_GET_FLOAT, interp, argv + i + 1, &shininess ) )
	  return TCL_ERROR;
      }
      else if ( strcmp(argv[i],"-ambient") == 0 ) {
	if ( !xcSplitList( XC_GET_RGBA, interp, argv + i + 1, &ambient ) )
	  return TCL_ERROR;
      }
      else if ( strcmp(argv[i],"-diffuse") == 0 ) {
	if ( !xcSplitList( XC_GET_RGBA, interp, argv + i + 1, &diffuse ) )
	  return TCL_ERROR;
      }
      else if ( strcmp(argv[i],"-emission") == 0 ) {
	if ( !xcSplitList( XC_GET_RGBA, interp, argv + i + 1, &emission ) )
	  return TCL_ERROR;
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", must be one of -what, -specular, -shininess, -ambient, -diffuse, or -emission, while executing %s %s %s %s ...",argv[i], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
  }
  /***************************************************/
  /* LIGHT                                           */
  else if ( strcmp(argv[1],"light") == 0 ) {
    glparam.what = GLPAR_LIGHT;    
    for (i=2; i<argc; i+=2) {      
      if ( strcmp(argv[i],"-light") == 0 ) {
	if ( !xcSplitList( XC_GET_INT, interp, argv + i + 1, &lightN ) )
	  return TCL_ERROR;
      } 
      else if ( strcmp(argv[i],"-ambient") == 0 ) {
	if ( !xcSplitList( XC_GET_RGBA, interp, argv + i + 1, &ambient ) )
	  return TCL_ERROR;
      }
      else if ( strcmp(argv[i],"-diffuse") == 0 ) {
	if ( !xcSplitList( XC_GET_RGBA, interp, argv + i + 1, &diffuse ) )
	  return TCL_ERROR;
      }
      else if ( strcmp(argv[i],"-specular") == 0 ) {
	if ( !xcSplitList( XC_GET_RGBA, interp, argv + i + 1, &specular ) )
	  return TCL_ERROR;     
      }
      else if ( strcmp(argv[i],"-fract_position") == 0 ) {
	if ( !xcSplitList( XC_GET_XYZW, interp, argv + i + 1, 
			   &fract_position ) )
	  return TCL_ERROR;   
      }
      else if ( strcmp(argv[i],"-spot_dir") == 0 ) {
	if ( !xcSplitList( XC_GET_XYZ, interp, argv + i + 1, &spot_dir ) )
	  return TCL_ERROR;   
      }
      else if ( strcmp(argv[i],"-spot_exp") == 0 ) {
	if ( !xcSplitList( XC_GET_FLOAT, interp, argv + i + 1, &spot_exp ) )
	  return TCL_ERROR;
      }
      else if ( strcmp(argv[i],"-spot_cutoff") == 0 ) {
	if ( !xcSplitList( XC_GET_FLOAT, interp, argv + i + 1, &spot_cutoff ) )
	  return TCL_ERROR;
      }
      else if ( strcmp(argv[i],"-const_atten") == 0 ) {
	if ( !xcSplitList( XC_GET_FLOAT, interp, argv + i + 1, &const_atten ) )
	  return TCL_ERROR;
      }
      else if ( strcmp(argv[i],"-lin_atten") == 0 ) {
	if ( !xcSplitList( XC_GET_FLOAT, interp, argv + i + 1, &lin_atten ) )
	  return TCL_ERROR;
      }
      else if ( strcmp(argv[i],"-quad_atten") == 0 ) {
	if ( !xcSplitList( XC_GET_FLOAT, interp, argv + i + 1, &quad_atten ) )
	  return TCL_ERROR;
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", must be one of -light -ambient, -diffuse, -specular, -fract_position, -spot_dir, -spot_exp, -spot_cutoff, -cont_atten, -lin_atten, -quad_atten, while executing %s %s %s %s ...",argv[i], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
  }
  /***************************************************/
  /* LIGHTMODEL                                      */
  else if ( strcmp(argv[1],"lightmodel") == 0 ) {
    glparam.what = GLPAR_LIGHTMODEL;
    for (i=2; i<argc; i+=2) {
      if ( strcmp(argv[i],"-two_side") == 0 ) {
	if ( !xcSplitList( XC_GET_FLOAT, interp, argv + i + 1, &two_side ) )
	  return TCL_ERROR;      
      }
      else if ( strcmp(argv[i],"-two_side_iso") == 0 ) {
	if ( !xcSplitList( XC_GET_FLOAT, interp, argv + i + 1, 
			   &two_side_iso ) )
	  return TCL_ERROR;      
      }
      else if ( strcmp(argv[i],"-ambient") == 0 ) {
	if ( !xcSplitList( XC_GET_RGBA, interp, argv + i + 1, &ambient ) )
	  return TCL_ERROR;
      }
      else if ( strcmp(argv[i],"-local_viewer") == 0 ) {
	if ( !xcSplitList( XC_GET_FLOAT, interp, argv + i + 1, 
			   &local_viewer ) )
	  return TCL_ERROR;
      } 
      else if ( strcmp(argv[i],"-enable_light") == 0 ) {
	if ( !xcSplitList( XC_GET_FLOAT, interp, argv + i + 1, 
			   &enable_light ) )
	  return TCL_ERROR;
      }
      else if ( strcmp(argv[i],"-disable_light") == 0 ) {
	if ( !xcSplitList( XC_GET_FLOAT, interp, argv + i + 1, 
			   &disable_light ) )
	  return TCL_ERROR;
      }
      else {	
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", must be one of -two_side, -two_side_iso, -ambient, -local_viewer, -enable_light -disable_light, while executing %s %s %s %s ...",argv[i], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;      
      }
    }
  }
  /***************************************************/
  /* FRONTFACE                                       */
  else if ( strcmp(argv[1],"frontface") == 0 ) {
    glparam.what = GLPAR_FRONTFACE;    
    for (i=2; i<argc; i+=2) {
      if ( strcmp(argv[i],"-what") == 0 ) {
	if ( strcmp(argv[i+1],"isosurf_one") == 0 || 
	     strcmp(argv[i+1],"isosurf_pos") == 0 ) 
	  glparam.face_what = GLPAR_FACE_WHAT_ISOSURF_POS;
	else if ( strcmp(argv[i+1],"isosurf_neg") == 0 )
	  glparam.face_what = GLPAR_FACE_WHAT_ISOSURF_NEG;
	else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"unknown value \"%s\", must be one of: isosurf_one, isosurf_pos or isosurf_neg; while executing %s %s %s %s", 
		   argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }
      else if ( strcmp(argv[i],"-frontface") == 0 ) {	
	if ( strcmp(argv[i+1],"CCW") == 0 )
	  glparam.face_frontface = GL_CCW;
	else if ( strcmp(argv[i+1],"CW") == 0 )
	  glparam.face_frontface = GL_CW;	
	else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"unknown value \"%s\", must be one of CCW, CW, while executing %s %s %s %s", 
		   argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", must be one of -what, -frontface, while executing %s %s %s %s", 
		 argv[i], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
  }
  /***************************************************/
  /* ISONORMAL                                       */
  else if ( strcmp(argv[1],"isonormal") == 0 ) {
    if ( strcmp(argv[2], "-what") == 0 ) {
      if ( strcmp(argv[3], "isosurf_one") == 0 )
	glparam.face_what = GLPAR_FACE_WHAT_ISOSURF_POS;
      else if ( strcmp(argv[3], "isosurf_pos") == 0 )
	glparam.face_what = GLPAR_FACE_WHAT_ISOSURF_POS;
      else if ( strcmp(argv[3], "isosurf_neg") == 0 )
	glparam.face_what = GLPAR_FACE_WHAT_ISOSURF_NEG;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss), "unknown value %s for -what option, must be one of isosurf_one, isosurf_pos or isosurf_neg, while executing %s %s %s",
		 argv[3], argv[0], argv[1], argv[2]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      /* now we must multiply the normals by -1.0 */
      for (i=0; i<surf.numV[glparam.face_what]; i++)
	RevertVectorfv( (float *) &(normal[glparam.face_what][i]) );
    } else {
      char rss[1024];
      snprintf(rss, sizeof(rss), "unknown option %s, must be -what, while executing %s %s %s", argv[2], argv[0], argv[1], argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }	
  /***************************************************/
  /* BLENDFUNC                                       */
  else if ( strcmp(argv[1],"blendfunc") == 0 ) {
    glparam.what = GLPAR_BLENDFUNC;
    for (i=2; i<argc; i+=2) {
      if ( strcmp(argv[i],"-what") == 0 ) {
	if ( strcmp(argv[i+1],"isosurf") == 0 ) 
	  glparam.blend_what = GLPAR_BLEND_WHAT_ISOSURF;
	else if ( strcmp(argv[i+1],"colorplane") == 0 )
	  glparam.blend_what = GLPAR_BLEND_WHAT_COLORPLANE;
	else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"unknown value \"%s\", must be one of isosurf, colorplane; while executing %s %s %s %s", 
		   argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }
      else if ( strcmp(argv[i],"-sfunc") == 0 ) {
	is_blend_sfunc = 1;
	if ( strcmp(argv[i+1],"GL_ZERO") == 0 ) 
	  glparam.blend_sfunc = GL_ZERO;
	else if ( strcmp(argv[i+1],"GL_ONE") == 0 )
	  glparam.blend_sfunc = GL_ONE;
	else if ( strcmp(argv[i+1],"GL_DST_COLOR") == 0 )
	  glparam.blend_sfunc = GL_DST_COLOR;
	else if ( strcmp(argv[i+1],"GL_ONE_MINUS_DST_COLOR") == 0 )
	  glparam.blend_sfunc = GL_ONE_MINUS_DST_COLOR;
	else if ( strcmp(argv[i+1],"GL_SRC_ALPHA") == 0 )
	  glparam.blend_sfunc = GL_SRC_ALPHA;
	else if ( strcmp(argv[i+1],"GL_ONE_MINUS_SRC_ALPHA") == 0 )
	  glparam.blend_sfunc = GL_ONE_MINUS_SRC_ALPHA;
	else if ( strcmp(argv[i+1],"GL_DST_ALPHA") == 0 )
	  glparam.blend_sfunc = GL_DST_ALPHA;
	else if ( strcmp(argv[i+1],"GL_ONE_MINUS_DST_ALPHA") == 0 )
	  glparam.blend_sfunc = GL_ONE_MINUS_DST_ALPHA;
	else if ( strcmp(argv[i+1],"GL_SRC_ALPHA_SATURATE") == 0 )
	  glparam.blend_sfunc = GL_SRC_ALPHA_SATURATE;
	else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"unknown SOURCE BLEND FUNCTION \"%s\", must be one of GL_ZERO, GL_ONE, GL_DST_COLOR, GL_ONE_MINUS_DST_COLOR, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, GL_SRC_ALPHA_SATURATE\", while executing %s %s %s %s", 
		   argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }
      else if ( strcmp(argv[i],"-dfunc") == 0 ) {
	is_blend_dfunc = 1;
	if ( strcmp(argv[i+1],"GL_ZERO") == 0 ) 
	  glparam.blend_dfunc = GL_ZERO;
	else if ( strcmp(argv[i+1],"GL_ONE") == 0 )
	  glparam.blend_dfunc = GL_ONE;
	else if ( strcmp(argv[i+1],"GL_SRC_COLOR") == 0 )
	  glparam.blend_dfunc = GL_SRC_COLOR;
	else if ( strcmp(argv[i+1],"GL_ONE_MINUS_SRC_COLOR") == 0 )
	  glparam.blend_dfunc = GL_ONE_MINUS_SRC_COLOR;
	else if ( strcmp(argv[i+1],"GL_SRC_ALPHA") == 0 )
	  glparam.blend_dfunc = GL_SRC_ALPHA;
	else if ( strcmp(argv[i+1],"GL_ONE_MINUS_SRC_ALPHA") == 0 )
	  glparam.blend_dfunc = GL_ONE_MINUS_SRC_ALPHA;
	else if ( strcmp(argv[i+1],"GL_DST_ALPHA") == 0 )
	  glparam.blend_dfunc = GL_DST_ALPHA;
	else if ( strcmp(argv[i+1],"GL_ONE_MINUS_DST_ALPHA") == 0 )
	  glparam.blend_dfunc = GL_ONE_MINUS_DST_ALPHA;
	else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"unknown DESTINATION BLEND FUNCTION \"%s\", must be one of GL_ZERO, GL_ONE, GL_SRC_COLOR, GL_ONE_MINUS_SCR_COLOR, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA\", while executing %s %s %s %s", 
		   argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }
    }
  }
  /***************************************************/
  /* FOG                                       */
  else if ( strcmp(argv[1],"fog") == 0 ) {
    glparam.what = GLPAR_FOG;
    if ( strcmp(argv[2],"-color") == 0 ) {
      if ( !xcSplitList( XC_GET_RGBA, interp, argv + 3, &fogcolor ) )
	return TCL_ERROR;
    } else {
      char rss[1024];
      snprintf(rss, sizeof(rss),"unknown option \"%s\", must be -color", argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  else {
    char rss[1024];
    snprintf(rss, sizeof(rss),"unknown argument \"%s\", must be one of material, light, lightmodel, frontface or isonormal or blendfunc or fog, while executing %s %s ...",
	     argv[1], argv[0], argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* now updata parameters */
  if ( glparam.what == GLPAR_MATERIAL ) {
    if ( glparam.mat_what == GLPAR_MAT_WHAT_STRUCTURE ) {
      if ( specular.is )
	CopyGLParVec( XC_GET_RGBA, (float *) &(mat_struct.specular), 
		      specular );
      if ( shininess.is )
	CopyGLParVec( XC_GET_FLOAT, (float *) &(mat_struct.shininess), 
		      shininess );
      if ( ambient.is ) 
	CopyGLParVec( XC_GET_RGBA, (float *) &(mat_struct.ambient), 
		      ambient );
      if ( diffuse.is ) 
	CopyGLParVec( XC_GET_RGBA, (float *) &(mat_struct.diffuse), diffuse );
      if ( emission.is ) 
	CopyGLParVec( XC_GET_RGBA, (float *) &(mat_struct.emission), emission );
    }
    else if ( glparam.mat_what == GLPAR_MAT_WHAT_ISOSURF ) {
      int tr;
      tr = isoDisp.transparent;
      if ( glparam.mat_what_isosurf == GLPAR_MAT_WHAT_ISOONE ) {
	if ( glparam.mat_what_isoside == GLPAR_MAT_WHAT_ISOFRONTSIDE ) {
	  if ( specular.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(front_mat_isosurf[tr].specular), 
			  specular );
	  if ( shininess.is ) 
	    CopyGLParVec( XC_GET_FLOAT, 
			  (float *) &(front_mat_isosurf[tr].shininess), 
			  shininess);
	  if ( ambient.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(front_mat_isosurf[tr].ambient), 
			  ambient );
	  if ( diffuse.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(front_mat_isosurf[tr].diffuse), 
			  diffuse );
	  if ( emission.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(front_mat_isosurf[tr].emission), 
			  emission );
	}
	else if ( glparam.mat_what_isoside == GLPAR_MAT_WHAT_ISOBACKSIDE ) {
	  if ( specular.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(back_mat_isosurf[tr].specular), 
			  specular );
	  if ( shininess.is ) 
	    CopyGLParVec( XC_GET_FLOAT, 
			  (float *) &(back_mat_isosurf[tr].shininess), 
			  shininess);
	  if ( ambient.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(back_mat_isosurf[tr].ambient), 
			  ambient );
	  if ( diffuse.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(back_mat_isosurf[tr].diffuse), 
			  diffuse );
	  if ( emission.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(back_mat_isosurf[tr].emission), 
			  emission );
	}
      }
      else if ( glparam.mat_what_isosurf == GLPAR_MAT_WHAT_ISOPOS ) {
	if ( glparam.mat_what_isoside == GLPAR_MAT_WHAT_ISOFRONTSIDE ) {
	  if ( specular.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(front_mat_pos_isosurf[tr].specular), 
			  specular );
	  if ( shininess.is ) 
	    CopyGLParVec( XC_GET_FLOAT, 
			  (float *) &(front_mat_pos_isosurf[tr].shininess), 
			  shininess);
	  if ( ambient.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(front_mat_pos_isosurf[tr].ambient), 
			  ambient );
	  if ( diffuse.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(front_mat_pos_isosurf[tr].diffuse), 
			  diffuse );
	  if ( emission.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(front_mat_pos_isosurf[tr].emission), 
			  emission );
	}
	else if ( glparam.mat_what_isoside == GLPAR_MAT_WHAT_ISOBACKSIDE ) {
	  if ( specular.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(back_mat_pos_isosurf[tr].specular), 
			  specular );
	  if ( shininess.is ) 
	    CopyGLParVec( XC_GET_FLOAT, 
			  (float *) &(back_mat_pos_isosurf[tr].shininess), 
			  shininess);
	  if ( ambient.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(back_mat_pos_isosurf[tr].ambient), 
			  ambient );
	  if ( diffuse.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(back_mat_pos_isosurf[tr].diffuse), 
			  diffuse );
	  if ( emission.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(back_mat_pos_isosurf[tr].emission), 
			  emission );
	}      
      }      
      else if ( glparam.mat_what_isosurf == GLPAR_MAT_WHAT_ISONEG ) {
	if ( glparam.mat_what_isoside == GLPAR_MAT_WHAT_ISOFRONTSIDE ) {
	  if ( specular.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(front_mat_neg_isosurf[tr].specular), 
			  specular );
	  if ( shininess.is ) 
	    CopyGLParVec( XC_GET_FLOAT, 
			  (float *) &(front_mat_neg_isosurf[tr].shininess), 
			  shininess);
	  if ( ambient.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(front_mat_neg_isosurf[tr].ambient), 
			  ambient );
	  if ( diffuse.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(front_mat_neg_isosurf[tr].diffuse), 
			  diffuse );
	  if ( emission.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(front_mat_neg_isosurf[tr].emission), 
			  emission );
	}
	else if ( glparam.mat_what_isoside == GLPAR_MAT_WHAT_ISOBACKSIDE ) {
	  if ( specular.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(back_mat_neg_isosurf[tr].specular), 
			  specular );
	  if ( shininess.is ) 
	    CopyGLParVec( XC_GET_FLOAT, 
			  (float *) &(back_mat_neg_isosurf[tr].shininess), 
			  shininess);
	  if ( ambient.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(back_mat_neg_isosurf[tr].ambient), 
			  ambient );
	  if ( diffuse.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(back_mat_neg_isosurf[tr].diffuse), 
			  diffuse );
	  if ( emission.is ) 
	    CopyGLParVec( XC_GET_RGBA, 
			  (float *) &(back_mat_neg_isosurf[tr].emission), 
			  emission );
	}
      }   
    }
  }
  else if ( glparam.what == GLPAR_LIGHT ) {
    /* if -light option wasn't present all data regard to GL_LIGHT0 */
    int il = 0;
    if ( lightN.is ) {
      il = (int) lightN.vec[0];
      if ( il < 0 || il > GLPAR_MAXLIGHT ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "light number %d is out of range, should be between 0 and %d", il, GLPAR_MAXLIGHT);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    if ( ambient.is ) 
      CopyGLParVec( XC_GET_RGBA, (float *) &(light[il].ambient), ambient );
    if ( diffuse.is )
      CopyGLParVec( XC_GET_RGBA, (float *) &(light[il].diffuse), diffuse );
    if ( specular.is )
      CopyGLParVec( XC_GET_RGBA, (float *) &(light[il].specular), specular );
    if ( fract_position.is ) 
      CopyGLParVec( XC_GET_XYZW, (float *) &(light[il].fract_position), 
		    fract_position );
    if ( spot_dir.is )
      CopyGLParVec( XC_GET_XYZ, (float *) &(light[il].spot_dir), spot_dir );
    if ( spot_exp.is )
      CopyGLParVec( XC_GET_FLOAT, (float *) &(light[il].spot_exp), spot_exp );
    if ( spot_cutoff.is )
      CopyGLParVec( XC_GET_FLOAT, (float *) &(light[il].spot_cutoff), 
		    spot_cutoff );
    if ( const_atten.is )
      CopyGLParVec( XC_GET_FLOAT, (float *) &(light[il].const_atten), 
		    const_atten );
    if ( lin_atten.is ) 
      CopyGLParVec( XC_GET_FLOAT, (float *) &(light[il].lin_atten), 
		    lin_atten );
    if ( quad_atten.is )
      CopyGLParVec( XC_GET_FLOAT, (float *) &(light[il].quad_atten), 
		    quad_atten );
  }      

  else if ( glparam.what == GLPAR_LIGHTMODEL ) {
    int   il;
    if ( two_side.is )
      CopyGLParVec( XC_GET_FLOAT, (float *) &(lightmodel.two_side), 
		    two_side );    
    if ( two_side_iso.is )
      CopyGLParVec( XC_GET_FLOAT, (float *) &(lightmodel.two_side_iso), 
		    two_side_iso );    
    if ( ambient.is ) 
      CopyGLParVec( XC_GET_RGBA, (float *) &(lightmodel.ambient), 
		    ambient );
    if ( local_viewer.is )
      CopyGLParVec( XC_GET_FLOAT, (float *) &(lightmodel.local_viewer), 
		    local_viewer );

    if ( enable_light.is || disable_light.is ) {
      /*CopyGLParVec( XC_GET_INT, (float *) &fl, enable_light );*/
      if ( enable_light.is  ) il = (int)enable_light.vec[0];
      else il = (int)disable_light.vec[0];
      if ( il < 0 || il > GLPAR_MAXLIGHT ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"light number %d is out of range, should be between 0 and %d", il, GLPAR_MAXLIGHT);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      if ( enable_light.is  ) light[il].isenabled = 1;
      if ( disable_light.is ) light[il].isenabled = 0;
    }
  }

  else if ( glparam.what == GLPAR_FRONTFACE ) {
    if ( glparam.face_what == GLPAR_FACE_WHAT_ISOSURF_POS || 
	 glparam.face_what == GLPAR_FACE_WHAT_ISOSURF_NEG )
      frontface_isosurf[glparam.face_what] = glparam.face_frontface;
  }

  else if ( glparam.what == GLPAR_BLENDFUNC ) {
    if ( glparam.blend_what == GLPAR_BLEND_WHAT_ISOSURF ) {
      if ( is_blend_sfunc )
	blend_isosurf.sfunc = glparam.blend_sfunc;
      if ( is_blend_dfunc )
	blend_isosurf.dfunc = glparam.blend_dfunc;
    } 
    else if ( glparam.what == GLPAR_BLEND_WHAT_COLORPLANE ) {
      if ( is_blend_sfunc )
	blend_colorplane.sfunc = glparam.blend_sfunc;
      if ( is_blend_dfunc )
	blend_colorplane.dfunc = glparam.blend_dfunc;
    }
  }

  else if (glparam.what == GLPAR_FOG ) {
    CopyGLParVec( XC_GET_RGBA, (float *) fog.color, fogcolor );
  }

  /* if we want  the changes to take place we must run 
   * the following functions
   */

  /* THIS IS NOT ENOUGH, BECAUSE in LOADLIGHTS, 
   * GL_POSITION & GL_SPOT DIRECTION are not set, 
   * but are rather set in xcMakeProjection3D
   * (see below !!!)
   */
  if (togl_exists) {
    LoadLights(); /* call this explicitly (needed by Fermi-Surfaces) */
    LoadStructMaterial();
  }

  /* t.k.: sometimes when some openGL parameters are changed, 
     the lights go crazy; maybe xcMakeProjection3D should be called ???? */
  if ( dimType == XC_3D && VPf.stropened ) {
    if (is.stickmode && !is.ballmode) xcMakeProjection3D("sticks");
    if (is.ballmode) xcMakeProjection3D("balls");
    if (is.spacefillmode) xcMakeProjection3D("space");
  }
  /* now update display */
  /* (*xcDisplay)();  this is now executed in Tcl script by .mesa render */

  return TCL_OK;
}



/*****************************************************************************/
/* XC_GetGLparamCmd --> inplementation of 'xc_getGLparam' custom Tcl command 
 * --------------------- 
 * Usage:  xc_getGLparam material \
 *			     -what structure|(isosurf pos|neg|one front|back) \
 *			     -get shininess|specular|ambient|diffuse|emission| \
 *                                def_shininess|def_specular|def_ambient| \
 *                                def_diffuse|def_emission
 *
 *	    or
 *	    
 *	    xc_getGLparam light \
 *			     -light {v} \
 *			     -get ambient|diffuse|fract_position|specular| \
 *                                spot_dir|spot_exp|spot_cutoff|const_atten \
 *                                lin_atten|quad_atten|
 *                                def_*  --> '*' means the same as without def_
 *                                            possix
 *
 *	    or
 *	    
 *	    xc_getGLparam lightmodel \
 *                           -get two_side|two_side_iso|ambient|local_viewer| \
 *                                {is_light_enabled <n>}| \
 *                                def_*  --> '*' means the same as without def_
 *                                            possix
 *
 *          or
 *
 *          xc_getGLparam frontface \
 *                           -what isosurf
 *
 *          or
 *
 *          xc_getGLparam blendfunc \
 *                           -what isosurf
 *
 *          or
 *
 *          xc_getGLparam fog -get fogcolor|def_fogcolor
 *
 * WARNING: xc_getGLparam material, xc_getGLparam light, 
 *          xc_getGLparam lightmodel return the requested parameters
 *          WHEREAS
 *          xc_getGLparam frontface & xcGLparam blendfunc return:
 *          
 *          xc_getGLparam frontface  --> {current_fronface, default_frontface}
 *          xc_getGLparam blendfunc  --> {current_sfunc, current_dfunc,
 *                                        default_sfunc, default_dfunc}
 */ 
int 
XC_GetGLparamCmd(ClientData clientData, Tcl_Interp *interp,
		 int argc, const char *argv[])
{
  int i;
  char *result = Tcl_Alloc( sizeof(char) * 120 );

  struct GLParam {
    int what;
    int mat_what;
    int mat_what_isosurf;
    int mat_what_isoside;
  } glparam;

  boolean specular      = 0;
  boolean shininess     = 0;
  boolean ambient       = 0;
  boolean diffuse       = 0;
  boolean emission      = 0;
  GetGlParam lightN     = { 0 };
  boolean fract_position= 0;
  boolean spot_dir      = 0;
  boolean spot_exp      = 0;
  boolean spot_cutoff   = 0;
  boolean const_atten   = 0;
  boolean lin_atten     = 0;
  boolean quad_atten    = 0;
  boolean fogcolor      = 0;

  boolean def_specular      = 0;
  boolean def_shininess     = 0;
  boolean def_ambient       = 0;
  boolean def_diffuse       = 0;
  boolean def_emission      = 0;
  boolean def_fract_position= 0;
  boolean def_spot_dir      = 0;
  boolean def_spot_exp      = 0;
  boolean def_spot_cutoff   = 0;
  boolean def_const_atten   = 0;
  boolean def_lin_atten     = 0;
  boolean def_quad_atten    = 0;
  boolean def_fogcolor      = 0;

  if ( (argc % 2) != 0 ) {
    Tcl_SetResult(interp,"Usage: xc_getGLparam material \n\
-what structure|(isosurf pos|neg|one front|back) \n\
-get shininess|specular|ambient|diffuse|emission\n\n\
or\n\n\
xc_getGLparam light -light {v} \n\
-get ambient|diffuse|fract_position|specular|spot_dir|spot_exp|spot_cutoff|const_atten|lin_atten|quad_atten\n\n\
or\n\n\
xc_getGLparam lightmodel -get two_side|two_side_iso|ambient|local_viewer|{is_light_enabled <n>}", TCL_STATIC);
    return TCL_ERROR;  
  }

  if ( strcmp(argv[1],"material") == 0 ) {
    glparam.what = GLPAR_MATERIAL;
    for (i=2; i<argc; i+=2) {
      if ( strcmp(argv[i],"-what") == 0 ) {
	if ( strcmp(argv[i+1],"structure") == 0 ) {
	  glparam.mat_what = GLPAR_MAT_WHAT_STRUCTURE;
	}
	else if ( strcmp(argv[i+1],"isosurf") == 0 ) {
	  i+=2;
	  glparam.mat_what = GLPAR_MAT_WHAT_ISOSURF;
	  if ( strcmp(argv[i],"pos") == 0 ) {
	    glparam.mat_what_isosurf = GLPAR_MAT_WHAT_ISOPOS;
	  }
	  else if ( strcmp(argv[i],"neg") == 0 ) {
	    glparam.mat_what_isosurf = GLPAR_MAT_WHAT_ISONEG;
	  }
	  else if ( strcmp(argv[i],"one") == 0 ) {
	    glparam.mat_what_isosurf = GLPAR_MAT_WHAT_ISOONE;
	  }
	  else {
	    char rss[1024];
	    snprintf(rss, sizeof(rss),"unknown value \"%s\", must be one of pos, neg, one, while executing %s %s %s %s ....", 
		     argv[i], argv[0], argv[1], argv[2], argv[3]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	  if ( strcmp(argv[i+1],"front") == 0 ) {
	    glparam.mat_what_isoside = GLPAR_MAT_WHAT_ISOFRONTSIDE;
	  }
	  else if ( strcmp(argv[i+1],"back") == 0 ) {
	    glparam.mat_what_isoside = GLPAR_MAT_WHAT_ISOBACKSIDE;
	  }
	  else {
	    char rss[1024];
	    snprintf(rss, sizeof(rss),"unknown value \"%s\", must be one of front, back, while executing %s %s %s %s ....", 
		     argv[i], argv[0], argv[1], argv[2], argv[3]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  } 
	}
	else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"unknown value \"%s\", must be one of structure, isosurf, while executing %s %s %s %s ...", 
		   argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      } 
      else if ( strcmp(argv[i],"-get") == 0 ) {
	if ( strcmp(argv[i+1],"specular") == 0 )
	  specular = 1;
	else if ( strcmp(argv[i+1],"shininess") == 0 ) 
	  shininess = 1;
	else if ( strcmp(argv[i+1],"ambient") == 0 ) 
	  ambient = 1;
	else if ( strcmp(argv[i+1],"diffuse") == 0 )
	  diffuse = 1;
	else if ( strcmp(argv[i+1],"emission") == 0 )
	  emission = 1;
	/* QUERY DEFAULTS */
	else if ( strcmp(argv[i+1],"def_specular") == 0 )
	  def_specular = 1;
	else if ( strcmp(argv[i+1],"def_shininess") == 0 ) 
	  def_shininess = 1;
	else if ( strcmp(argv[i+1],"def_ambient") == 0 ) 
	  def_ambient = 1;
	else if ( strcmp(argv[i+1],"def_diffuse") == 0 )
	  def_diffuse = 1;
	else if ( strcmp(argv[i+1],"def_emission") == 0 )
	  def_emission = 1;
	else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"unknown value \"%s\" for -get option - must be one of \"specular, shininess, ambient, diffuse, emission\"def_*\"\", while executing %s %s %s %s ...",argv[i+1],argv[0],argv[1],argv[2],argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      } else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", must be one of -what, -get, while executing %s %s %s %s ...",argv[i], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
  }
  else if ( strcmp(argv[1],"light") == 0 ) {
    glparam.what = GLPAR_LIGHT;    
    for (i=2; i<argc; i+=2) {      
      if ( strcmp(argv[i],"-light") == 0 ) {
	if ( !xcSplitList( XC_GET_INT, interp, argv + i + 1, &lightN ) )
	  return TCL_ERROR;
      }
      else if ( strcmp(argv[i],"-get") == 0 ) {
	if ( strcmp(argv[i+1],"ambient") == 0 ) 
	  ambient = 1;
	else if ( strcmp(argv[i+1],"diffuse") == 0 ) 
	  diffuse = 1;
	else if ( strcmp(argv[i+1],"specular") == 0 ) 
	  specular = 1;
	else if ( strcmp(argv[i+1],"fract_position") == 0 )
	  fract_position = 1;
	else if ( strcmp(argv[i+1],"spot_dir") == 0 ) 
	  spot_dir = 1;
	else if ( strcmp(argv[i+1],"spot_exp") == 0 )
	  spot_exp = 1;
	else if ( strcmp(argv[i+1],"spot_cutoff") == 0 )
	  spot_cutoff = 1;
	else if ( strcmp(argv[i+1],"const_atten") == 0 )
	  const_atten = 1;
	else if ( strcmp(argv[i+1],"lin_atten") == 0 )
	  lin_atten = 1;
	else if ( strcmp(argv[i+1],"quad_atten") == 0 )
	  quad_atten = 1;
	/* QUERY DEFAULTS */
	else if ( strcmp(argv[i+1],"def_ambient") == 0 ) 
	  def_ambient = 1;
	else if ( strcmp(argv[i+1],"def_diffuse") == 0 ) 
	  def_diffuse = 1;
	else if ( strcmp(argv[i+1],"def_specular") == 0 ) 
	  def_specular = 1;
	else if ( strcmp(argv[i+1],"def_fract_position") == 0 )
	  def_fract_position = 1;
	else if ( strcmp(argv[i+1],"def_spot_dir") == 0 ) 
	  def_spot_dir = 1;
	else if ( strcmp(argv[i+1],"def_spot_exp") == 0 )
	  def_spot_exp = 1;
	else if ( strcmp(argv[i+1],"def_spot_cutoff") == 0 )
	  def_spot_cutoff = 1;
	else if ( strcmp(argv[i+1],"def_const_atten") == 0 )
	  def_const_atten = 1;
	else if ( strcmp(argv[i+1],"def_lin_atten") == 0 )
	  def_lin_atten = 1;
	else if ( strcmp(argv[i+1],"def_quad_atten") == 0 )
	  def_quad_atten = 1;
	else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"unknown value \"%s\" for -get option - must be one of \"specular, ambient, diffuse, fract_position, spot_dir, spot_exp, spot_cutoff, const_atten, lin_ateen, quad_atten, \"def_*\"\", while executing %s %s %s %s ...",argv[i+1],argv[0],argv[1],argv[2],argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", must be one of -light or -get, while executing %s %s %s %s ...",argv[i], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
  }
  else if ( strcmp(argv[1],"lightmodel") == 0 ) {
    glparam.what = GLPAR_LIGHTMODEL;
    if ( strcmp(argv[2],"-get") == 0 ) {
      if ( strcmp(argv[3],"two_side") == 0 ) 	
	sprintf(result,"%d", iroundf(lightmodel.two_side[0]));
      else if ( strcmp(argv[3],"two_side_iso") == 0 ) 	
	sprintf(result,"%d", iroundf(lightmodel.two_side_iso[0]));
      else if ( strcmp(argv[3],"ambient") == 0 ) 
	sprintf(result,"%f %f %f %f", 
		lightmodel.ambient[0], lightmodel.ambient[1],
		lightmodel.ambient[2], lightmodel.ambient[3]);
      else if ( strcmp(argv[3],"local_viewer") == 0 ) 
	sprintf(result,"%d", iroundf(lightmodel.local_viewer[0]));	
      /* now here comes list ????? */
      else if ( strncmp(argv[3],"is_light_enabled",4) == 0 ) {
	int   argcList;
	const char  **argvList;
	int   light_number;
	char  *cln = Tcl_Alloc( sizeof(char) * 10 );
	Tcl_SplitList(interp, argv[3], &argcList, &argvList);
	if ( argcList != 2 ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"expected \"is_light_enabled <n>\", but got \"%s\", while executing %s %s %s %s ...",
		   argv[3], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	if ( Tcl_GetInt(interp, argvList[1], &light_number) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"expected integer, but got \"%s\", while executing %s %s %s %s ...",
		   argv[3], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	} else if ( light_number < 0 || light_number > GLPAR_MAXLIGHT ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"only up to %d lights can be defined, but you try to query light number %d, while executing %s %s %s %s ...",
		   MAX_NUMBER_OF_LIGHTS, light_number, 
		   argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	} else {
	  sprintf(cln,"%d",light[light_number].isenabled);
	  Tcl_SetResult(interp, cln, TCL_DYNAMIC);
	  return TCL_OK;
	}
      }
      /**********************/
      /* GET DEFAULT VALUES */
      /**********************/
      else if ( strcmp(argv[3],"def_two_side") == 0 )
	sprintf(result,"%d", iroundf(def_lightmodel_two_side[0]));
      else if ( strcmp(argv[3],"def_two_side_iso") == 0 )
	sprintf(result,"%d", iroundf(def_lightmodel_two_side_iso[0]));
      else if ( strcmp(argv[3],"def_ambient") == 0 )
	sprintf(result,"%f %f %f %f",
		def_lightmodel_ambient[0], def_lightmodel_ambient[1],
		def_lightmodel_ambient[2], def_lightmodel_ambient[3]);
      else if ( strcmp(argv[3],"def_local_viewer") == 0 )
	sprintf(result,"%d", iroundf(def_lightmodel_local_viewer[0]));
      /* now here comes list ????? */
      else if ( strncmp(argv[3],"def_is_light_enabled",4) == 0 ) {
	int   argcList;
	const char  **argvList;
	int   light_number;
	char  *cln = Tcl_Alloc( sizeof(char) * 10 );
	Tcl_SplitList(interp, argv[3], &argcList, &argvList);
	if ( argcList != 2 ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"expected \"is_light_enabled <n>\", but got \"%s\", while executing %s %s %s %s ...",
		   argv[3], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	if ( Tcl_GetInt(interp, argvList[1], &light_number) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"expected integer, but got \"%s\", while executing %s %s %s %s ...",
		   argv[3], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	} else if ( light_number < 0 || light_number > GLPAR_MAXLIGHT ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"only up to %d lights can be defined, but you try to query light number %d, while executing %s %s %s %s ...",
		   MAX_NUMBER_OF_LIGHTS, light_number, 
		   argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	} else {
	  if ( light_number == 0 )
	    sprintf(cln,"%d", 1);
	  else 
	    sprintf(cln,"%d", 0);
	  Tcl_SetResult(interp, cln, TCL_DYNAMIC);
	  return TCL_OK;
	}
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown value \"%s\" for -get option, must be one of \"two_side. two_side_iso, ambient, local_viewer, {is_light_enabled <n>}\", while executing %s %s %s %s ....",argv[3], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else {	
      char rss[1024];
      snprintf(rss, sizeof(rss),"unknown option \"%s\", must be one of -get, while executing %s %s %s %s ...",argv[3], argv[0], argv[1], argv[2], argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;            
    }
  }
  else if ( strcmp(argv[1],"frontface") == 0 ) {
    char res1[10];
    char res2[10];

    /* default res1, res2 values */
    sprintf(res1,"NULL");
    sprintf(res2,"NULL");

    if ( strcmp(argv[2],"-what") == 0 ) {
      if ( strncmp(argv[3],"isosurf_",8) == 0 ) {
	if ( strcmp(argv[3],"isosurf_one") == 0 ||
	     strcmp(argv[3],"isosurf_pos") == 0 ) i = 0;
	else if ( strcmp(argv[3],"isosurf_neg") == 0 ) i = 1;
	else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"unknown value \"%s\", must be one of: isosurf_one, isosurf_pos, isosurf_neg; while executing %s %s %s %s", argv[3], argv[0], argv[2], argv[3], argv[4]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	if ( frontface_isosurf[i] == GL_CCW )  sprintf(res1,"CCW");
	if ( frontface_isosurf[i] == GL_CW )   sprintf(res1,"CW");
	if ( def_frontface_isosurf == GL_CCW ) sprintf(res2,"CCW");
	if ( def_frontface_isosurf == GL_CW )  sprintf(res2,"CW");
	sprintf(result,"%s %s", res1, res2);
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown value \"%s\", must be one of: isosurf_one, isosurf_pos, isosurf_neg; while executing %s %s %s %s", argv[3], argv[0], argv[2], argv[3], argv[4]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    } else {
      char rss[1024];
      snprintf(rss, sizeof(rss),"unknown option \"%s\", must be -what, while executing %s %s %s %s", argv[2], argv[0], argv[2], argv[3], argv[4]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }
  else if ( strcmp(argv[1],"blendfunc") == 0 ) {
    char *res1 = (char *) malloc( sizeof(char) * 30 );
    char *res2 = (char *) malloc( sizeof(char) * 30 );    
    char *res3 = (char *) malloc( sizeof(char) * 30 );
    char *res4 = (char *) malloc( sizeof(char) * 30 );    

    /* default resX values */
    sprintf(res1,"NULL");
    sprintf(res2,"NULL");
    sprintf(res3,"NULL");
    sprintf(res4,"NULL");

    if ( strcmp(argv[2],"-what") == 0 ) {
      if ( strcmp(argv[3],"isosurf") == 0 ) {
	/* CURRENT SOURCE FUNCTION */
	if ( blend_isosurf.sfunc == GL_ZERO ) 
	  sprintf(res1,"GL_ZERO");
	else if ( blend_isosurf.sfunc == GL_ONE )
	  sprintf(res1,"GL_ONE");
	else if ( blend_isosurf.sfunc == GL_DST_COLOR )
	  sprintf(res1,"GL_DST_COLOR");
	else if ( blend_isosurf.sfunc == GL_ONE_MINUS_DST_COLOR )
	  sprintf(res1,"GL_ONE_MINUS_DST_COLOR");
	else if ( blend_isosurf.sfunc == GL_SRC_ALPHA )
	  sprintf(res1,"GL_SRC_ALPHA");
	else if ( blend_isosurf.sfunc == GL_ONE_MINUS_SRC_ALPHA )
	  sprintf(res1,"GL_ONE_MINUS_SRC_ALPHA");
	else if ( blend_isosurf.sfunc == GL_DST_ALPHA )
	  sprintf(res1,"GL_DST_ALPHA");
	else if ( blend_isosurf.sfunc == GL_ONE_MINUS_DST_ALPHA )
	  sprintf(res1,"GL_ONE_MINUS_DST_ALPHA");
	else if ( blend_isosurf.sfunc == GL_SRC_ALPHA_SATURATE )
	  sprintf(res1,"GL_SRC_ALPHA_SATURATE");
	/* CURRENT DESTINATION FUNCTION */
	if ( blend_isosurf.dfunc == GL_ZERO ) 
	  sprintf(res2,"GL_ZERO");
	else if ( blend_isosurf.dfunc == GL_ONE )
	  sprintf(res2,"GL_ONE");
	else if ( blend_isosurf.dfunc == GL_SRC_COLOR )
	  sprintf(res2,"GL_SRC_COLOR");
	else if ( blend_isosurf.dfunc == GL_ONE_MINUS_SRC_COLOR )
	  sprintf(res2,"GL_ONE_MINUS_SRC_COLOR");
	else if ( blend_isosurf.dfunc == GL_SRC_ALPHA )
	  sprintf(res2,"GL_SRC_ALPHA");
	else if ( blend_isosurf.dfunc == GL_ONE_MINUS_SRC_ALPHA )
	  sprintf(res2,"GL_ONE_MINUS_SRC_ALPHA");
	else if ( blend_isosurf.dfunc == GL_DST_ALPHA )
	  sprintf(res2,"GL_DST_ALPHA");
	else if ( blend_isosurf.dfunc == GL_ONE_MINUS_DST_ALPHA )
	  sprintf(res2,"GL_ONE_MINUS_DST_ALPHA");
	/* DEFAULT SOURCE FUNCTION */
	if ( def_blend_isosurf.sfunc == GL_ZERO ) 
	  sprintf(res3,"GL_ZERO");
	else if ( def_blend_isosurf.sfunc == GL_ONE )
	  sprintf(res3,"GL_ONE");
	else if ( def_blend_isosurf.sfunc == GL_DST_COLOR )
	  sprintf(res3,"GL_DST_COLOR");
	else if ( def_blend_isosurf.sfunc == GL_ONE_MINUS_DST_COLOR )
	  sprintf(res3,"GL_ONE_MINUS_DST_COLOR");
	else if ( def_blend_isosurf.sfunc == GL_SRC_ALPHA )
	  sprintf(res3,"GL_SRC_ALPHA");
	else if ( def_blend_isosurf.sfunc == GL_ONE_MINUS_SRC_ALPHA )
	  sprintf(res3,"GL_ONE_MINUS_SRC_ALPHA");
	else if ( def_blend_isosurf.sfunc == GL_DST_ALPHA )
	  sprintf(res3,"GL_DST_ALPHA");
	else if ( def_blend_isosurf.sfunc == GL_ONE_MINUS_DST_ALPHA )
	  sprintf(res3,"GL_ONE_MINUS_DST_ALPHA");
	else if ( def_blend_isosurf.sfunc == GL_SRC_ALPHA_SATURATE )
	  sprintf(res3,"GL_SRC_ALPHA_SATURATE");
	/* DEFAULT DESTINATION FUNCTION */
	if ( blend_isosurf.dfunc == GL_ZERO ) 
	  sprintf(res4,"GL_ZERO");
	else if ( def_blend_isosurf.dfunc == GL_ONE )
	  sprintf(res4,"GL_ONE");
	else if ( def_blend_isosurf.dfunc == GL_SRC_COLOR )
	  sprintf(res4,"GL_SRC_COLOR");
	else if ( def_blend_isosurf.dfunc == GL_ONE_MINUS_SRC_COLOR )
	  sprintf(res4,"GL_ONE_MINUS_SRC_COLOR");
	else if ( def_blend_isosurf.dfunc == GL_SRC_ALPHA )
	  sprintf(res4,"GL_SRC_ALPHA");
	else if ( def_blend_isosurf.dfunc == GL_ONE_MINUS_SRC_ALPHA )
	  sprintf(res4,"GL_ONE_MINUS_SRC_ALPHA");
	else if ( def_blend_isosurf.dfunc == GL_DST_ALPHA )
	  sprintf(res4,"GL_DST_ALPHA");
	else if ( def_blend_isosurf.dfunc == GL_ONE_MINUS_DST_ALPHA )
	  sprintf(res4,"GL_ONE_MINUS_DST_ALPHA");
	/* NOW PRINT ALL RESx to "result" */
	sprintf(result,"%s %s %s %s", res1, res2, res3, res4);	
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown value \"%s\", must be isosurf, while executing %s %s %s %s", argv[3], argv[0], argv[2], argv[3], argv[4]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      free( (void *) res1);
      free( (void *) res2);
      free( (void *) res3);
      free( (void *) res4);
    }
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss),"unknown option \"%s\", must be -what, while executing %s %s %s %s", argv[3], argv[0], argv[2], argv[3], argv[4]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  /* FOG */
  else if ( strcmp(argv[1],"fog") == 0 ) {
    glparam.what = GLPAR_FOG;
    if ( strcmp(argv[2],"-get") == 0 ) {
      if ( strcmp(argv[3],"fogcolor") == 0 ) {
	fogcolor = 1;
      } 
      else if ( strcmp(argv[3],"def_fogcolor") == 0 ) {
	def_fogcolor = 1;
      } 
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss), "unknown value \"%s\", must be fogcolor or def_fogcolor", argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss), "unknown option %s, must be -get", argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }        

  else {
    char rss[1024];
    snprintf(rss, sizeof(rss),"unknown argument \"%s\", must be one of material, light, lightmodel, blendfunc or fog, while executing %s %s ...",
	     argv[1], argv[0], argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /*************************/
  /* now updata parameters */
  /*************************/
  if ( glparam.what == GLPAR_MATERIAL ) {
    if ( glparam.mat_what == GLPAR_MAT_WHAT_STRUCTURE ) {
      if ( specular )
	sprintf(result,"%f %f %f %f", 
		mat_struct.specular[0], mat_struct.specular[1],
		mat_struct.specular[2], mat_struct.specular[3]);      
      else if ( shininess ) 
	sprintf(result,"%f",mat_struct.shininess[0]);
      else if ( ambient )
	sprintf(result,"%f %f %f %f",
		mat_struct.ambient[0], mat_struct.ambient[1],
		mat_struct.ambient[2], mat_struct.ambient[3]);
      else if ( diffuse ) 
	sprintf(result,"%f %f %f %f", 
		mat_struct.diffuse[0], mat_struct.diffuse[1],
		mat_struct.diffuse[2], mat_struct.diffuse[3]);
      else if ( emission ) 
	sprintf(result,"%f %f %f %f", 
		mat_struct.emission[0], mat_struct.emission[1],
		mat_struct.emission[2], mat_struct.emission[3]);
      /* GET DEFAULTS */
      else if ( def_specular )
	sprintf(result,"%f %f %f %f", 
		def_mat_specular[0], def_mat_specular[1],
		def_mat_specular[2], def_mat_specular[3]);      
      else if ( def_shininess ) 
	sprintf(result,"%f",def_mat_shininess[0]);
      else if ( def_ambient )
	sprintf(result,"%f %f %f %f",
		def_mat_ambient[0], def_mat_ambient[1],
		def_mat_ambient[2], def_mat_ambient[3]);
      else if ( def_diffuse ) 
	sprintf(result,"%f %f %f %f", 
		def_mat_diffuse[0], def_mat_diffuse[1],
		def_mat_diffuse[2], def_mat_diffuse[3]);
      else if ( def_emission ) 
	sprintf(result,"%f %f %f %f", 
		def_mat_emission[0], def_mat_emission[1],
		def_mat_emission[2], def_mat_emission[3]);
      /* END -- GET DEFAULTS */
    }
    else if ( glparam.mat_what == GLPAR_MAT_WHAT_ISOSURF ) {
      int tr;
      tr = isoDisp.transparent;
      if ( glparam.mat_what_isosurf == GLPAR_MAT_WHAT_ISOONE ) {
	if ( glparam.mat_what_isoside == GLPAR_MAT_WHAT_ISOFRONTSIDE ) {
	  if ( specular ) 
	    sprintf(result,"%f %f %f %f", 
		    front_mat_isosurf[tr].specular[0], 
		    front_mat_isosurf[tr].specular[1],
		    front_mat_isosurf[tr].specular[2], 
		    front_mat_isosurf[tr].specular[3]);
	  if ( shininess ) 
	    sprintf(result,"%f", 
		    front_mat_isosurf[tr].shininess[0]);
	  if ( ambient ) 
	    sprintf(result,"%f %f %f %f", 
		    front_mat_isosurf[tr].ambient[0], 
		    front_mat_isosurf[tr].ambient[1],
		    front_mat_isosurf[tr].ambient[2], 
		    front_mat_isosurf[tr].ambient[3]);	    
	  if ( diffuse ) 
	    sprintf(result,"%f %f %f %f", 
		    front_mat_isosurf[tr].diffuse[0], 
		    front_mat_isosurf[tr].diffuse[1],
		    front_mat_isosurf[tr].diffuse[2], 
		    front_mat_isosurf[tr].diffuse[3]);
	  if ( emission ) 
	    sprintf(result,"%f %f %f %f", 
		    front_mat_isosurf[tr].emission[0], 
		    front_mat_isosurf[tr].emission[1],
		    front_mat_isosurf[tr].emission[2], 
		    front_mat_isosurf[tr].emission[3]);
	  /* GET DEFAULTS */
	  if ( tr == 0 ) {
	    /* TRANSPARENCY IS TURNED OFF */
	    if ( def_specular ) 
	      sprintf(result,"%f %f %f %f", 
		      def_front_mat_specular[0], 
		      def_front_mat_specular[1],
		      def_front_mat_specular[2], 
		      def_front_mat_specular[3]);
	    if ( def_shininess ) 
	      sprintf(result,"%f", 
		      def_front_mat_shininess[0]);
	    if ( def_ambient ) 
	      sprintf(result,"%f %f %f %f", 
		      def_front_mat_ambient[0], 
		      def_front_mat_ambient[1],
		      def_front_mat_ambient[2], 
		      def_front_mat_ambient[3]);	    
	    if ( def_diffuse ) 
	      sprintf(result,"%f %f %f %f", 
		      def_front_mat_diffuse[0], 
		      def_front_mat_diffuse[1],
		      def_front_mat_diffuse[2], 
		      def_front_mat_diffuse[3]);
	    if ( def_emission ) 
	      sprintf(result,"%f %f %f %f", 
		      def_front_mat_emission[0], 
		      def_front_mat_emission[1],
		      def_front_mat_emission[2], 
		      def_front_mat_emission[3]);
	  } else {
	    /* TRANSPARENCY IS TURNED OFF */
	    if ( def_specular ) 
	      sprintf(result,"%f %f %f %f", 
		      def_blend_front_mat_specular[0], 
		      def_blend_front_mat_specular[1],
		      def_blend_front_mat_specular[2], 
		      def_blend_front_mat_specular[3]);
	    if ( def_shininess ) 
	      sprintf(result,"%f", 
		      def_blend_front_mat_shininess[0]);
	    if ( def_ambient ) 
	      sprintf(result,"%f %f %f %f", 
		      def_blend_front_mat_ambient[0], 
		      def_blend_front_mat_ambient[1],
		      def_blend_front_mat_ambient[2], 
		      def_blend_front_mat_ambient[3]);	    
	    if ( def_diffuse ) 
	      sprintf(result,"%f %f %f %f", 
		      def_blend_front_mat_diffuse[0], 
		      def_blend_front_mat_diffuse[1],
		      def_blend_front_mat_diffuse[2], 
		      def_blend_front_mat_diffuse[3]);
	    if ( def_emission ) 
	      sprintf(result,"%f %f %f %f", 
		      def_blend_front_mat_emission[0], 
		      def_blend_front_mat_emission[1],
		      def_blend_front_mat_emission[2], 
		      def_blend_front_mat_emission[3]);
	  }
	  /* END -- GET DEFAULTS */
	}
	else if ( glparam.mat_what_isoside == GLPAR_MAT_WHAT_ISOBACKSIDE ) {
	  if ( specular ) 
	    sprintf(result,"%f %f %f %f", 
		    back_mat_isosurf[tr].specular[0], 
		    back_mat_isosurf[tr].specular[1],
		    back_mat_isosurf[tr].specular[2], 
		    back_mat_isosurf[tr].specular[3]);
	  if ( shininess ) 
	    sprintf(result,"%f", back_mat_isosurf[tr].shininess[0]);
	  if ( ambient ) 
	    sprintf(result,"%f %f %f %f", 
		    back_mat_isosurf[tr].ambient[0], 
		    back_mat_isosurf[tr].ambient[1],
		    back_mat_isosurf[tr].ambient[2], 
		    back_mat_isosurf[tr].ambient[3]);
	  if ( diffuse ) 
	    sprintf(result,"%f %f %f %f", 
		    back_mat_isosurf[tr].diffuse[0], 
		    back_mat_isosurf[tr].diffuse[1],
		    back_mat_isosurf[tr].diffuse[2], 
		    back_mat_isosurf[tr].diffuse[3]);
	  if ( emission ) 
	    sprintf(result,"%f %f %f %f", 
		    back_mat_isosurf[tr].emission[0], 
		    back_mat_isosurf[tr].emission[1],
		    back_mat_isosurf[tr].emission[2], 
		    back_mat_isosurf[tr].emission[3]);
	  /* GET DEFAULTS */
	  if ( tr == 0 ) {
	    /* TRANSPARENCY IS TURNED OFF */
	    if ( def_specular ) 
	      sprintf(result,"%f %f %f %f", 
		      def_back_mat_specular[0], 
		      def_back_mat_specular[1],
		      def_back_mat_specular[2], 
		      def_back_mat_specular[3]);
	    if ( def_shininess ) 
	      sprintf(result,"%f", def_back_mat_shininess[0]);
	    if ( def_ambient ) 
	      sprintf(result,"%f %f %f %f", 
		      def_back_mat_ambient[0], 
		      def_back_mat_ambient[1],
		      def_back_mat_ambient[2], 
		      def_back_mat_ambient[3]);
	    if ( def_diffuse ) 
	      sprintf(result,"%f %f %f %f", 
		      def_back_mat_diffuse[0], 
		      def_back_mat_diffuse[1],
		      def_back_mat_diffuse[2], 
		      def_back_mat_diffuse[3]);
	    if ( def_emission ) 
	      sprintf(result,"%f %f %f %f", 
		      def_back_mat_emission[0], 
		      def_back_mat_emission[1],
		      def_back_mat_emission[2], 
		      def_back_mat_emission[3]);
	  } else {
	    /* TRANSPARENCY IS TURNED ON */
	    if ( def_specular ) 
	      sprintf(result,"%f %f %f %f", 
		      def_blend_back_mat_specular[0], 
		      def_blend_back_mat_specular[1],
		      def_blend_back_mat_specular[2], 
		      def_blend_back_mat_specular[3]);
	    if ( def_shininess ) 
	      sprintf(result,"%f", def_blend_back_mat_shininess[0]);
	    if ( def_ambient ) 
	      sprintf(result,"%f %f %f %f", 
		      def_blend_back_mat_ambient[0], 
		      def_blend_back_mat_ambient[1],
		      def_blend_back_mat_ambient[2], 
		      def_blend_back_mat_ambient[3]);
	    if ( def_diffuse ) 
	      sprintf(result,"%f %f %f %f", 
		      def_blend_back_mat_diffuse[0], 
		      def_blend_back_mat_diffuse[1],
		      def_blend_back_mat_diffuse[2], 
		      def_blend_back_mat_diffuse[3]);
	    if ( def_emission ) 
	      sprintf(result,"%f %f %f %f", 
		      def_blend_back_mat_emission[0], 
		      def_blend_back_mat_emission[1],
		      def_blend_back_mat_emission[2], 
		      def_blend_back_mat_emission[3]);
	  }
	  /* END -- GET DEFAULTS */
	}
      }
      else if ( glparam.mat_what_isosurf == GLPAR_MAT_WHAT_ISOPOS ) {
	if ( glparam.mat_what_isoside == GLPAR_MAT_WHAT_ISOFRONTSIDE ) {
	  if ( specular )
	    sprintf(result,"%f %f %f %f", 
		    front_mat_pos_isosurf[tr].specular[0], 
		    front_mat_pos_isosurf[tr].specular[1],
		    front_mat_pos_isosurf[tr].specular[2], 
		    front_mat_pos_isosurf[tr].specular[3]);
	  if ( shininess )
	    sprintf(result,"%f", front_mat_pos_isosurf[tr].shininess[0]); 
	  if ( ambient )
	    sprintf(result,"%f %f %f %f", 
		    front_mat_pos_isosurf[tr].ambient[0], 
		    front_mat_pos_isosurf[tr].ambient[1],
		    front_mat_pos_isosurf[tr].ambient[2], 
		    front_mat_pos_isosurf[tr].ambient[3]);
	  if ( diffuse )
	    sprintf(result,"%f %f %f %f", 
		    front_mat_pos_isosurf[tr].diffuse[0], 
		    front_mat_pos_isosurf[tr].diffuse[1],
		    front_mat_pos_isosurf[tr].diffuse[2], 
		    front_mat_pos_isosurf[tr].diffuse[3]);
	  if ( emission )
	    sprintf(result,"%f %f %f %f", 
		    front_mat_pos_isosurf[tr].emission[0], 
		    front_mat_pos_isosurf[tr].emission[1],
		    front_mat_pos_isosurf[tr].emission[2], 
		    front_mat_pos_isosurf[tr].emission[3]);
	  /* GET DEFAULTS */
	  if ( tr == 0 ) {
	    /* TRANSPARENCY IS TURNED OFF */
	    /* for SPECULAR & SHININESS componenets there are no 
	     * special-designed parameters, use def_front_mat_*
	     */
	    if ( def_specular )
	      sprintf(result,"%f %f %f %f", 
		      def_front_mat_specular[0], 
		      def_front_mat_specular[1],
		      def_front_mat_specular[2], 
		      def_front_mat_specular[3]);
	    if ( def_shininess )
	      sprintf(result,"%f", def_front_mat_shininess[0]); 
	    if ( def_ambient )
	      sprintf(result,"%f %f %f %f", 
		      def_pos_front_mat_ambient[0], 
		      def_pos_front_mat_ambient[1],
		      def_pos_front_mat_ambient[2], 
		      def_pos_front_mat_ambient[3]);
	    if ( def_diffuse )
	      sprintf(result,"%f %f %f %f", 
		      def_pos_front_mat_diffuse[0], 
		      def_pos_front_mat_diffuse[1],
		      def_pos_front_mat_diffuse[2], 
		      def_pos_front_mat_diffuse[3]);
	    if ( def_emission )
	      sprintf(result,"%f %f %f %f", 
		      def_pos_front_mat_emission[0], 
		      def_pos_front_mat_emission[1],
		      def_pos_front_mat_emission[2], 
		      def_pos_front_mat_emission[3]);
	  } else {
	    /* TRANSPARENCY IS TURNED ON */
	    /* for SPECULAR & SHININESS componenets there are no 
	     * special-designed parameters, use def_front_mat_*
	     */
	    if ( def_specular )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_front_mat_specular[0], 
		      def_blend_front_mat_specular[1],
		      def_blend_front_mat_specular[2], 
		      def_blend_front_mat_specular[3]);
	    if ( def_shininess )
	      sprintf(result,"%f", def_blend_front_mat_shininess[0]); 
	    if ( def_ambient )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_pos_front_mat_ambient[0], 
		      def_blend_pos_front_mat_ambient[1],
		      def_blend_pos_front_mat_ambient[2], 
		      def_blend_pos_front_mat_ambient[3]);
	    if ( def_diffuse )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_pos_front_mat_diffuse[0], 
		      def_blend_pos_front_mat_diffuse[1],
		      def_blend_pos_front_mat_diffuse[2], 
		      def_blend_pos_front_mat_diffuse[3]);	  
	    if ( def_emission )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_pos_front_mat_emission[0], 
		      def_blend_pos_front_mat_emission[1],
		      def_blend_pos_front_mat_emission[2], 
		      def_blend_pos_front_mat_emission[3]);	  
	  }
	  /* END -- GET DEFAULTS */
	}
	else if ( glparam.mat_what_isoside == GLPAR_MAT_WHAT_ISOBACKSIDE ) {
	  if ( specular )
	    sprintf(result,"%f %f %f %f", 
		    back_mat_pos_isosurf[tr].specular[0], 
		    back_mat_pos_isosurf[tr].specular[1],
		    back_mat_pos_isosurf[tr].specular[2], 
		    back_mat_pos_isosurf[tr].specular[3]);
	  if ( shininess ) 
	    sprintf(result,"%f", front_mat_pos_isosurf[tr].shininess[0]); 
	  if ( ambient ) 
	    sprintf(result,"%f %f %f %f", 
		    back_mat_pos_isosurf[tr].ambient[0], 
		    back_mat_pos_isosurf[tr].ambient[1],
		    back_mat_pos_isosurf[tr].ambient[2], 
		    back_mat_pos_isosurf[tr].ambient[3]);
	  if ( diffuse )
	    sprintf(result,"%f %f %f %f", 
		    back_mat_pos_isosurf[tr].diffuse[0], 
		    back_mat_pos_isosurf[tr].diffuse[1],
		    back_mat_pos_isosurf[tr].diffuse[2], 
		    back_mat_pos_isosurf[tr].diffuse[3]);
	  if ( emission )
	    sprintf(result,"%f %f %f %f", 
		    back_mat_pos_isosurf[tr].emission[0], 
		    back_mat_pos_isosurf[tr].emission[1],
		    back_mat_pos_isosurf[tr].emission[2], 
		    back_mat_pos_isosurf[tr].emission[3]);
	  /* GET DEFAULTS */
	  if ( tr == 0 ) {
	    /* TRANSPARENCY IS TURNED OFF */
	    /* for SPECULAR & SHININESS componenets there are no 
	     * special-designed parameters, use def_back_mat_*
	     */
	    if ( def_specular )
	      sprintf(result,"%f %f %f %f", 
		      def_back_mat_specular[0], 
		      def_back_mat_specular[1],
		      def_back_mat_specular[2], 
		      def_back_mat_specular[3]);
	    if ( def_shininess )
	      sprintf(result,"%f", def_back_mat_shininess[0]); 
	    if ( def_ambient )
	      sprintf(result,"%f %f %f %f", 
		      def_pos_back_mat_ambient[0], 
		      def_pos_back_mat_ambient[1],
		      def_pos_back_mat_ambient[2], 
		      def_pos_back_mat_ambient[3]);
	    if ( def_diffuse )
	      sprintf(result,"%f %f %f %f", 
		      def_pos_back_mat_diffuse[0], 
		      def_pos_back_mat_diffuse[1],
		      def_pos_back_mat_diffuse[2], 
		      def_pos_back_mat_diffuse[3]);
	    if ( def_emission )
	      sprintf(result,"%f %f %f %f", 
		      def_pos_back_mat_emission[0], 
		      def_pos_back_mat_emission[1],
		      def_pos_back_mat_emission[2], 
		      def_pos_back_mat_emission[3]);
	  } else {
	    /* TRANSPARENCY IS TURNED ON */
	    /* for SPECULAR & SHININESS componenets there are no 
	     * special-designed parameters, use def_back_mat_*
	     */
	    if ( def_specular )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_back_mat_specular[0], 
		      def_blend_back_mat_specular[1],
		      def_blend_back_mat_specular[2], 
		      def_blend_back_mat_specular[3]);
	    if ( def_shininess )
	      sprintf(result,"%f", def_blend_back_mat_shininess[0]); 
	    if ( def_ambient )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_pos_back_mat_ambient[0], 
		      def_blend_pos_back_mat_ambient[1],
		      def_blend_pos_back_mat_ambient[2], 
		      def_blend_pos_back_mat_ambient[3]);
	    if ( def_diffuse )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_pos_back_mat_diffuse[0], 
		      def_blend_pos_back_mat_diffuse[1],
		      def_blend_pos_back_mat_diffuse[2], 
		      def_blend_pos_back_mat_diffuse[3]);	  
	    if ( def_emission )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_pos_back_mat_emission[0], 
		      def_blend_pos_back_mat_emission[1],
		      def_blend_pos_back_mat_emission[2], 
		      def_blend_pos_back_mat_emission[3]);	  
	  }
	  /* END -- GET DEFAULTS */
	}      
      }      
      else if ( glparam.mat_what_isosurf == GLPAR_MAT_WHAT_ISONEG ) {
	if ( glparam.mat_what_isoside == GLPAR_MAT_WHAT_ISOFRONTSIDE ) {
	  if ( specular ) 
	    sprintf(result,"%f %f %f %f", 
		    front_mat_neg_isosurf[tr].specular[0], 
		    front_mat_neg_isosurf[tr].specular[1],
		    front_mat_neg_isosurf[tr].specular[2], 
		    front_mat_neg_isosurf[tr].specular[3]); 
	  if ( shininess )
	    sprintf(result,"%f", front_mat_neg_isosurf[tr].shininess[0]);
	  if ( ambient )
	    sprintf(result,"%f %f %f %f", 
		    front_mat_neg_isosurf[tr].ambient[0], 
		    front_mat_neg_isosurf[tr].ambient[1],
		    front_mat_neg_isosurf[tr].ambient[2], 
		    front_mat_neg_isosurf[tr].ambient[3]); 
	  if ( diffuse )
	    sprintf(result,"%f %f %f %f", 
		    front_mat_neg_isosurf[tr].diffuse[0], 
		    front_mat_neg_isosurf[tr].diffuse[1],
		    front_mat_neg_isosurf[tr].diffuse[2], 
		    front_mat_neg_isosurf[tr].diffuse[3]); 
	  if ( emission )
	    sprintf(result,"%f %f %f %f", 
		    front_mat_neg_isosurf[tr].emission[0], 
		    front_mat_neg_isosurf[tr].emission[1],
		    front_mat_neg_isosurf[tr].emission[2], 
		    front_mat_neg_isosurf[tr].emission[3]); 
	  /* GET DEFAULTS */
	  if ( tr == 0 ) {
	    /* TRANSPARENCY IS TURNED OFF */
	    /* for SPECULAR & SHININESS componenets there are no 
	     * special-designed parameters, use def_front_mat_*
	     */
	    if ( def_specular )
	      sprintf(result,"%f %f %f %f", 
		      def_front_mat_specular[0], 
		      def_front_mat_specular[1],
		      def_front_mat_specular[2], 
		      def_front_mat_specular[3]);
	    if ( def_shininess )
	      sprintf(result,"%f", def_front_mat_shininess[0]); 
	    if ( def_ambient )
	      sprintf(result,"%f %f %f %f", 
		      def_neg_front_mat_ambient[0], 
		      def_neg_front_mat_ambient[1],
		      def_neg_front_mat_ambient[2], 
		      def_neg_front_mat_ambient[3]);
	    if ( def_diffuse )
	      sprintf(result,"%f %f %f %f", 
		      def_neg_front_mat_diffuse[0], 
		      def_neg_front_mat_diffuse[1],
		      def_neg_front_mat_diffuse[2], 
		      def_neg_front_mat_diffuse[3]);
	    if ( def_emission )
	      sprintf(result,"%f %f %f %f", 
		      def_neg_front_mat_emission[0], 
		      def_neg_front_mat_emission[1],
		      def_neg_front_mat_emission[2], 
		      def_neg_front_mat_emission[3]);
	  } else {
	    /* TRANSPARENCY IS TURNED ON */
	    /* for SPECULAR & SHININESS componenets there are no 
	     * special-designed parameters, use def_front_mat_*
	     */
	    if ( def_specular )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_front_mat_specular[0], 
		      def_blend_front_mat_specular[1],
		      def_blend_front_mat_specular[2], 
		      def_blend_front_mat_specular[3]);
	    if ( def_shininess )
	      sprintf(result,"%f", def_blend_front_mat_shininess[0]); 
	    if ( def_ambient )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_neg_front_mat_ambient[0], 
		      def_blend_neg_front_mat_ambient[1],
		      def_blend_neg_front_mat_ambient[2], 
		      def_blend_neg_front_mat_ambient[3]);
	    if ( def_diffuse )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_neg_front_mat_diffuse[0], 
		      def_blend_neg_front_mat_diffuse[1],
		      def_blend_neg_front_mat_diffuse[2], 
		      def_blend_neg_front_mat_diffuse[3]);	  
	    if ( def_emission )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_neg_front_mat_emission[0], 
		      def_blend_neg_front_mat_emission[1],
		      def_blend_neg_front_mat_emission[2], 
		      def_blend_neg_front_mat_emission[3]);	  
	  }
	  /* END -- GET DEFAULTS */
	}
	else if ( glparam.mat_what_isoside == GLPAR_MAT_WHAT_ISOBACKSIDE ) {
	  if ( specular ) 
	    sprintf(result,"%f %f %f %f", 
		    back_mat_neg_isosurf[tr].specular[0], 
		    back_mat_neg_isosurf[tr].specular[1],
		    back_mat_neg_isosurf[tr].specular[2], 
		    back_mat_neg_isosurf[tr].specular[3]); 
	  if ( shininess )
	    sprintf(result,"%f", 
		    back_mat_neg_isosurf[tr].shininess[0]);
	  if ( ambient )
	    sprintf(result,"%f %f %f %f", 
		    back_mat_neg_isosurf[tr].ambient[0], 
		    back_mat_neg_isosurf[tr].ambient[1],
		    back_mat_neg_isosurf[tr].ambient[2], 
		    back_mat_neg_isosurf[tr].ambient[3]); 
	  if ( diffuse )
	    sprintf(result,"%f %f %f %f", 
		    back_mat_neg_isosurf[tr].diffuse[0], 
		    back_mat_neg_isosurf[tr].diffuse[1],
		    back_mat_neg_isosurf[tr].diffuse[2], 
		    back_mat_neg_isosurf[tr].diffuse[3]); 
	  if ( emission )
	    sprintf(result,"%f %f %f %f", 
		    back_mat_neg_isosurf[tr].emission[0], 
		    back_mat_neg_isosurf[tr].emission[1],
		    back_mat_neg_isosurf[tr].emission[2], 
		    back_mat_neg_isosurf[tr].emission[3]); 
	  /* GET DEFAULTS */
	  if ( tr == 0 ) {
	    /* TRANSPARENCY IS TURNED OFF */
	    /* for SPECULAR & SHININESS componenets there are no 
	     * special-designed parameters, use def_back_mat_*
	     */
	    if ( def_specular )
	      sprintf(result,"%f %f %f %f", 
		      def_back_mat_specular[0], 
		      def_back_mat_specular[1],
		      def_back_mat_specular[2], 
		      def_back_mat_specular[3]);
	    if ( def_shininess )
	      sprintf(result,"%f", 
		      def_back_mat_shininess[0]); 
	    if ( def_ambient )
	      sprintf(result,"%f %f %f %f", 
		      def_neg_back_mat_ambient[0], 
		      def_neg_back_mat_ambient[1],
		      def_neg_back_mat_ambient[2], 
		      def_neg_back_mat_ambient[3]);
	    if ( def_diffuse )
	      sprintf(result,"%f %f %f %f", 
		      def_neg_back_mat_diffuse[0], 
		      def_neg_back_mat_diffuse[1],
		      def_neg_back_mat_diffuse[2], 
		      def_neg_back_mat_diffuse[3]);
	    if ( def_emission )
	      sprintf(result,"%f %f %f %f", 
		      def_neg_back_mat_emission[0], 
		      def_neg_back_mat_emission[1],
		      def_neg_back_mat_emission[2], 
		      def_neg_back_mat_emission[3]);
	  } else {
	    /* TRANSPARENCY IS TURNED ON */
	    /* for SPECULAR & SHININESS componenets there are no 
	     * special-designed parameters, use def_back_mat_*
	     */
	    if ( def_specular )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_back_mat_specular[0], 
		      def_blend_back_mat_specular[1],
		      def_blend_back_mat_specular[2], 
		      def_blend_back_mat_specular[3]);
	    if ( def_shininess )
	      sprintf(result,"%f", 
		      def_blend_back_mat_shininess[0]); 
	    if ( def_ambient )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_neg_back_mat_ambient[0], 
		      def_blend_neg_back_mat_ambient[1],
		      def_blend_neg_back_mat_ambient[2], 
		      def_blend_neg_back_mat_ambient[3]);
	    if ( def_diffuse )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_neg_back_mat_diffuse[0], 
		      def_blend_neg_back_mat_diffuse[1],
		      def_blend_neg_back_mat_diffuse[2], 
		      def_blend_neg_back_mat_diffuse[3]);	  
	    if ( def_emission )
	      sprintf(result,"%f %f %f %f", 
		      def_blend_neg_back_mat_emission[0], 
		      def_blend_neg_back_mat_emission[1],
		      def_blend_neg_back_mat_emission[2], 
		      def_blend_neg_back_mat_emission[3]);	  
	  }
	  /* END -- GET DEFAULTS */
	}
      }   
    }
  }
  else if ( glparam.what == GLPAR_LIGHT ) {
    /* if -light option wasn't present all data regard to GL_LIGHT0 */
    int il = 0;
    if ( lightN.is ) {
      il = (int) lightN.vec[0];
      if ( il < 0 || il > GLPAR_MAXLIGHT ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"light number %d is out of range, should be between 0 and %d",il, GLPAR_MAXLIGHT);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    if ( ambient )
      sprintf(result,"%f %f %f %f", 
	      light[il].ambient[0], light[il].ambient[1],
	      light[il].ambient[2], light[il].ambient[3]);
    if ( diffuse )
      sprintf(result,"%f %f %f %f", 
	      light[il].diffuse[0], light[il].diffuse[1],
	      light[il].diffuse[2], light[il].diffuse[3]);
    if ( specular )
      sprintf(result,"%f %f %f %f", 
	      light[il].specular[0], light[il].specular[1],
	      light[il].specular[2], light[il].specular[3]);
    if ( fract_position )
      sprintf(result,"%f %f %f %f", 
	      light[il].fract_position[0], light[il].fract_position[1],
	      light[il].fract_position[2], light[il].fract_position[3]);    
    if ( spot_dir )
      sprintf(result,"%f %f %f", 
	      light[il].spot_dir[0], light[il].spot_dir[1],
	      light[il].spot_dir[2]);
    if ( spot_exp )
      sprintf(result,"%f", light[il].spot_exp[0]);
    if ( spot_cutoff )
      sprintf(result,"%f", light[il].spot_cutoff[0]);
    if ( const_atten )
      sprintf(result,"%f", light[il].const_atten[0]);
    if ( lin_atten ) 
      sprintf(result,"%f", light[il].lin_atten[0]);
    if ( quad_atten )
      sprintf(result,"%f", light[il].quad_atten[0]);
    /* GET DEFAULTS */
    if ( def_ambient )
      sprintf(result,"%f %f %f %f", 
	      def_light_ambient[il][0], def_light_ambient[il][1],
	      def_light_ambient[il][2], def_light_ambient[il][3]);
    if ( def_diffuse )
      sprintf(result,"%f %f %f %f", 
	      def_light_diffuse[il][0], def_light_diffuse[il][1],
	      def_light_diffuse[il][2], def_light_diffuse[il][3]);
    if ( def_specular )
      sprintf(result,"%f %f %f %f", 
	      def_light_specular[il][0], def_light_specular[il][1],
	      def_light_specular[il][2], def_light_specular[il][3]);
    if ( def_fract_position )
      sprintf(result,"%f %f %f %f", 
	      def_light_fract_position[il][0], def_light_fract_position[il][1],
	      def_light_fract_position[il][2], def_light_fract_position[il][3]);    
    if ( def_spot_dir )
      sprintf(result,"%f %f %f", 
	      def_light_spotdir[il][0], def_light_spotdir[il][1],
	      def_light_spotdir[il][2]);
    if ( def_spot_exp )
      sprintf(result,"%f", def_light_spot_exp[il][0]);
    if ( def_spot_cutoff )
      sprintf(result,"%f", def_light_spot_cutoff[il][0]);
    if ( def_const_atten )
      sprintf(result,"%f", def_light_const_atten[il][0]);
    if ( def_lin_atten ) 
      sprintf(result,"%f", def_light_lin_atten[il][0]);
    if ( def_quad_atten )
      sprintf(result,"%f", def_light_quad_atten[il][0]);
  }

  else if ( glparam.what == GLPAR_FOG ) {
    if ( fogcolor ) {
      sprintf(result,"%f %f %f %f", 
	      fog.color[0], fog.color[1], fog.color[2], fog.color[3]);
    }
    else if ( def_fogcolor ) {
      sprintf(result,"%f %f %f %f",
	      def_fog_color[0], def_fog_color[1],
	      def_fog_color[2], def_fog_color[3]);
    }
  }

  /* now set the result */
  Tcl_SetResult(interp, result, TCL_DYNAMIC);

  return TCL_OK;
}



int
xcSplitList( int code, Tcl_Interp *interp, const char *argv[], GetGlParam *var) 
{
  int i, inte, imax = 1;
  int argcList;
  const char **argvList;
  double val[20];

  if ( code == XC_GET_RGBA )    imax = 4;
  if ( code == XC_GET_XYZW )    imax = 4;
  if ( code == XC_GET_XYZ )     imax = 3;
  if ( code == XC_GET_RGB )     imax = 3;
  if ( code == XC_GET_FLOAT )   imax = 1;
  if ( code == XC_GET_FLOAT2 )  imax = 2;
  if ( code == XC_GET_FLOAT4 )  imax = 4;
  if ( code == XC_GET_FLOAT6 )  imax = 6;
  if ( code == XC_GET_FLOAT12 ) imax = 12;
  if ( code == XC_GET_INT )     imax = 1;
  if ( code == XC_GET_INT3 )    imax = 3;
  if ( code == XC_GET_INT12 )   imax = 12;

  Tcl_SplitList(interp, argv[0], &argcList, &argvList);
  if ( argcList != imax ) {
    xcTclListError( code, interp, argv );
    return XC_ERROR;
  }
  if ( code != XC_GET_INT &&
       code != XC_GET_INT3 &&
       code != XC_GET_INT12 ) {
    for (i=0; i<imax; i++) {
      if ( Tcl_GetDouble(interp, argvList[i], &(val[i])) == TCL_ERROR ) {
	xcTclListError( XC_GET_FLOAT, interp, argv );
	return XC_ERROR;
      }
      var->vec[i] = (float) val[i];
    }
  } else {
    for (i=0; i<imax; i++) {
      if ( Tcl_GetInt(interp, argvList[i], &inte) == TCL_ERROR ) {
	xcTclListError( XC_GET_INT, interp, argv );
	return XC_ERROR;
      }
      var->vec[i] = (float) inte;
    }
  }
  /* var->vec has been assigned, so put vec->is to 1 */
  var->is = 1;

  return XC_OK;
}


void
xcTclListError( int code, Tcl_Interp *interp, const char **argv )
{
  if ( code == XC_GET_RGBA ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"invalid format of RGBA list, should be {r g b a}; while evaluating option %s %s", argv[0], argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
  }

  if ( code == XC_GET_XYZ ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"invalid format of XYZ list, should be {x y z}; while evaluating option %s %s", argv[0], argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
  }

  if ( code == XC_GET_XYZW ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"invalid format of XYZW list, should be {x y z w}; while evaluating option %s %s", argv[0], argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
  }

  if ( code == XC_GET_FLOAT ) {    
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted float, but got %s, while evaluating %s %s",
	     argv[1], argv[0], argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
  }

  if ( code == XC_GET_INT ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, \
but got %s, while evaluating %s %s", argv[1], argv[0], argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
  }
}


static void
CopyGLParVec(int code, float *mat, GetGlParam par)
{
  int i, imax = 1;

  if ( code == XC_GET_RGBA )   imax = 4;
  if ( code == XC_GET_XYZ )    imax = 3;
  if ( code == XC_GET_XYZW)    imax = 4;
  if ( code == XC_GET_FLOAT )  imax = 1;
  if ( code == XC_GET_INT )    imax = 1;
  if ( code == XC_GET_INT3 )   imax = 3;
  if ( code == XC_GET_INT12 )  imax = 12;

  for (i=0; i<imax; i++)
    mat[i] = par.vec[i];
}
