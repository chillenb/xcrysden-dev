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
 * Source: $XCRYSDEN_TOPDIR/C/xcGLparam.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

/* this is preprocessor option tomark that we have include this file */
#define XC_CPP_GLPARAM

#ifndef XC_CPP_GLH
#   include <GL/gl.h>
#   define  XC_CPP_GLH
#endif

/***********************************************/
#define GLPAR_MATERIAL                 0
#define GLPAR_LIGHT                    1
#define GLPAR_LIGHTMODEL               3
#define GLPAR_FRONTFACE                4
#define GLPAR_BLENDFUNC                5
#define GLPAR_FOG                      6

#define GLPAR_MAT_WHAT_STRUCTURE       0
#define GLPAR_MAT_WHAT_ISOSURF         1

#define GLPAR_MAT_WHAT_ISOPOS          0
#define GLPAR_MAT_WHAT_ISONEG          1
#define GLPAR_MAT_WHAT_ISOONE          3

#define GLPAR_MAT_WHAT_ISOBACKSIDE     0
#define GLPAR_MAT_WHAT_ISOFRONTSIDE    1

/* please, do not change these two values */
#define GLPAR_FACE_WHAT_ISOSURF_POS    0
#define GLPAR_FACE_WHAT_ISOSURF_NEG    1

#define GLPAR_BLEND_WHAT_ISOSURF       0
#define GLPAR_BLEND_WHAT_COLORPLANE    1

#define GLPAR_MAXLIGHT                 5
/***********************************************/

/* this is for displaying vorious vectors and cages */
#define GLPAR_N_OF_VEC_AND_CAGE 8
#define GLPAR_PRIMCAGE    0
#define GLPAR_CONVCAGE    1
#define GLPAR_SCCAGE      2
#define GLPAR_PRIMVEC     3
#define GLPAR_CONVVEC     4
#define GLPAR_SCVEC       5
#define GLPAR_ISOPARALLE  6
#define GLPAR_ISOPARAPIP  7

#define XC_GET_RGBA    0
#define XC_GET_XYZ     1
#define XC_GET_XYZW    2
#define XC_GET_FLOAT   3
#define XC_GET_FLOAT4  4 
#define XC_GET_FLOAT6  5
#define XC_GET_FLOAT12 6
#define XC_GET_INT     7
#define XC_GET_INT3    8
#define XC_GET_INT12   9
#define XC_GET_FLOAT2  10
#define XC_GET_RGB     11

typedef struct {
  float ambient[4];
  float diffuse[4];
  float specular[4];
  float shininess[1];
  float emission[4];
} xcMATERIAL;

#ifndef XC_CPP_BOOLEAN
   typedef int boolean;
#  define XC_CPP_BOOLEAN
#endif

typedef struct {
  boolean isenabled;
  float   ambient[4];
  float   diffuse[4];
  float   specular[4];
  float   fract_position[4];  /* position in fractional coordinates */
  float   position[4];
  float   spot_dir[3];
  float   spot_exp[1];
  float   spot_cutoff[1];
  float   const_atten[1];
  float   lin_atten[1];
  float   quad_atten[1];
} LIGHT;

#define LIGHTMODEL_STRUCTURE   0 /* this is used by GL_LIGHT_MODEL_TWO_SIDE */
#define LIGHTMODEL_ISOSURF     1
typedef struct {
  float   two_side[1];              /* two_side for structure            */
  float   two_side_iso[1];          /* two_side for isosufrace & related */
  float   ambient[4];
  float   local_viewer[1];
} LIGHTMODEL;

typedef struct {
  GLenum sfunc;
  GLenum dfunc;
} BLENDFUNC;

typedef GLenum FRONTFACE;

xcMATERIAL mat_struct, mat_colorplane[2];
xcMATERIAL back_mat_isosurf[2], front_mat_isosurf[2];
xcMATERIAL back_mat_pos_isosurf[2], front_mat_pos_isosurf[2];
xcMATERIAL back_mat_neg_isosurf[2], front_mat_neg_isosurf[2];
xcMATERIAL mat_voronoi[2]; 
xcMATERIAL mat_vec_and_cage[GLPAR_N_OF_VEC_AND_CAGE];

#define MAX_NUMBER_OF_LIGHTS      GLPAR_MAXLIGHT+1
LIGHT light[MAX_NUMBER_OF_LIGHTS];
LIGHTMODEL lightmodel;

BLENDFUNC blend_isosurf, blend_colorplane, blend_voronoi, blend_cellcage;
FRONTFACE frontface_isosurf[2];

typedef struct {
  boolean is;
  float vec[20];
} GetGlParam;

