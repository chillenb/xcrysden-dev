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
 * Source: $XCRYSDEN_TOPDIR/C/lighting.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#define  XC_CPP_NO_STDIO    /* look in struct.h file */
#include <GL/gl.h>
#include "struct.h"
#include "xcGLparam.h"
#include "isosurf.h"

void CalcLightPosition (int il, GLdouble size);
void LoadLights(void);
void LoadLightModelTwoSide(int type);
void LoadBlendfunc_And_Frontface(void);
void LoadStructMaterial(void);
void LoadIsoMaterial(int type);
void LoadVoronoiMaterial(int voronoi_type);
void LoadCageOrVecMaterial(int type);

extern ISO_ATTRIB isoDisp;
extern OrthoProj ort;
extern int togl_exists;

/* ================================================ */
/* EMITTED light == barva, ki izvira iz predmeta; neodvisna od barvnega 
   izvora */
/* AMBIENT light == se nakjucno odbija po prostoru */
/* DIFFUSE light == ce neka svetloba pade na nek predmet in se razprsi v 
   vse smeri potem je ta razprsena svetloba DIFFUSE light */
/* SPECULAR light == pride iz dolocene smeri in se odbije, 
   npr. bela pikica na kroglici */

GLenum _LIGHT[8] = { GL_LIGHT0,
		     GL_LIGHT1,
		     GL_LIGHT2,
		     GL_LIGHT3,
		     GL_LIGHT4,
		     GL_LIGHT5,
		     GL_LIGHT6,
		     GL_LIGHT7 };

/*****************************************************************************/
/* this is for structure's material                                          */
GLfloat def_mat_shininess[] = {128.0};
GLfloat def_mat_specular[]  = {0.97, 0.97, 0.97, 1.0};
GLfloat def_mat_ambient[]   = {0.3, 0.3, 0.3, 1.0};
GLfloat def_mat_diffuse[]   = {0.8, 0.8, 0.8, 1.0};
/*GLfloat def_mat_emission[]  = {0.05, 0.05, 0.08, 1.0};*/
GLfloat def_mat_emission[]  = {0.0, 0.0, 0.1, 1.0};

/*****************************************************************************/
/* this is for colorplanes's material                                        */
GLfloat def_cpl_shininess[] = {0.0};
GLfloat def_cpl_specular[]  = {0.0, 0.0, 0.0, 1.0};
GLfloat def_cpl_ambient[]   = {0.2, 0.2, 0.2, 1.0};
GLfloat def_cpl_diffuse[]   = {0.8, 0.8, 0.8, 1.0};

GLfloat def_blend_cpl_shininess[] = {0.0};
GLfloat def_blend_cpl_specular[]  = {0.0, 0.0, 0.0, 1.0};
GLfloat def_blend_cpl_ambient[]   = {0.2, 0.2, 0.2, .6};
GLfloat def_blend_cpl_diffuse[]   = {0.8, 0.8, 0.8, .6};


/***************************************************************************
 * for ISOSURF_TRANSP_ON                                                   */
float def_blend_front_mat_shininess[] = {30.0}; 
float def_blend_front_mat_specular[]  = {1.0, 1.0, 1.0, 1.0}; 
/* s tem kontroliram prednjo barvo, manjse vrednoti ->  */
float def_blend_front_mat_diffuse[]   = {0.7, 0.8, 0.7, 0.6};
float def_blend_front_mat_ambient[]   = {0.5, 0.6, 0.5, 1.0};
float def_blend_front_mat_emission[]  = {0.0, 0.0, 0.0, 1.0};

float def_blend_pos_front_mat_diffuse[] = {1.0, 0.3, 0.3, 0.6};
float def_blend_pos_front_mat_ambient[] = {0.7, 0.5, 0.5, 1.0};
float def_blend_pos_front_mat_emission[]= {0.0, 0.0, 0.0, 1.0};
float def_blend_neg_front_mat_diffuse[] = {0.3, 0.3, 1.0, 0.6}; 
float def_blend_neg_front_mat_ambient[] = {0.5, 0.5, 0.6, 1.0};
float def_blend_neg_front_mat_emission[]= {0.0, 0.0, 0.0, 1.0};

/* s tem kotroliram koliko se vidi steklo, ce za njim ni nobenega */
float def_blend_back_mat_shininess[] = {128.0};
float def_blend_back_mat_specular[]  = {1.0, 1.0, 1.0, 1.0};
/* s tem kontroliram, koliko se vidi zadnjo stran stekla */
float def_blend_back_mat_diffuse[]   = {0.7, .7, .7, .5};
float def_blend_back_mat_ambient[]   = { .7, .7, .7, 1.0};
float def_blend_back_mat_emission[]  = {0.0, 0.0, 0.0, 1.0};

float def_blend_pos_back_mat_diffuse[] = {0.40, .15, .15, 0.5}; 
float def_blend_pos_back_mat_ambient[] = { .25, .20, .20, 1.0};
float def_blend_pos_back_mat_emission[]= { .00, .00, .00, 1.0};
float def_blend_neg_back_mat_diffuse[] = {0.15, .15, .40, 0.5};  
float def_blend_neg_back_mat_ambient[] = {0.20, .20, .25, 1.0};
float def_blend_neg_back_mat_emission[]= { .00, .00, .00, 1.0};


/***************************************************************************
 * for ISOSURF_TRANSP_OFF                                                  */
float def_front_mat_shininess[] = {104.0}; 
float def_front_mat_specular[] = {0.8, 0.8, 0.8, 1.0}; 
/*float def_front_mat_diffuse[] = {1.0, 1.0, 1.0, 1.0};*/
float def_front_mat_diffuse[] = {1.0, 1.0, 0.2, 1.0};
float def_front_mat_ambient[] = {1.0, 1.0, 1.0, 1.0};
float def_front_mat_emission[]= {0.0, 0.0, 0.0, 1.0};

float def_pos_front_mat_diffuse[] = {1.0, 0.3, 0.3, 1.0};
float def_pos_front_mat_ambient[] = {0.7, 0.5, 0.5, 1.0};
float def_pos_front_mat_emission[]= {0.0, 0.0, 0.0, 1.0};
float def_neg_front_mat_diffuse[] = {0.3, 0.3, 1.0, 1.0}; 
float def_neg_front_mat_ambient[] = {0.5, 0.5, 0.6, 1.0};
float def_neg_front_mat_emission[]= {0.0, 0.0, 0.0, 1.0};

/* s tem kotroliram koliko se vidi steklo, ce za njim ni nobenega */
float def_back_mat_shininess[] = {32.0};
float def_back_mat_specular[]  = {0.2, 0.2, 0.2, 1.0};
/*float def_back_mat_diffuse[] = {0.347, .347, .347, 1.0};*/
float def_back_mat_diffuse[]   = {1.0, 0.8, 0.1, 1.0};
float def_back_mat_ambient[]   = { .8, .8, .8, 1.0};
float def_back_mat_emission[]  = {0.0, 0.0, 0.0, 1.0};

float def_pos_back_mat_diffuse[] = {0.40, .15, .15, 0.5}; 
float def_pos_back_mat_ambient[] = { .25, .20, .20, 1.0};
float def_pos_back_mat_emission[]= { .00, .00, .00, 1.0};
float def_neg_back_mat_diffuse[] = {0.15, .15, .40, 0.5};  
float def_neg_back_mat_ambient[] = {0.20, .20, .25, 1.0};
float def_neg_back_mat_emission[]= { .00, .00, .00, 1.0};


/*****************************************************************************/
/* this is for voronoi's material                                            */
GLfloat def_voronoi_shininess[] = {94.0};
GLfloat def_voronoi_specular[]  = {0.40, 0.40, 0.40, 1.0};
GLfloat def_voronoi_ambient[]   = {0.00, 0.95, 0.95, 0.8};
GLfloat def_voronoi_diffuse[]   = {0.00, 0.95, 0.95, 0.8};

GLfloat def_wigner_shininess[] = {94.0};
GLfloat def_wigner_specular[]  = {0.40, 0.40, 0.40, 1.0};
GLfloat def_wigner_ambient[]   = {0.00, 0.95, 0.95, 0.8};
GLfloat def_wigner_diffuse[]   = {0.00, 0.95, 0.95, 0.8};

GLfloat def_bz_shininess[] = {94.0};
GLfloat def_bz_specular[]  = {0.40, 0.40, 0.40, 1.0};
GLfloat def_bz_ambient[]   = {0.00, 0.95, 0.95, 0.4};
GLfloat def_bz_diffuse[]   = {0.00, 0.95, 0.95, 0.4};


/*****************************************************************************/
/* this is for various vectors & cages                                       */
GLfloat def_primcage_shininess[] = {1.0};
GLfloat def_primcage_specular[]  = {0.00, 0.00, 0.00, 1.0};
GLfloat def_primcage_ambient[]   = {1.00, 1.00, 1.00, 1.0};
GLfloat def_primcage_diffuse[]   = {0.00, 1.00, 0.37, 0.43};

GLfloat def_convcage_shininess[] = {1.0};
GLfloat def_convcage_specular[]  = {0.00, 0.00, 0.00, 1.0};
GLfloat def_convcage_ambient[]   = {1.00, 1.00, 1.00, 1.0};
GLfloat def_convcage_diffuse[]   = {1.00, 0.00, 0.37, 0.43};

GLfloat def_sccage_shininess[]   = {1.0};
GLfloat def_sccage_specular[]    = {0.00, 0.00, 0.00, 1.0};
GLfloat def_sccage_ambient[]     = {1.00, 1.00, 1.00, 1.0};
GLfloat def_sccage_diffuse[]     = {0.00, 0.37, 1.00, 0.43};

GLfloat def_primvec_shininess[] = {1.0};
GLfloat def_primvec_specular[]  = {0.00, 0.00, 0.00, 1.0};
GLfloat def_primvec_ambient[]   = {1.00, 1.00, 1.00, 1.0};
GLfloat def_primvec_diffuse[]   = {0.00, 1.00, 0.37, 1.0};

GLfloat def_convvec_shininess[] = {1.0};
GLfloat def_convvec_specular[]  = {0.00, 0.00, 0.00, 1.0};
GLfloat def_convvec_ambient[]   = {1.00, 1.00, 1.00, 1.0};
GLfloat def_convvec_diffuse[]   = {1.00, 0.00, 0.37, 1.0};

GLfloat def_scvec_shininess[] = {1.0};
GLfloat def_scvec_specular[]  = {0.00, 0.00, 0.00, 1.0};
GLfloat def_scvec_ambient[]   = {1.00, 1.00, 1.00, 1.0};
GLfloat def_scvec_diffuse[]   = {0.00, 0.37, 1.00, 1.0};

/*
  for parallelogram IsoSpaceSel
*/
GLfloat def_isoparalle_shininess[]   = {1.0};
GLfloat def_isoparalle_specular[]    = {0.00, 0.00, 0.00, 0.8};
GLfloat def_isoparalle_ambient[]     = {0.80, 0.80, 0.90, 0.8};
GLfloat def_isoparalle_diffuse[]     = {0.80, 0.80, 1.00, 0.6};

/*****************************************************************************/
/* LIGHTS                                                                    */
#define MNOL MAX_NUMBER_OF_LIGHTS     /* MAX_NUMBER_OF_LIGHTS IS 6 */
boolean def_light_isenabled[MNOL]      = { 1, 1, 0, 0, 0, 0 };

/* 
XCRYSDEN: 0.7.2 settings 
------------------------
GLfloat def_light_diffuse[MNOL][4]     = { {0.5, 0.5, 0.6, 1.0},
					   {0.1, 0.1, 0.1, 0.5},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0} };
GLfloat def_light_specular[MNOL][4]    = { {0.8, 0.8, 0.8, 1.0},
					   {0.2, 0.2, 0.2, 0.5},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0} };
GLfloat def_light_fract_position[MNOL][4]={{  5.0,   5.0,  3.0, 0.01},
					   { -5.0,  -5.0, -1.0, 0.01},
					   {  0.0,   0.0,  1.0, 0.0},
					   {  0.0,   0.0,  1.0, 0.0},
					   {  0.0,   0.0,  1.0, 0.0},
					   {  0.0,   0.0,  1.0, 0.0} };
 XCRYSDEN_END */

/*
XCRYSDEN pre 0.8.0 settings:
----------------------------
GLfloat def_light_ambient[MNOL][4]     = { {0.2, 0.2, 0.2, 1.0},
					   {0.01, 0.01, 0.01, 0.5},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0} };

GLfloat def_light_diffuse[MNOL][4]     = { {0.45, 0.45, 0.56, 1.0},
					   {0.1, 0.1, 0.1, 0.5},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0} };
GLfloat def_light_specular[MNOL][4]    = { {0.8, 0.8, 0.8, 1.0},
					   {0.3, 0.3, 0.3, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0} };
GLfloat def_light_fract_position[MNOL][4]={{  5.0,   5.0,  6.0, 0.01},
					   { -5.0,  -5.0, -1.0, 0.01},
					   {  0.0,   0.0,  1.0, 0.0},
					   {  0.0,   0.0,  1.0, 0.0},
					   {  0.0,   0.0,  1.0, 0.0},
					   {  0.0,   0.0,  1.0, 0.0} };
GLfloat def_light_spotdir[MNOL][3]     = { { -5.0,  -5.0,  -3.0},
					   {  5.0,   5.0,   0.5},
					   {  0.0,   0.0,  -1.0},
					   {  0.0,   0.0,  -1.0},
					   {  0.0,   0.0,  -1.0},
					   {  0.0,   0.0,  -1.0} };
*/

GLfloat def_light_ambient[MNOL][4]     = { {0.01, 0.01, 0.01, 1.0},
					   {0.01, 0.01, 0.01, 0.5},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0} };

GLfloat def_light_diffuse[MNOL][4]     = { {0.76, 0.76, 0.76, 1.0},
					   {0.76, 0.76, 0.76, 0.5},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0} };
GLfloat def_light_specular[MNOL][4]    = { {0.6, 0.6, 0.6, 1.0},
					   {0.3, 0.3, 0.3, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0},
					   {0.0, 0.0, 0.0, 1.0} };
/* before light1 was: 5.0,   5.0,  3.0 */
GLfloat def_light_fract_position[MNOL][4]={{  5.0,   8.0,  6.0, 0.00},
					   { -5.0,  -1.0, -1.6, 0.00},
					   {  0.0,   0.0,  1.0, 1.0},
					   {  0.0,   0.0,  1.0, 0.0},
					   {  0.0,   0.0,  1.0, 0.0},
					   {  0.0,   0.0,  1.0, 0.0} };
GLfloat def_light_spotdir[MNOL][3]     = { { -5.0,  -8.0,  -6.0},
					   {  5.0,   1.0,   1.6},
					   {  0.0,   0.0,  -1.0},
					   {  0.0,   0.0,  -1.0},
					   {  0.0,   0.0,  -1.0},
					   {  0.0,   0.0,  -1.0} };
GLfloat def_light_spot_exp[MNOL][1]    = { {128.0}, {115.0}, {10.0},
					   {0.0}, {0.0}, {0.0} };
GLfloat def_light_spot_cutoff[MNOL][1] = { {180.0}, {180.0}, {180.0},
					   {180.0}, {180.0},{180.0} };
GLfloat def_light_const_atten[MNOL][1] = { {1.0}, {1.0}, {1.0},
					   {1.0}, {1.0}, {1.0} };
GLfloat def_light_lin_atten[MNOL][1]   = { {0.0}, {0.0}, {0.0},
					   {0.0}, {0.0}, {0.0} };
GLfloat def_light_quad_atten[MNOL][1]  = { {0.0}, {0.0}, {0.0},
					   {0.0}, {0.0}, {0.0} };


/*****************************************************************************/
/*   LIGHTMODEL --- LIGHTMODEL --- LIGHTMODEL --- LIGHTMODEL --- LIGHTMODEL  */
/*                           this are OpenGL defaults                        */
GLfloat def_lightmodel_two_side[1]     = {0.0};
GLfloat def_lightmodel_two_side_iso[1] = {0.0};
GLfloat def_lightmodel_ambient[4]      = {0.2, 0.2, 0.2, 1.0}; 
GLfloat def_lightmodel_local_viewer[1] = {0.0};


/*****************************************************************************/
/*   FRONTFACE & BLENDFUNC   ---   FRONTFACE & BLENDFUNC                     */
FRONTFACE def_frontface_isosurf = GL_CCW;
BLENDFUNC def_blend_isosurf    = { GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA };
BLENDFUNC def_blend_colorplane = { GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA };
BLENDFUNC def_blend_voronoi    = { GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA };
BLENDFUNC def_blend_cellcage   = { GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA };

/*****************************************************************************/

void CalcLightPosition (int il, GLdouble size) 
{
  GLfloat _spot_dir[3];

  if ( ! togl_exists ) return;

  light[il].position[0] = size * light[il].fract_position[0];
  light[il].position[1] = size * light[il].fract_position[1];
  light[il].position[2] = size * light[il].fract_position[2];
  light[il].position[3] = light[il].fract_position[3];
  glLightfv(_LIGHT[il], GL_POSITION, light[il].position);
  
  _spot_dir[0] = size * light[il].spot_dir[0];
  _spot_dir[1] = size * light[il].spot_dir[1];
  _spot_dir[2] = size * light[il].spot_dir[2];
  glLightfv(_LIGHT[il], GL_SPOT_DIRECTION, _spot_dir);

  /*
	fprintf(stderr,"DEBUG> Light No.%d: position = (%8.3f,%8.3f,%8.3f,%8.3f)\n", il,
		light[il].position[0], light[il].position[1],
		light[il].position[2], light[il].position[3]);
	fprintf(stderr,"DEBUG> Light No.%d: spot_dir = (%8.3f,%8.3f,%8.3f)\n", il, 
		_spot_dir[0], _spot_dir[1], _spot_dir[2]);
  */
}

void
LoadLights(void)
{
  register int i, il;
  static int first_time = 1;
  GLdouble size;

  if ( ! togl_exists ) return;

  /* if I will add more default light parameters, I need to change a little
     the following loop */
  if ( first_time ) {
    /* disable all but the first light */
    for (il=0; il<MAX_NUMBER_OF_LIGHTS; il++) {
      light[il].isenabled      = def_light_isenabled[il];
      light[il].spot_exp[0]    = def_light_spot_exp[il][0];
      light[il].spot_cutoff[0] = def_light_spot_cutoff[il][0];
      light[il].const_atten[0] = def_light_const_atten[il][0];
      light[il].lin_atten[0]   = def_light_lin_atten[il][0];
      light[il].quad_atten[0]  = def_light_quad_atten[il][0];
      for (i=0; i<3; i++) 
	light[il].spot_dir[i]  = def_light_spotdir[il][i];
      for (i=0; i<4; i++) {
	light[il].ambient[i]   = def_light_ambient[il][i];
	light[il].diffuse[i]   = def_light_diffuse[il][i];
	light[il].specular[i]  = def_light_specular[il][i];
	light[il].fract_position[i]  = def_light_fract_position[il][i];
      }      
    }
    /* enable first light */
    light[0].isenabled = 1;
    /* THIS IS TEMPO */

    /* load default lightmodel values */
    lightmodel.two_side[0]      = def_lightmodel_two_side[0];
    lightmodel.two_side_iso[0]  = def_lightmodel_two_side_iso[0];
    lightmodel.ambient[0]       = def_lightmodel_ambient[0]; 
    lightmodel.ambient[1]       = def_lightmodel_ambient[1];
    lightmodel.ambient[2]       = def_lightmodel_ambient[2];
    lightmodel.ambient[3]       = def_lightmodel_ambient[3];
    lightmodel.local_viewer[0]  = def_lightmodel_local_viewer[0];
    first_time = 0;
  }
  
  /* 
   *  this is instead ort.maxz::
   *  so far maxx == maxy, but in future this maight be changed, so:
   */
/*   size = ort.maxx; */
/*   if ( ort.maxy > size ) */
/*     size = ort.maxy; */

  size = ort.size;

  /* take care of lights */
  for ( il=0; il<GLPAR_MAXLIGHT; il++) {
    if ( light[il].isenabled ) {

      /*fprintf(stderr,"\nDEBUG>---------------glLight*** BEGIN; struct.size = %f\n", size);*/

      glLightfv(_LIGHT[il], GL_AMBIENT,  light[il].ambient);
      glLightfv(_LIGHT[il], GL_DIFFUSE,  light[il].diffuse);
      glLightfv(_LIGHT[il], GL_SPECULAR, light[il].specular);
      if ( VPf.stropened ) {
	/*GLfloat _spot_dir[3];*/
	/* */
	/*light[il].position[0] = ort.maxx * light[il].fract_position[0]; */
	/*light[il].position[1] = ort.maxy * light[il].fract_position[1]; */
	/*light[il].position[2] = size     * light[il].fract_position[2]; */
	/*light[il].position[3] = light[il].fract_position[3]; */
	/*glLightfv(_LIGHT[il], GL_POSITION, light[il].position); */
	
	/*_spot_dir[0] = ort.maxx * light[il].spot_dir[0]; */
	/*_spot_dir[1] = ort.maxy * light[il].spot_dir[1]; */
	/*_spot_dir[2] = size     * light[il].spot_dir[2]; */

	CalcLightPosition (il, size);
      }
      
      glLightfv(_LIGHT[il], GL_SPOT_EXPONENT,  light[il].spot_exp);
      glLightfv(_LIGHT[il], GL_SPOT_CUTOFF,    light[il].spot_cutoff);

      glLightfv(_LIGHT[il], GL_CONSTANT_ATTENUATION, light[il].const_atten );
      glLightfv(_LIGHT[il], GL_LINEAR_ATTENUATION, light[il].lin_atten );
      glLightfv(_LIGHT[il], GL_QUADRATIC_ATTENUATION, light[il].quad_atten );

      /*if ( light[il].spot_exp[0] != 0.0 ) { */
      /*glLightfv(_LIGHT[il], GL_SPOT_EXPONENT, light[il].spot_exp); */
      /*} */
      /*if ( light[il].const_atten[0] != 1.0 )*/
      /*if ( light[il].lin_atten[0] != 0.0 ) */
      /*if ( light[il].quad_atten[0] != 0.0 ) */	

      /*
      fprintf(stderr,"DEBUG> Light No.%d: spot_exp = (%8.3f)\n", il, 
	      light[il].spot_exp[0]);
      fprintf(stderr,"DEBUG> Light No.%d: spot_cutoff = (%8.3f)\n", il, 
	      light[il].spot_cutoff[0]);
      fprintf(stderr,"DEBUG> Light No.%d: const_atten = (%8.3f)\n", il, 
	      light[il].const_atten[0]);
      fprintf(stderr,"DEBUG> Light No.%d: lin_atten = (%8.3f)\n", il, 
	      light[il].lin_atten[0]);
      fprintf(stderr,"DEBUG> Light No.%d: quad_atten = (%8.3f)\n", il, 
	      light[il].quad_atten[0]);
      fprintf(stderr,"DEBUG>---------------glLight*** END\n\n");
      */

      /* now enable light */
      if (!glIsEnabled(_LIGHT[il])) glEnable( _LIGHT[il] );
    } else {
      /* disable light */
      glDisable( _LIGHT[il] );
    }
  }

  /* take care for whole ModelLight but GL_LIGHT_MODEL_TWO_SIDE */
  glLightModelfv( GL_LIGHT_MODEL_AMBIENT, lightmodel.ambient );
  glLightModelfv( GL_LIGHT_MODEL_LOCAL_VIEWER, lightmodel.local_viewer );
}



void
LoadLightModelTwoSide(int type)
{
  if ( ! togl_exists ) return;

  if ( type == LIGHTMODEL_STRUCTURE )
    glLightModelfv( GL_LIGHT_MODEL_TWO_SIDE, lightmodel.two_side );
  if ( type == LIGHTMODEL_ISOSURF )
    glLightModelfv( GL_LIGHT_MODEL_TWO_SIDE, lightmodel.two_side_iso );
}



void
LoadBlendfunc_And_Frontface(void)
{
  static int first_time = 1;

  if ( first_time ) {
    /* BLEND_FUNCTION */
    blend_isosurf.sfunc = def_blend_isosurf.sfunc;
    blend_isosurf.dfunc = def_blend_isosurf.dfunc;
    
    blend_colorplane.sfunc = def_blend_colorplane.sfunc;
    blend_colorplane.dfunc = def_blend_colorplane.dfunc;
    
    blend_voronoi.sfunc = def_blend_voronoi.sfunc;
    blend_voronoi.dfunc = def_blend_voronoi.dfunc;

    blend_cellcage.sfunc = def_blend_voronoi.sfunc;
    blend_cellcage.dfunc = def_blend_voronoi.dfunc;

    /* FRONT_FACE */
    frontface_isosurf[0] = def_frontface_isosurf;
    frontface_isosurf[1] = def_frontface_isosurf;

    first_time = 0;
  }

  /* insert code here in the future (MAYBE) */
}



/*****************************************************************************/
/* This func. load a material properties for structure;                      */
void
LoadStructMaterial(void)
{
  static int first_time = 1;
  int i;
  /* EMITTED light == barva, ki izvira iz predmeta */
  /* AMBIENT light == se nakjucno odbija po prostoru */
  /* DIFFUSE light == ce neka svetloba pade na nek predmet in se razprsi v 
     vse smeri potem je ta razprsena svetloba DIFFUSE light */
  /* SPECULAR light == pride iz dolocene smeri in se odbije, 
     npr. bela pikica na kroglici */

  /* ALL ATOMS HAVE THE SAME specular & shininess MATTERIAL PROPERTIES */

  if ( ! togl_exists ) return;
  
  if ( first_time ) {
    mat_struct.shininess[0] = def_mat_shininess[0];
    for (i=0; i<4; i++) {
      mat_struct.specular[i] = def_mat_specular[i]; 
      mat_struct.ambient[i]  = def_mat_ambient[i];
      mat_struct.diffuse[i]  = def_mat_diffuse[i];
      mat_struct.emission[i] = def_mat_emission[i];
    }
    first_time = 0;
  }

  glMaterialfv( GL_FRONT, GL_SHININESS, mat_struct.shininess );
  glMaterialfv( GL_FRONT, GL_SPECULAR,  mat_struct.specular ); 
  glMaterialfv( GL_FRONT, GL_AMBIENT,   mat_struct.ambient );
  glMaterialfv( GL_FRONT, GL_DIFFUSE,   mat_struct.diffuse );
  glMaterialfv( GL_FRONT, GL_EMISSION,  mat_struct.emission );

}



/*===========================================================================*/
void
LoadIsoMaterial(int type)
{
  int i;
  static int first_time = 1;

  if ( ! togl_exists ) return;

  if ( first_time ) {
    /* ISOSURF_TRANSP_OFF = 0 */
    front_mat_isosurf[0].shininess[0] = def_front_mat_shininess[0];
    back_mat_isosurf[0].shininess[0]  = def_back_mat_shininess[0];
    
    front_mat_pos_isosurf[0].shininess[0] = def_front_mat_shininess[0];
    back_mat_pos_isosurf[0].shininess[0]  = def_back_mat_shininess[0];

    front_mat_neg_isosurf[0].shininess[0] = def_front_mat_shininess[0];
    back_mat_neg_isosurf[0].shininess[0]  = def_back_mat_shininess[0];

    /* ISOSURF_TRANSP_ON  = 1 */    
    front_mat_isosurf[1].shininess[0] = def_blend_front_mat_shininess[0];
    back_mat_isosurf[1].shininess[0]  = def_blend_back_mat_shininess[0];
    
    front_mat_pos_isosurf[1].shininess[0] = def_blend_front_mat_shininess[0];
    back_mat_pos_isosurf[1].shininess[0]  = def_blend_back_mat_shininess[0];

    front_mat_neg_isosurf[1].shininess[0] = def_blend_front_mat_shininess[0];
    back_mat_neg_isosurf[1].shininess[0]  = def_blend_back_mat_shininess[0];

    /* COLORPLANE - TRANSP_OFF */
    mat_colorplane[0].shininess[0] = def_cpl_shininess[0];
    /*            - TRANSP_ON */
    mat_colorplane[1].shininess[0] = def_blend_cpl_shininess[0];

    for (i=0; i<4; i++) {
      /* ISOSURF_TRANSP_OFF = 0 */
      front_mat_isosurf[0].specular[i] = def_front_mat_specular[i];
      front_mat_isosurf[0].diffuse[i]  = def_front_mat_diffuse[i];
      front_mat_isosurf[0].ambient[i]  = def_front_mat_ambient[i];
      front_mat_isosurf[0].emission[i] = def_front_mat_emission[i];

      back_mat_isosurf[0].specular[i]  = def_back_mat_specular[i];
      back_mat_isosurf[0].diffuse[i]   = def_back_mat_diffuse[i];
      back_mat_isosurf[0].ambient[i]   = def_back_mat_ambient[i];
      back_mat_isosurf[0].emission[i]  = def_back_mat_emission[i];

      
      front_mat_pos_isosurf[0].specular[i] = def_front_mat_specular[i];
      front_mat_pos_isosurf[0].diffuse[i]  = def_pos_front_mat_diffuse[i];
      front_mat_pos_isosurf[0].ambient[i]  = def_pos_front_mat_ambient[i];
      front_mat_pos_isosurf[0].emission[i] = def_pos_front_mat_emission[i];

      back_mat_pos_isosurf[0].specular[i]  = def_back_mat_specular[i];
      back_mat_pos_isosurf[0].diffuse[i]   = def_pos_back_mat_diffuse[i];
      back_mat_pos_isosurf[0].ambient[i]   = def_pos_back_mat_ambient[i];
      back_mat_pos_isosurf[0].emission[i]  = def_pos_back_mat_emission[i];
      

      front_mat_neg_isosurf[0].specular[i] = def_front_mat_specular[i];
      front_mat_neg_isosurf[0].diffuse[i]  = def_neg_front_mat_diffuse[i];
      front_mat_neg_isosurf[0].ambient[i]  = def_neg_front_mat_ambient[i];
      front_mat_neg_isosurf[0].emission[i] = def_neg_front_mat_emission[i];

      back_mat_neg_isosurf[0].specular[i]  = def_back_mat_specular[i];
      back_mat_neg_isosurf[0].diffuse[i]   = def_neg_back_mat_diffuse[i];
      back_mat_neg_isosurf[0].ambient[i]   = def_neg_back_mat_ambient[i];
      back_mat_neg_isosurf[0].emission[i]  = def_neg_back_mat_emission[i];

      /* ISOSURF_TRANSP_ON = 1 */
      front_mat_isosurf[1].specular[i] = def_blend_front_mat_specular[i];
      front_mat_isosurf[1].diffuse[i]  = def_blend_front_mat_diffuse[i];
      front_mat_isosurf[1].ambient[i]  = def_blend_front_mat_ambient[i];
      front_mat_isosurf[1].emission[i] = def_blend_front_mat_emission[i];

      back_mat_isosurf[1].specular[i]  = def_blend_back_mat_specular[i];
      back_mat_isosurf[1].diffuse[i]   = def_blend_back_mat_diffuse[i];
      back_mat_isosurf[1].ambient[i]   = def_blend_back_mat_ambient[i];
      back_mat_isosurf[1].emission[i]  = def_blend_back_mat_emission[i];
      

      front_mat_pos_isosurf[1].specular[i]= def_blend_front_mat_specular[i];
      front_mat_pos_isosurf[1].diffuse[i] = def_blend_pos_front_mat_diffuse[i];
      front_mat_pos_isosurf[1].ambient[i] = def_blend_pos_front_mat_ambient[i];
      front_mat_pos_isosurf[1].emission[i] = def_blend_pos_front_mat_emission[i];

      back_mat_pos_isosurf[1].specular[i] = def_blend_back_mat_specular[i];
      back_mat_pos_isosurf[1].diffuse[i]  = def_blend_pos_back_mat_diffuse[i];
      back_mat_pos_isosurf[1].ambient[i]  = def_blend_pos_back_mat_ambient[i];
      back_mat_pos_isosurf[1].emission[i]  = def_blend_pos_back_mat_emission[i];
      

      front_mat_neg_isosurf[1].specular[i]= def_blend_front_mat_specular[i];
      front_mat_neg_isosurf[1].diffuse[i] = def_blend_neg_front_mat_diffuse[i];
      front_mat_neg_isosurf[1].ambient[i] = def_blend_neg_front_mat_ambient[i];
      front_mat_neg_isosurf[1].emission[i] = def_blend_neg_front_mat_emission[i];

      back_mat_neg_isosurf[1].specular[i] = def_blend_back_mat_specular[i];
      back_mat_neg_isosurf[1].diffuse[i]  = def_blend_neg_back_mat_diffuse[i];
      back_mat_neg_isosurf[1].ambient[i]  = def_blend_neg_back_mat_ambient[i];
      back_mat_neg_isosurf[1].emission[i]  = def_blend_neg_back_mat_emission[i];


      /* COLORPLANE - TRANSP OFF */
      mat_colorplane[0].specular[i] = def_cpl_specular[i];
      mat_colorplane[0].ambient[i]  = def_cpl_ambient[i];
      mat_colorplane[0].diffuse[i]  = def_cpl_diffuse[i];
      /*            - TRANSP ON */
      mat_colorplane[1].specular[i] = def_blend_cpl_specular[i];
      mat_colorplane[1].ambient[i]  = def_blend_cpl_ambient[i];
      mat_colorplane[1].diffuse[i]  = def_blend_cpl_diffuse[i];
    }
    first_time = 0;
  } 
  
  if ( type == (MAT_ONELEVEL | MAT_ISOSURF) ) {
    GLenum front = GL_FRONT, back  = GL_BACK;
    if ( frontface_isosurf[GLPAR_FACE_WHAT_ISOSURF_POS] == GL_CW ) {
      back  = GL_FRONT;
      front = GL_BACK;
    }
    glMaterialfv(front, GL_SHININESS, 
		 front_mat_isosurf[isoDisp.transparent].shininess);
    glMaterialfv(front, GL_SPECULAR,
		 front_mat_isosurf[isoDisp.transparent].specular);
    glMaterialfv(front, GL_DIFFUSE,
		 front_mat_isosurf[isoDisp.transparent].diffuse);
    glMaterialfv(front, GL_AMBIENT,
		 front_mat_isosurf[isoDisp.transparent].ambient);
    glMaterialfv(front, GL_EMISSION,
		 front_mat_isosurf[isoDisp.transparent].emission);

    glMaterialfv(back, GL_SHININESS,
		 back_mat_isosurf[isoDisp.transparent].shininess);
    glMaterialfv(back, GL_SPECULAR,
		 back_mat_isosurf[isoDisp.transparent].specular);
    glMaterialfv(back, GL_DIFFUSE,
		 back_mat_isosurf[isoDisp.transparent].diffuse);
    glMaterialfv(back, GL_AMBIENT,
		 back_mat_isosurf[isoDisp.transparent].ambient);
    glMaterialfv(back, GL_EMISSION,
		 back_mat_isosurf[isoDisp.transparent].emission);
  }

  if ( type == (MAT_POSLEVEL | MAT_ISOSURF) ) {
    GLenum front = GL_FRONT, back  = GL_BACK;
    if ( frontface_isosurf[GLPAR_FACE_WHAT_ISOSURF_POS] == GL_CW ) {
      back  = GL_FRONT;
      front = GL_BACK;
    }
    glMaterialfv(front, GL_SHININESS, 		   
		 front_mat_pos_isosurf[isoDisp.transparent].shininess);
    glMaterialfv(front, GL_SPECULAR,
		 front_mat_pos_isosurf[isoDisp.transparent].specular);
    glMaterialfv(front, GL_DIFFUSE,
		 front_mat_pos_isosurf[isoDisp.transparent].diffuse);
    glMaterialfv(front, GL_AMBIENT,
		 front_mat_pos_isosurf[isoDisp.transparent].ambient);
    glMaterialfv(front, GL_EMISSION,
		 front_mat_pos_isosurf[isoDisp.transparent].emission);

    glMaterialfv(back, GL_SHININESS,
		 back_mat_pos_isosurf[isoDisp.transparent].shininess);
    glMaterialfv(back, GL_SPECULAR,
		 back_mat_pos_isosurf[isoDisp.transparent].specular);
    glMaterialfv(back, GL_DIFFUSE,
		 back_mat_pos_isosurf[isoDisp.transparent].diffuse);
    glMaterialfv(back, GL_AMBIENT,
		 back_mat_pos_isosurf[isoDisp.transparent].ambient);
    glMaterialfv(back, GL_EMISSION,
		 back_mat_pos_isosurf[isoDisp.transparent].emission);
  }
  
  if ( type == (MAT_NEGLEVEL | MAT_ISOSURF)
       || type == (MAT_NEGLEVEL | MAT_ISOSURF | MAT_SOLID) ) {
    GLenum front = GL_FRONT, back  = GL_BACK;
    GLint twoside_lighting;
    glGetIntegerv (GL_LIGHT_MODEL_TWO_SIDE, &twoside_lighting);

    if ( twoside_lighting ) {
      if ( type & MAT_SOLID ) {
        if ( frontface_isosurf[GLPAR_FACE_WHAT_ISOSURF_NEG] == GL_CCW ) {
          back  = GL_FRONT;
          front = GL_BACK;
        }
      }
    } else {
      if ( frontface_isosurf[GLPAR_FACE_WHAT_ISOSURF_NEG] == GL_CW ) {
	back  = GL_FRONT;
	front = GL_BACK;
      }
    }
    
    glMaterialfv(front, GL_SHININESS, 		   
		 front_mat_neg_isosurf[isoDisp.transparent].shininess);
    glMaterialfv(front, GL_SPECULAR,
		 front_mat_neg_isosurf[isoDisp.transparent].specular);
    glMaterialfv(front, GL_DIFFUSE,
		 front_mat_neg_isosurf[isoDisp.transparent].diffuse);
    glMaterialfv(front, GL_AMBIENT,
		 front_mat_neg_isosurf[isoDisp.transparent].ambient);
    glMaterialfv(front, GL_EMISSION,
		 front_mat_neg_isosurf[isoDisp.transparent].emission);

    glMaterialfv(back, GL_SHININESS,
		 back_mat_neg_isosurf[isoDisp.transparent].shininess);
    glMaterialfv(back, GL_SPECULAR,
		 back_mat_neg_isosurf[isoDisp.transparent].specular);
    glMaterialfv(back, GL_DIFFUSE,
		 back_mat_neg_isosurf[isoDisp.transparent].diffuse);
    glMaterialfv(back, GL_AMBIENT,
		 back_mat_neg_isosurf[isoDisp.transparent].ambient);
    glMaterialfv(back, GL_EMISSION,
		 back_mat_neg_isosurf[isoDisp.transparent].emission);
  }

  if ( type == MAT_COLORPLANE ) {
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 		   
		 mat_colorplane[isoDisp.transparent].shininess);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,
		 mat_colorplane[isoDisp.transparent].specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,
		 mat_colorplane[isoDisp.transparent].diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,
		 mat_colorplane[isoDisp.transparent].ambient);
  }
}


/*===========================================================================*/
/* WARNING: this routine was written very fast and should be rewriten        */
void
LoadVoronoiMaterial(int voronoi_type)
{
  int i;
  
  if ( ! togl_exists ) return;

  /* for voronoi polyhedra we don't need GL_BACK -> use just GL_FRONT and
   * perform polygon cooling
   */

  if (voronoi_type == VORONOI_WIGNERSEITZ) {
    mat_voronoi[0].shininess[0] = def_wigner_shininess[0];
    for (i=0; i<4; i++) {
      mat_voronoi[0].specular[i] = def_wigner_specular[i];
      mat_voronoi[0].ambient[i]  = def_wigner_ambient[i];  
      mat_voronoi[0].diffuse[i]  = def_wigner_diffuse[i]; 
    }
  } else if (voronoi_type == VORONOI_BZ) {
    mat_voronoi[0].shininess[0] = def_bz_shininess[0];
    for (i=0; i<4; i++) {
      mat_voronoi[0].specular[i] = def_bz_specular[i];
      mat_voronoi[0].ambient[i]  = def_bz_ambient[i];  
      mat_voronoi[0].diffuse[i]  = def_bz_diffuse[i]; 
    }
  } else {
    mat_voronoi[0].shininess[0] = def_voronoi_shininess[0];
    for (i=0; i<4; i++) {
      mat_voronoi[0].specular[i] = def_voronoi_specular[i];
      mat_voronoi[0].ambient[i]  = def_voronoi_ambient[i];  
      mat_voronoi[0].diffuse[i]  = def_voronoi_diffuse[i]; 
    }
  }

  glMaterialfv(GL_FRONT, GL_SHININESS, 		   
	       mat_voronoi[0].shininess);
  glMaterialfv(GL_FRONT, GL_SPECULAR,
	       mat_voronoi[0].specular);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,
	       mat_voronoi[0].diffuse);
  glMaterialfv(GL_FRONT, GL_AMBIENT,
	       mat_voronoi[0].ambient);
}


/*===========================================================================*/
/* WARNING: this routine was written very fast and should be rewriten        */
void
LoadCageOrVecMaterial(int type)
{
  int i;
  static int first_time = 1;
  
  if ( ! togl_exists ) return;

  /* for voronoi polyhedra we don't need GL_BACK -> use just GL_FRONT and
   * perform polygon cooling
   */
  if (first_time) {
    mat_vec_and_cage[GLPAR_PRIMCAGE].shininess[0] = def_primcage_shininess[0];
    mat_vec_and_cage[GLPAR_CONVCAGE].shininess[0] = def_convcage_shininess[0];
    mat_vec_and_cage[GLPAR_SCCAGE].shininess[0]   = def_sccage_shininess[0];  
    mat_vec_and_cage[GLPAR_PRIMVEC].shininess[0]  = def_primvec_shininess[0]; 
    mat_vec_and_cage[GLPAR_CONVVEC].shininess[0]  = def_convvec_shininess[0]; 
    mat_vec_and_cage[GLPAR_SCVEC].shininess[0]    = def_scvec_shininess[0];
    mat_vec_and_cage[GLPAR_ISOPARALLE].shininess[0] = 
      def_isoparalle_shininess[0];

    for (i=0; i<4; i++) {
      mat_vec_and_cage[GLPAR_PRIMCAGE].specular[i] = def_primcage_specular[i]; 
      mat_vec_and_cage[GLPAR_PRIMCAGE].ambient[i]  = def_primcage_ambient[i];  
      mat_vec_and_cage[GLPAR_PRIMCAGE].diffuse[i]  = def_primcage_diffuse[i];  

      mat_vec_and_cage[GLPAR_CONVCAGE].specular[i] = def_convcage_specular[i]; 
      mat_vec_and_cage[GLPAR_CONVCAGE].ambient[i]  = def_convcage_ambient[i];  
      mat_vec_and_cage[GLPAR_CONVCAGE].diffuse[i]  = def_convcage_diffuse[i];  

      mat_vec_and_cage[GLPAR_SCCAGE].specular[i] = def_sccage_specular[i];   
      mat_vec_and_cage[GLPAR_SCCAGE].ambient[i]  = def_sccage_ambient[i];    
      mat_vec_and_cage[GLPAR_SCCAGE].diffuse[i]  = def_sccage_diffuse[i];    

      mat_vec_and_cage[GLPAR_PRIMVEC].specular[i] = def_primvec_specular[i];  
      mat_vec_and_cage[GLPAR_PRIMVEC].ambient[i]  = def_primvec_ambient[i];   
      mat_vec_and_cage[GLPAR_PRIMVEC].diffuse[i]  = def_primvec_diffuse[i];   

      mat_vec_and_cage[GLPAR_CONVVEC].specular[i] = def_convvec_specular[i];  
      mat_vec_and_cage[GLPAR_CONVVEC].ambient[i]  = def_convvec_ambient[i];   
      mat_vec_and_cage[GLPAR_CONVVEC].diffuse[i]  = def_convvec_diffuse[i];   

      mat_vec_and_cage[GLPAR_SCVEC].specular[i] = def_scvec_specular[i];
      mat_vec_and_cage[GLPAR_SCVEC].ambient[i]  = def_scvec_ambient[i];
      mat_vec_and_cage[GLPAR_SCVEC].diffuse[i]  = def_scvec_diffuse[i];

      mat_vec_and_cage[GLPAR_ISOPARALLE].specular[i] = 
	def_isoparalle_specular[i];
      mat_vec_and_cage[GLPAR_ISOPARALLE].ambient[i]  = 
	def_isoparalle_ambient[i];
      mat_vec_and_cage[GLPAR_ISOPARALLE].diffuse[i]  = 
	def_isoparalle_diffuse[i];
    }
    first_time = 0;
  }

  if ( type != GLPAR_ISOPARALLE && xcr.dim == 3 ) {
    glMaterialfv(GL_FRONT, GL_SHININESS, 		   
		 mat_vec_and_cage[type].shininess);
    glMaterialfv(GL_FRONT, GL_SPECULAR,
		 mat_vec_and_cage[type].specular);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,
		 mat_vec_and_cage[type].diffuse);
    glMaterialfv(GL_FRONT, GL_AMBIENT,
		 mat_vec_and_cage[type].ambient);
  } else {
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 		   
		 mat_vec_and_cage[type].shininess);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,
		 mat_vec_and_cage[type].specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,
		 mat_vec_and_cage[type].diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,
		 mat_vec_and_cage[type].ambient);
  }

}
