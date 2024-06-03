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
 * Source: $XCRYSDEN_TOPDIR/C/lighting.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

/*****************************************************************************/
/* this is for structure's material                                          */
extern GLfloat def_mat_shininess[];
extern GLfloat def_mat_specular[];
extern GLfloat def_mat_ambient[];
extern GLfloat def_mat_diffuse[];
extern GLfloat def_mat_emission[];

/*****************************************************************************/
/* this is for colorplanes's material                                        */
extern GLfloat def_cpl_shininess[];
extern GLfloat def_cpl_specular[];
extern GLfloat def_cpl_ambient[];
extern GLfloat def_cpl_diffuse[];

extern GLfloat def_blend_cpl_shininess[];
extern GLfloat def_blend_cpl_specular[];
extern GLfloat def_blend_cpl_ambient[];
extern GLfloat def_blend_cpl_diffuse[];

/***************************************************************************/
/* for ISOSURF_TRANSP_ON                                                   */
extern float def_blend_front_mat_shininess[];
extern float def_blend_front_mat_specular[];
/* s tem kontroliram prednjo barvo, manjse vrednoti ->  */
extern float def_blend_front_mat_diffuse[];
extern float def_blend_front_mat_ambient[];
extern float def_blend_front_mat_emission[];

extern float def_blend_pos_front_mat_diffuse[];
extern float def_blend_pos_front_mat_ambient[];
extern float def_blend_pos_front_mat_emission[];
extern float def_blend_neg_front_mat_diffuse[];
extern float def_blend_neg_front_mat_ambient[];
extern float def_blend_neg_front_mat_emission[];

/* s tem kotroliram koliko se vidi steklo, ce za njim ni nobenega */
extern float def_blend_back_mat_shininess[];
extern float def_blend_back_mat_specular[];
/* s tem kontroliram, koliko se vidi zadnjo stran stekla */
extern float def_blend_back_mat_diffuse[];
extern float def_blend_back_mat_ambient[];
extern float def_blend_back_mat_emission[];

extern float def_blend_pos_back_mat_diffuse[];
extern float def_blend_pos_back_mat_ambient[];
extern float def_blend_pos_back_mat_emission[];
extern float def_blend_neg_back_mat_diffuse[];
extern float def_blend_neg_back_mat_ambient[];
extern float def_blend_neg_back_mat_emission[];

/***************************************************************************/
/* for ISOSURF_TRANSP_OFF                                                  */
extern float def_front_mat_shininess[];
extern float def_front_mat_specular[];
/* s tem kontroliram prednjo barvo, manjse vrednoti ->  */
extern float def_front_mat_diffuse[];
extern float def_front_mat_ambient[];
extern float def_front_mat_emission[];

extern float def_pos_front_mat_diffuse[];
extern float def_pos_front_mat_ambient[];
extern float def_pos_front_mat_emission[];
extern float def_neg_front_mat_diffuse[];
extern float def_neg_front_mat_ambient[];
extern float def_neg_front_mat_emission[];

/* s tem kotroliram koliko se vidi steklo, ce za njim ni nobenega */
extern float def_back_mat_shininess[];
extern float def_back_mat_specular[];
/* s tem kontroliram, koliko se vidi zadnjo stran stekla */
extern float def_back_mat_diffuse[];
extern float def_back_mat_ambient[];
extern float def_back_mat_emission[];

extern float def_pos_back_mat_diffuse[];
extern float def_pos_back_mat_ambient[];
extern float def_pos_back_mat_emission[];
extern float def_neg_back_mat_diffuse[];
extern float def_neg_back_mat_ambient[];
extern float def_neg_back_mat_emission[];

/*****************************************************************************/
/* LIGHTS                                                                    */
#define MNOL MAX_NUMBER_OF_LIGHTS
extern boolean def_light_isenabled[MNOL];
extern GLfloat def_light_ambient[MNOL][4];
extern GLfloat def_light_diffuse[MNOL][4];
extern GLfloat def_light_specular[MNOL][4];
extern GLfloat def_light_fract_position[MNOL][4];
extern GLfloat def_light_spotdir[MNOL][3];
extern GLfloat def_light_spot_exp[MNOL][1];
extern GLfloat def_light_spot_cutoff[MNOL][1];
extern GLfloat def_light_const_atten[MNOL][1];
extern GLfloat def_light_lin_atten[MNOL][1];
extern GLfloat def_light_quad_atten[MNOL][1];

/*****************************************************************************/
/*   LIGHTMODEL --- LIGHTMODEL --- LIGHTMODEL --- LIGHTMODEL --- LIGHTMODEL  */
extern GLfloat def_lightmodel_two_side[1];
extern GLfloat def_lightmodel_two_side_iso[1];
extern GLfloat def_lightmodel_ambient[4];
extern GLfloat def_lightmodel_local_viewer[1];

/*****************************************************************************/
/*   FRONTFACE & BLENDFUNC   ---   FRONTFACE & BLENDFUNC                     */
extern FRONTFACE def_frontface_isosurf;
extern BLENDFUNC def_blend_isosurf;
extern BLENDFUNC def_blend_voronoi;


