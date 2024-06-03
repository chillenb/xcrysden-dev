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
 * Source: $XCRYSDEN_TOPDIR/C/setOpenGLState.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include "xcGLparam.h"

void SetCageOGLState( int dim, BLENDFUNC blend, GLint *shade_model, 
		      GLboolean *two_side, GLboolean *cull_face );

void DisableCageOGLState( int dim, GLint shade_model, 
			  GLboolean two_side, GLboolean cull_face );

/* SetDefaultOGLState will have these:
   ----------------------------------
   GL_BLEND
   GL_SHADE_MODEL
   GL_LIGHT_MODEL_TWO_SIDE
   GL_CULL_FACE

*/

void
SetCageOGLState( int dim, BLENDFUNC blend, GLint *shade_model, 
		 GLboolean *two_side, GLboolean *cull_face )
{
  /* now get the right state for rendering the cellcage */

  /* enable transparency */
  glEnable( GL_BLEND );
  glDepthMask( GL_FALSE );
  glBlendFunc( blend.sfunc, blend.dfunc );
  
  /* SHADEMODEL */
  glGetIntegerv( GL_SHADE_MODEL, shade_model );
  glShadeModel( GL_FLAT );

  glGetBooleanv( GL_LIGHT_MODEL_TWO_SIDE, two_side );
  glGetBooleanv( GL_CULL_FACE, cull_face );
  if ( dim == 3 ) {
    if ( *two_side == GL_TRUE )
      /* disable GL_LIGHT_MODEL_TWO_SIDE */
      glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE );
    
    /* disable CULL FACE */
    if ( (*cull_face) ) 
      glDisable( GL_CULL_FACE ); 
  } else if ( dim == 2 ) {
    if ( *two_side == GL_FALSE )
      /* enable GL_LIGHT_MODEL_TWO_SIDE */
      glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );
    
    /* disable CULL FACE */
    if ( (*cull_face) ) 
      glDisable( GL_CULL_FACE ); 
  }
}


/* 
   in future this function will no longer exists -> when structure will 
   have its own OGL State function, i.e. when DEFAULT OGL State function
   will exists
*/
void
DisableCageOGLState( int dim, GLint shade_model, 
		     GLboolean two_side, GLboolean cull_face )
{
  GLint get_shade_model;
  GLboolean get_two_side, get_cull_face;

  /* DISABLE TRANSPARENCY */  
  glDisable(GL_BLEND); 
  glDepthMask(GL_TRUE);  
  
  /* go back to old shademodel */
  glGetIntegerv( GL_SHADE_MODEL, &get_shade_model );
  if ( shade_model != get_shade_model )
    glShadeModel( shade_model );
    
  /* go to old GL_LIGHT_MODEL_TWO_SIDE */
  /* obstaja funkcija LoadLightModelTwoSide(int type) */
  glGetBooleanv( GL_LIGHT_MODEL_TWO_SIDE, &get_two_side );
  if ( two_side != get_two_side )
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, two_side );

  glGetBooleanv( GL_CULL_FACE, &get_cull_face );  
  if ( cull_face != get_cull_face ) {
    if ( cull_face == GL_TRUE )
      glEnable( GL_CULL_FACE );
    else if ( cull_face == GL_FALSE )
      glDisable( GL_CULL_FACE );
  }
}

