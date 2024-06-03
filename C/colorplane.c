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
 * Source: $XCRYSDEN_TOPDIR/C/colorplane.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <stdlib.h>
#include <math.h>
#include <GL/gl.h>
#include <tk.h>
#include "struct.h"
#include "isosurf.h"
#include "xcfunc.h"

extern PLANEVERTEX ***plvertex;

typedef struct {
  float plvec[2][3];
  float plnml[3];
} PLANEVEC;


PLANEVEC   pln;
static BASECOLOR  bcol;

float monocol[2][4] = {
  { 0.1, 0.1, 0.1, 0.95 },
  { 1.0, 1.0, 1.0, 0.95 }};

float rainbow[6][4] = {
  { 1.0, 0.0, 0.0, 0.95 },
  { 1.0, 1.0, 0.0, 0.95 },
  { 0.0, 1.0, 0.0, 0.95 },
  { 0.0, 1.0, 1.0, 0.95 },
  { 0.0, 0.0, 1.0, 0.95 },
  { 1.0, 0.0, 1.0, 0.95 }};

float rgbcol[3][4] = {
  { 1.0, 0.0, 0.0, 0.95 },
  { 0.0, 1.0, 0.0, 0.95 },
  { 0.0, 0.0, 1.0, 0.95 }};

float geographic[11][4]= {
  {   0./255.,   0./255., 204./255., 0.95 },
  { 150./255., 153./255., 255./255., 0.95 },
  { 192./255., 216./255., 255./255., 0.95 },
  {   0./255., 192./255.,   3./255., 0.95 },
  { 103./255., 255./255.,  51./255., 0.95 },
  { 206./255., 255./255.,  51./255., 0.95 },
  { 255./255., 255./255., 124./255., 0.95 },
  { 223./255., 191./255., 143./255., 0.95 },
  { 157./255.,  48./255., 152./255., 0.95 },
  { 255./255.,   0./255., 255./255., 0.95 },
  { 255./255., 255./255., 255./255., 0.95 }};

float blue_white_red[5][4] = {
  {0.0, 0.0, 1.0, 0.95},
  {0.6, 0.6, 1.0, 0.95},
  {1.0, 1.0, 1.0, 0.95},
  {1.0, 0.6, 0.6, 0.95},
  {1.0, 0.0, 0.0, 0.95}
};

float black_brown_white[4][4] = {
  {0.00, 0.00, 0.00, 0.95},
  {0.48, 0.14, 0.02, 0.95},
  {0.78, 0.63, 0.07, 0.95},
  {1.00, 1.00, 1.00, 0.95}
};

/* 
 * if grid value is out of prescribed range, make it in outcolor 
 */
float outcolor[4] = { 0.0, 0.0, 0.0, 0.95 };

static unsigned short line_stipple[10] = {
  LINESTIPPLE0,
  LINESTIPPLE1,
  LINESTIPPLE2,
  LINESTIPPLE3,
  LINESTIPPLE4,
  LINESTIPPLE5,
  LINESTIPPLE6,
  LINESTIPPLE7,
  LINESTIPPLE8,
  LINESTIPPLE9};
static int line_stipple_factor[4] = {
  LINESTIPPLE_FACTOR0,
  LINESTIPPLE_FACTOR1,
  LINESTIPPLE_FACTOR2,
  LINESTIPPLE_FACTOR3};

/* --- function prototypes --- */
static void ScaleFunctionType( float (*Func)(float value));
static float (*ScaleFunc)(float value) = 0;
float xcLinf(float value);
float xcLogf(float value);
float xcLog10f(float value);
float xcSqrtf(float value);
float xcRoot3f(float value);
float xcGaussf(float value);
float xcSlaterf(float value);
float xcExpf(float value);
float xcExp2f(float value);
/*
  extern float logf(float value);
  extern float log10f(float value);
  extern float sqrtf(float value);
  extern float expf(float value);
*/

void
ColorPlane(int obj, int cb, int fn, int nx, int ny)
{
  int i,j,k;
  float min, max, fac;
  float shift = 0.0;

  min = isodata.min_allowed[obj];
  max = isodata.max_allowed[obj];

  if ( fn == SCALE_FUNC_LIN )    ScaleFunctionType(xcLinf);
  if ( fn == SCALE_FUNC_LOG )    ScaleFunctionType(xcLogf);
  if ( fn == SCALE_FUNC_LOG10 )  ScaleFunctionType(xcLog10f);
  if ( fn == SCALE_FUNC_SQRT )   ScaleFunctionType(xcSqrtf);
  if ( fn == SCALE_FUNC_ROOT3 )  ScaleFunctionType(xcRoot3f);
  if ( fn == SCALE_FUNC_GAUSS )  ScaleFunctionType(xcGaussf);
  if ( fn == SCALE_FUNC_SLATER ) ScaleFunctionType(xcSlaterf);

  if ( (ScaleFunc == xcLogf || ScaleFunc == xcLog10f ||
	ScaleFunc == xcSqrtf || ScaleFunc == xcRoot3f ||
	ScaleFunc == xcGaussf ) && min < 1.e-7 ) {
    /* shift = -(min + 1e-7); OLD */
    shift = -min + 1e-7;
    min   = 1.e-7;
    max = max + shift;
  } 
    
  printf("shift, min, max:: %f %f %f\n",shift,min,max);
  if ( cb == COLORBASE_MONO ) {
    /* monocolor */
    bcol.ncol = 2;
    bcol.baseintv[0] = (*ScaleFunc)(min);
    bcol.baseintv[1] = (*ScaleFunc)(max);
    for(i=0; i<2; i++)
      for(j=0; j<4; j++)
	bcol.basecol[i][j] = monocol[i][j];
  } else if ( cb == COLORBASE_RAINBOW ) {
    float intv;
    /* rainbow */
    bcol.ncol = 6;
    intv = ((*ScaleFunc)(max) - (*ScaleFunc)(min)) / 5.0;
    printf("intv = %f\n",intv);
    for (i=0;i<6;i++) {
      bcol.baseintv[i] = (*ScaleFunc)(min) + ((float) i * intv); 
      printf("bcol.baseintv[i] = %f\n",bcol.baseintv[i]);
    }
    for(i=0; i<6; i++)
      for(j=0; j<4; j++)
	bcol.basecol[i][j] = rainbow[i][j];
  } else if ( cb == COLORBASE_RGB ) {
    float intv;
    /* rgb */
    bcol.ncol = 3;
    intv = ((*ScaleFunc)(max) - (*ScaleFunc)(min)) / 2.0;
    printf("intv = %f\n",intv);
    for (i=0;i<3;i++) {
      bcol.baseintv[i] = (*ScaleFunc)(min) + ((float) i * intv); 
      printf("bcol.baseintv[i] = %f\n",bcol.baseintv[i]);
    }
    for(i=0; i<3; i++)
      for(j=0; j<4; j++)
	bcol.basecol[i][j] = rgbcol[i][j];
  } else if ( cb == COLORBASE_GEOGRAPHIC ) {
    float intv;
    /* rainbow */
    bcol.ncol = 11;
    intv = ((*ScaleFunc)(max) - (*ScaleFunc)(min)) / 10.0;
    printf("intv = %f\n",intv);
    for (i=0;i<11;i++) {
      bcol.baseintv[i] = (*ScaleFunc)(min) + ((float) i * intv); 
      printf("bcol.baseintv[i] = %f\n",bcol.baseintv[i]);
    }
    for(i=0; i<11; i++)
      for(j=0; j<4; j++)
	bcol.basecol[i][j] = geographic[i][j];
  } else if ( cb == COLORBASE_BLUE_WHITE_RED ) {
    float intv;
    /* blue-white-red */
    bcol.ncol = 5;
    intv = ((*ScaleFunc)(max) - (*ScaleFunc)(min)) / (bcol.ncol - 1);
    printf("intv = %f\n",intv);
    for (i=0;i<bcol.ncol;i++) {
      bcol.baseintv[i] = (*ScaleFunc)(min) + ((float) i * intv); 
      printf("bcol.baseintv[i] = %f\n",bcol.baseintv[i]);
    }
    for(i=0; i<bcol.ncol; i++)
      for(j=0; j<4; j++)
	bcol.basecol[i][j] = blue_white_red[i][j];
  } else if ( cb == COLORBASE_BLACK_BROWN_WHITE ) {
    float intv;
    /* black-brown-white */
    bcol.ncol = 4;
    intv = ((*ScaleFunc)(max) - (*ScaleFunc)(min)) / (bcol.ncol - 1);
    printf("intv = %f\n",intv);
    for (i=0;i<bcol.ncol;i++) {
      bcol.baseintv[i] = (*ScaleFunc)(min) + ((float) i * intv); 
      printf("bcol.baseintv[i] = %f\n",bcol.baseintv[i]);
    }
    for(i=0; i<bcol.ncol; i++)
      for(j=0; j<4; j++)
	bcol.basecol[i][j] = black_brown_white[i][j];
  }

  /* now assign color */
  if ( ScaleFunc == xcGaussf || ScaleFunc == xcSlaterf ) {
    /* printf("SLATER OR GAUSS"); */
    for (j=0; j<nx; j++) 
      for (k=0; k<ny; k++) {
	if ( plvertex[obj][j][k].val > isodata.min_allowed[obj] &&
	     plvertex[obj][j][k].val < isodata.max_allowed[obj] ) {
	  for (i=0; i<bcol.ncol-1; i++)
	    if ((*ScaleFunc)(plvertex[obj][j][k].val + shift) <= bcol.baseintv[i  ] &&
		(*ScaleFunc)(plvertex[obj][j][k].val + shift) >= bcol.baseintv[i+1] ) 
	      {

		fac = (bcol.baseintv[i] - (*ScaleFunc)(plvertex[obj][j][k].val + shift)) / (bcol.baseintv[i] - bcol.baseintv[i+1]);	    

		plvertex[obj][j][k].col[0] = (1. - fac) * bcol.basecol[i][0] +  fac * bcol.basecol[i+1][0];
		plvertex[obj][j][k].col[1] = (1. - fac) * bcol.basecol[i][1] +  fac * bcol.basecol[i+1][1];
		plvertex[obj][j][k].col[2] = (1. - fac) * bcol.basecol[i][2] +  fac * bcol.basecol[i+1][2];
		plvertex[obj][j][k].col[3] = (1. - fac) * bcol.basecol[i][3] +  fac * bcol.basecol[i+1][3];
		
		break;	    
	      }  
	}
	else if ( plvertex[obj][j][k].val < isodata.min_allowed[obj] ) {
	  plvertex[obj][j][k].col[0] = bcol.basecol[0][0];
	  plvertex[obj][j][k].col[1] = bcol.basecol[0][1]; 
	  plvertex[obj][j][k].col[2] = bcol.basecol[0][2]; 
	  plvertex[obj][j][k].col[3] = bcol.basecol[0][3];
	} 
	else {
	  plvertex[obj][j][k].col[0] = bcol.basecol[bcol.ncol-1][0];
	  plvertex[obj][j][k].col[1] = bcol.basecol[bcol.ncol-1][1];
	  plvertex[obj][j][k].col[2] = bcol.basecol[bcol.ncol-1][2];
	  plvertex[obj][j][k].col[3] = bcol.basecol[bcol.ncol-1][3];	  
	}
      }
  } else { 
    for (j=0; j<nx; j++) 
      for (k=0; k<ny; k++) {
	if ( plvertex[obj][j][k].val > isodata.min_allowed[obj] &&
	     plvertex[obj][j][k].val < isodata.max_allowed[obj] ) {
	  for (i=0; i<bcol.ncol-1; i++)
	    if ((*ScaleFunc)(plvertex[obj][j][k].val + shift) >= bcol.baseintv[i  ] &&
		(*ScaleFunc)(plvertex[obj][j][k].val + shift) <=  bcol.baseintv[i+1] )
	      {
		fac = ((*ScaleFunc)(plvertex[obj][j][k].val + shift) - bcol.baseintv[i]) / (bcol.baseintv[i+1] - bcol.baseintv[i]);	    

		plvertex[obj][j][k].col[0] = (1. - fac) * bcol.basecol[i][0] + fac * bcol.basecol[i+1][0];
		plvertex[obj][j][k].col[1] = (1. - fac) * bcol.basecol[i][1] + fac * bcol.basecol[i+1][1];
		plvertex[obj][j][k].col[2] = (1. - fac) * bcol.basecol[i][2] + fac * bcol.basecol[i+1][2];
		plvertex[obj][j][k].col[3] = (1. - fac) * bcol.basecol[i][3] + fac * bcol.basecol[i+1][3];
		break;
		/*
		  printf("val:: %f;   fac:: %f\n",plvertex[obj][j][k].val,fac);
		  printf("i, i+1, basecol:: %d %d; %f %f %f\n                 %f %f %f\n",i,i+1,
		  bcol.basecol[i][1],bcol.basecol[i][2],bcol.basecol[i][3],
		  bcol.basecol[i+1][0],bcol.basecol[i+1][1],bcol.basecol[i+1][2]);
		  printf(">> color: %f %f %f %f\n",
		  plvertex[obj][j][k].col[0],
		  plvertex[obj][j][k].col[1],
		  plvertex[obj][j][k].col[2],
		  plvertex[obj][j][k].col[3]);
		*/
	      }
	}	
    	else if ( plvertex[obj][j][k].val < isodata.min_allowed[obj] ) {
	  plvertex[obj][j][k].col[0] = bcol.basecol[0][0];
	  plvertex[obj][j][k].col[1] = bcol.basecol[0][1]; 
	  plvertex[obj][j][k].col[2] = bcol.basecol[0][2]; 
	  plvertex[obj][j][k].col[3] = bcol.basecol[0][3];
	} 
	else {
	  plvertex[obj][j][k].col[0] = bcol.basecol[bcol.ncol-1][0];
	  plvertex[obj][j][k].col[1] = bcol.basecol[bcol.ncol-1][1];
	  plvertex[obj][j][k].col[2] = bcol.basecol[bcol.ncol-1][2];
	  plvertex[obj][j][k].col[3] = bcol.basecol[bcol.ncol-1][3];	  
	}
      }
  }
}


/*****************************************************************************/
/* this is for IsoLine2D --- IsoLine2D --- IsoLine2D --- IsoLine2D           */
/*****************************************************************************/
void
GetIsoLine2D_Attributes(int obj, int cb, int fn)
{
  register int i, j, k, l;
  register int nlevel1;
  float min, max, d, f, fac;
  float x, shift = 0.0;

  nlevel1= isoline2D[obj].nlevel - 1;
  min    = isodata.min_allowed[obj];
  max    = isodata.max_allowed[obj];

  if ( fn == SCALE_FUNC_LIN )    ScaleFunctionType(xcLinf);
  if ( fn == SCALE_FUNC_LOG )    ScaleFunctionType(xcLogf);
  if ( fn == SCALE_FUNC_LOG10 )  ScaleFunctionType(xcLog10f);
  if ( fn == SCALE_FUNC_SQRT )   ScaleFunctionType(xcSqrtf);
  if ( fn == SCALE_FUNC_ROOT3 )  ScaleFunctionType(xcRoot3f);
  if ( fn == SCALE_FUNC_GAUSS )  ScaleFunctionType(xcGaussf);
  if ( fn == SCALE_FUNC_SLATER ) ScaleFunctionType(xcSlaterf);
  

  /*******************/
  /*  assign LEVELS  */
  /*******************/
  /* this is just in case */
  if (max < min) {
    d = max;
    max = min;
    min = d;
  }

  if ( (ScaleFunc == xcLogf || ScaleFunc == xcLog10f ||
	ScaleFunc == xcSqrtf || ScaleFunc == xcRoot3f ||
	ScaleFunc == xcGaussf) && min < 1.e-7) {
    shift = -(min - 1.e-7);
    min   = 1.e-7;
    max = max + shift;
  }

  d = ScaleFunc(max) - ScaleFunc(min);  
  /* get the levels, i.e. isoline2D.level[] */
  for (i=0; i<=nlevel1; i++) {
    x = (*ScaleFunc)(min) + (float) i / nlevel1 * d;
    if ( fn == SCALE_FUNC_LIN)
      isoline2D[obj].level[i] = x - shift;
    if ( fn == SCALE_FUNC_LOG )
      isoline2D[obj].level[i] = exp(x) - shift;
    if ( fn == SCALE_FUNC_LOG10 )
      isoline2D[obj].level[i] = pow(10.0, x) - shift;
    if ( fn == SCALE_FUNC_SQRT )
      isoline2D[obj].level[i] = pow( x, 2.0) - shift;
    if ( fn == SCALE_FUNC_ROOT3 )
      isoline2D[obj].level[i] = pow( x, 3.0) - shift;
    if ( fn == SCALE_FUNC_GAUSS ) {
      /* d is negative && x is always positive */
      if ( x < 1.e-8 ) x = 1.e-8;
      f = log(x);
      if ( f < 0 ) f *= -1.0;
      isoline2D[obj].level[i] = sqrt( f ) - shift;
    }
    if ( fn == SCALE_FUNC_SLATER ) {
      /* d is negative -> x is or positive or negative */
      if ( x < 1.e-8 && x > -1.e-8 ) x = 1.e-8;
      if ( x < 0.0 ) x *= -1.0;
      f = log(x);
      isoline2D[obj].level[i] = f - shift;
    }    
  }

  /***************************/
  /* assign DASH & MONOCOLOR */
  /***************************/
  j = 0; /* this is for FULLDASH */
  k = 0; /* this is for FULLDASH */
  for (i=0; i<isoline2D[obj].nlevel; i++) {
    /* dash */
    if ( isoline2D[obj].linedash == ISOLINE_NODASH ) {
      isoline2D[obj].dash[i]       = LINESTIPPLE_SOLID;
      isoline2D[obj].dashfactor[1] = 1;
    }
    else if ( isoline2D[obj].linedash == ISOLINE_NEGDASH ) {
      if ( isoline2D[obj].level[i] < 0 ) {
	isoline2D[obj].dash[i]       = LINESTIPPLE0;
	isoline2D[obj].dashfactor[1] = LINESTIPPLE_FACTOR0; 
      } else {
	isoline2D[obj].dash[i]       = LINESTIPPLE_SOLID;
	isoline2D[obj].dashfactor[1] = 1;
      }
    } 
    else if ( isoline2D[obj].linedash == ISOLINE_FULLDASH ) {
      /* insert the code */
      isoline2D[obj].dash[i]       = line_stipple[j];
      isoline2D[obj].dashfactor[i] = line_stipple_factor[k];
      k++;
      if ( k == 4 ) {
	k = 0;
	j++;
      }
    }
    /* monocolor */
    if ( isoline2D[obj].linecolor == ISOLINE_MONOCOLOR ) {
      for (l=0; l<4; l++)
	isoline2D[obj].color[i][l] = isoline2D[obj].monocolor[l];
    }
  }


  /********************/
  /* POLYCOLOR option */
  /********************/	
  if ( isoline2D[obj].linecolor == ISOLINE_POLYCOLOR ) {
    if ( cb == COLORBASE_MONO ) {
      /* monocolor base */
      bcol.ncol = 2;
      bcol.baseintv[0] = (*ScaleFunc)(min);
      bcol.baseintv[1] = (*ScaleFunc)(max);
      for(i=0; i<2; i++)
	for(j=0; j<4; j++)
	  bcol.basecol[i][j] = monocol[i][j];
    } else if ( cb == COLORBASE_RAINBOW ) {
      float intv;
      /* rainbow */
      bcol.ncol = 6;
      intv = ((*ScaleFunc)(max) - (*ScaleFunc)(min)) / 5.0;
      printf("intv = %f\n",intv);
      for (i=0;i<6;i++) {
	bcol.baseintv[i] = (*ScaleFunc)(min) + ((float) i * intv); 
	printf("bcol.baseintv[i] = %f\n",bcol.baseintv[i]);
      }
      for(i=0; i<6; i++)
	for(j=0; j<4; j++)
	  bcol.basecol[i][j] = rainbow[i][j];
    } else if ( cb == COLORBASE_RGB ) {
      float intv;
      /* rainbow */
      bcol.ncol = 3;
      intv = ((*ScaleFunc)(max) - (*ScaleFunc)(min)) / 2.0;
      printf("intv = %f\n",intv);
      for (i=0;i<3;i++) {
	bcol.baseintv[i] = (*ScaleFunc)(min) + ((float) i * intv); 
	printf("bcol.baseintv[i] = %f\n",bcol.baseintv[i]);
      }
      for(i=0; i<3; i++)
	for(j=0; j<4; j++)
	  bcol.basecol[i][j] = rgbcol[i][j];
    } else if ( cb == COLORBASE_GEOGRAPHIC ) {
      float intv;
      /* rainbow */
      bcol.ncol = 11;
      intv = ((*ScaleFunc)(max) - (*ScaleFunc)(min)) / 10.0;
      printf("intv = %f\n",intv);
      for (i=0;i<11;i++) {
	bcol.baseintv[i] = (*ScaleFunc)(min) + ((float) i * intv); 
	printf("bcol.baseintv[i] = %f\n",bcol.baseintv[i]);
      }
      for(i=0; i<11; i++)
	for(j=0; j<4; j++)
	  bcol.basecol[i][j] = geographic[i][j];
    } else if ( cb == COLORBASE_BLUE_WHITE_RED ) {
      float intv;
      /* blue-white-red */
      bcol.ncol = 5;
      intv = ((*ScaleFunc)(max) - (*ScaleFunc)(min)) / (bcol.ncol - 1);
      printf("intv = %f\n",intv);
      for (i=0;i<bcol.ncol;i++) {
	bcol.baseintv[i] = (*ScaleFunc)(min) + ((float) i * intv); 
	printf("bcol.baseintv[i] = %f\n",bcol.baseintv[i]);
      }
      for(i=0; i<bcol.ncol; i++)
	for(j=0; j<4; j++)
	  bcol.basecol[i][j] = blue_white_red[i][j];
    }  else if ( cb == COLORBASE_BLACK_BROWN_WHITE ) {
      float intv;
    /* black-brown-white */
    bcol.ncol = 4;
    intv = ((*ScaleFunc)(max) - (*ScaleFunc)(min)) / (bcol.ncol - 1);
    printf("intv = %f\n",intv);
    for (i=0;i<bcol.ncol;i++) {
      bcol.baseintv[i] = (*ScaleFunc)(min) + ((float) i * intv); 
      printf("bcol.baseintv[i] = %f\n",bcol.baseintv[i]);
    }
    for(i=0; i<bcol.ncol; i++)
      for(j=0; j<4; j++)
	bcol.basecol[i][j] = black_brown_white[i][j];
  }
    
    /* now assign color */
    i=0;
    if ( ScaleFunc == xcGaussf || ScaleFunc == xcSlaterf ) {
      /* printf("SLATER OR GAUSS"); */
      for (j=0; j<isoline2D[obj].nlevel; j++) 
	for (; i<bcol.ncol-1; i++)
	  if ((*ScaleFunc)(isoline2D[obj].level[j]+shift) <= 
	      bcol.baseintv[i] &&
	      (*ScaleFunc)(isoline2D[obj].level[j]+shift) >= 
	      bcol.baseintv[i+1])
	    {
	      fac = ( bcol.baseintv[i] - 
		      (*ScaleFunc)(isoline2D[obj].level[j] + shift) ) /
		(bcol.baseintv[i] - bcol.baseintv[i+1]);	    
	      isoline2D[obj].color[j][0] = (1. - fac) * bcol.basecol[i][0] + 
		fac * bcol.basecol[i+1][0];
	      isoline2D[obj].color[j][1] = (1. - fac) * bcol.basecol[i][1] + 
		fac * bcol.basecol[i+1][1];
	      isoline2D[obj].color[j][2] = (1. - fac) * bcol.basecol[i][2] + 
		fac * bcol.basecol[i+1][2];
	      isoline2D[obj].color[j][3] = (1. - fac) * bcol.basecol[i][3] + 
		fac * bcol.basecol[i+1][3];
	      break;	    
	    } 
    } else { 
      for (j=0; j<isoline2D[obj].nlevel; j++) 
	for (; i<bcol.ncol-1; i++)
	  if ((*ScaleFunc)(isoline2D[obj].level[j]+shift) >= 
	      bcol.baseintv[i] &&
	      (*ScaleFunc)(isoline2D[obj].level[j]+shift) <= 
	      bcol.baseintv[i+1] )
	    {
	      fac = ((*ScaleFunc)(isoline2D[obj].level[j] + shift) - 
		     bcol.baseintv[i])/
		(bcol.baseintv[i+1] - bcol.baseintv[i]);	    
	      isoline2D[obj].color[j][0] = (1. - fac) * bcol.basecol[i][0] + 
		fac * bcol.basecol[i+1][0];
	      isoline2D[obj].color[j][1] = (1. - fac) * bcol.basecol[i][1] + 
		fac * bcol.basecol[i+1][1];
	      isoline2D[obj].color[j][2] = (1. - fac) * bcol.basecol[i][2] + 
		fac * bcol.basecol[i+1][2];
	      isoline2D[obj].color[j][3] = (1. - fac) * bcol.basecol[i][3] + 
		fac * bcol.basecol[i+1][3];
	      break;
	    }
    }
  }
}


/*****************************************************************************/
static void
ScaleFunctionType( float (*Func)(float value))
{
  ScaleFunc = Func;
}

/* --- all this function accept only positive values --- */
float xcLinf(float value)
{
  return value;
}

float xcLogf(float value)
{
  return log(value);
}
 
float xcLog10f(float value)
{
  return log10(value);
}

float xcSqrtf(float value)
{
  return sqrt(value);
}

float xcRoot3f(float value)
{
  return exp((1./3.) * log(value));
}

float xcGaussf(float value)
{
  if (value > 10.0) return 0.0;
  else return exp(-(value * value));
}

float xcSlaterf(float value)
{
  if (value > 100.0) return 0.0;
  else return exp(-value);
}

float xcExpf(float value)
{
  return exp(value);
}


float xcExp2f(float value)
{
  return exp(value*value);
}


