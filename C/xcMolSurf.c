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
* Source: $XCRYSDEN_TOPDIR/C/xcMolSurf.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <togl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "xcGLparam.h"
#include "struct.h"
#include "isosurf.h"
#include "memory.h"
#include "molsurf.h"
#include "xcfunc.h"

#define DIST3(x,y,z)   ( sqrt((x)*(x) + (y)*(y) + (z)*(z)) )

extern StructSize ss;

MOL_SURF *mols=NULL, *molsPtr = NULL;

static void  MolSurfGauss_Init(MOL_SURF *mols);
static void  MolSurfExp_Init(MOL_SURF *mols);
static float MolSurfGauss(MOL_SURF *mols, int nat, float dis);
static float MolSurfExp(MOL_SURF *mols, int nat, float dis);
static float MolSurfGaussPlain(MOL_SURF *mols, int nat, float dis);
static float MolSurfExpPlain(MOL_SURF *mols, int nat, float dis);

extern void MarchCubeNew(ISOSURFACE *iso,
			 float ***gridvalue,
			 float vec[][3],
			 float origin[],
			 int ix, int iy, int iz,
			 float isolevel, int algorithm, int shade_model, int normals_model);
extern void xcSurfSmoothing( ISOSURFACE *iso, float isolevel );
extern void MolSurfGridAtomicColor(ISOSURFACE *iso, float mols_d);

MOL_SURF *FindMolSurf(int index)
{
  MOL_SURF *g = molsPtr;
  while (g) {
    if (index == g->index)
      return g;
    g = g->ptr;
  }
  return NULL;
}
static void AddToMolSurfList(MOL_SURF *g)
{
  g->ptr  = molsPtr;
  molsPtr = g;
}
int NewMolSurf(void) 
{
  static int index = -1;
  mols = (MOL_SURF *) calloc (1, sizeof(MOL_SURF) );
  mols->index=++index;
  mols->isosurf_index = -1;

  /* set to DEFAULTS */
  mols->type             = MOLS_GAUSSIAN;
  mols->rad              = rcov; /* alternative is rvdw */
  mols->colorscheme      = MOLS_ATOMIC;
  mols->dispt.drawstyle  = ISOSURF_SOLID;
  mols->dispt.shademodel = GL_SMOOTH;
  mols->dispt.transparent= ISOSURF_TRANSP_OFF;
  mols->d                = 0.4;
  /* allocate somewhat bigger vector to acount for 
     primcell -> convcell conversion */
  mols->atmindex         = (int *) malloc((size_t)(2*natoms) * sizeof(int));
  mols->smooth_nstep     = 0;
  mols->smooth_weight    = 0.5;

  mols->isoexpand.irepvec[0] = 1;
  mols->isoexpand.irepvec[1] = 1;
  mols->isoexpand.irepvec[2] = 1;

  mols->interp_molgrid = (float ***)NULL;
  mols->interp_n[0] = 1;
  mols->interp_n[1] = 1;
  mols->interp_n[2] = 1;

  mols->frontface     = GL_CW;
  mols->revertnormals = 0;

  mols->tessellation_algorithm = ISOSURF_MARCHING_CUBES;
  mols->normals_algorithm      = ISOSURF_NORMALS_GRADIENT;    

  AddToMolSurfList( mols );
  return index;
}


static float ln_of_2 = 0.69314718055994530941;
static void MolSurfGauss_Init(MOL_SURF *mols)
{
  int i;
  mols->c[0] = 0.0;
  for(i=1; i<=MAXNAT; i++)
    mols->c[i] = ln_of_2 / (mols->rad[i]*mols->rad[i]);
}
static void MolSurfExp_Init(MOL_SURF *mols)
{
  int i;
  mols->c[0] = 0.0;
  for(i=1; i<=MAXNAT; i++)
    mols->c[i] = ln_of_2 / mols->rad[i];
}
static float MolSurfGauss(MOL_SURF *mols, int nat, float dis)
{
  return 2.0 * exp(-mols->c[nat] * dis*dis);
}
static float MolSurfExp(MOL_SURF *mols, int nat, float dis)
{
  return 2.0 * exp(-mols->c[nat] * dis);
}
static float MolSurfGaussPlain(MOL_SURF *mols, int nat, float dis)
{
  return 2.0 * exp(-ln_of_2 * dis*dis);
}
static float MolSurfExpPlain(MOL_SURF *mols, int nat, float dis)
{
  return 2.0 * exp(-ln_of_2 * dis);
}


int 
XC_MolSurfRegCmd(ClientData clientData,Tcl_Interp *interp,
		 int argc, char *argv[])
{
  int ident;
  char *result = (char *) Tcl_Alloc((size_t) 10 * sizeof(char));

  /* <togl> so far ignored */
  ident=NewMolSurf();
  sprintf(result,"%d",ident);
  Tcl_SetResult(interp, result, TCL_DYNAMIC);

  VPf.surface    = (int *) realloc( VPf.surface, 
				    (size_t) (VPf.nsurface+1) * sizeof(int) );
  VPf.surfacePtr = (void **) realloc( VPf.surfacePtr, 
				      (size_t) (VPf.nsurface+1) * sizeof(void*) );

  VPf.surface[VPf.nsurface]      = 0;
  VPf.surfacePtr[VPf.nsurface++] = (void *) FindMolSurf( ident );
  return TCL_OK;
}

/*****************************************************************************
xc_molsurf <togl> -ident        <identifier>
-type         gauss|exp|unigauss|uniexp \
-radius       cov|VdW \
-level        <isoleve> \
-cutoff       <cutoff> \
-colorscheme  atomic|monochrome \
-drawstyle    solid|wire|dot \
-transparent  0|1 \
-shademodel   smooth|flat \
-monocolor    {<r> <g> <b>} \
-surfacetype  molsurf|gap
-resolution   <x> \
-smoothsteps  <steps> \
-smoothweight <weight> \
-tessellation cube|tetrahedron \
-normals      gradient|triangle
ADD:

xc_molsurfreg <togl>

xc_molsurfconfig <tolg> -ident <identifier> \
-render 0|1 \
-level <isolevel >
-colorscheme atomic|monochrome \
-drawstyle solid|wire|dot \
-transparent  0|1 \
-shademodel   smooth|flat
-monocolor    {<r> <g> <b>} \			  
-surfacetype  molsurf|gap \
-smoothsteps  <steps> \
-smoothweight <weight> \
-tessellation cube|tetrahedron \
-normals      gradient|triangle

write the display-function
complete the colorschemes
the distance colorcheme, can be made separate, which is better that grasping
from atomic-colorscheme !!!

ce imamo crystal, potem so drugi mols->lowcoor and mols->size; 
change that !!!

preveri se pointerje pri iso. Vzacetku jih postavi na NULL, 
potem pa v MarchCube, ce niso NULL jih najprej free();
*****************************************************************************/
int
XC_MolSurfConfigCmd(ClientData clientData,Tcl_Interp *interp,
		    int argc, char *argv[])
{
  register int i, j;
  int   s_nstep, s_algorithm, s_normals; 
  float s_weight, s_level;
  struct Togl *togl;
  ISOSURFACE *iso;

  /* <togl> so far ignored !!! */
  /* find togl associated with toglName */
  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* -IDENTIFIER */
  if ( strncmp(argv[2], "-iden", 5) == 0 ) {
    int ident;
    if ( Tcl_GetInt(interp, argv[3], &ident) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing %s %s %s %s ...", argv[3], argv[0], argv[1], argv[2], argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    mols = FindMolSurf( ident );
    if (!mols) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"couldn't find MolSurface %d !!!\n",
	       ident);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if (mols->isosurf_index < 0) {
      Tcl_SetResult(interp, "xc_molsurf must be called before xc_molsurfconfig!!!", TCL_STATIC);
      return TCL_ERROR;
    }
  } else {
    Tcl_SetResult(interp, "Usage: xc_molsurfconfig <tolg> -ident indent options ...", TCL_STATIC);
    return TCL_ERROR;
  }

  s_level  = mols->isolevel;
  s_nstep  = mols->smooth_nstep;
  s_weight = mols->smooth_weight;
  s_algorithm = mols->tessellation_algorithm;
  s_normals   = mols->normals_algorithm;

  for (i=4; i<argc; i+=2) {
    /* -RENDER */
    if ( strncmp(argv[i], "-rend", 5) == 0 ) {
      if ( strncmp(argv[i+1], "0", 1) == 0 ) 
	VPf.surface[mols->isosurf_index] = 0;
      else VPf.surface[mols->isosurf_index] = 1;
    }     
    /* -LEVEL */
    else if ( strncmp(argv[i], "-leve", 5) == 0 ) {
      double level;
      if ( Tcl_GetDouble(interp, argv[i+1], &level) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      mols->isolevel = (float) level;
    }
    /* -COLORSHEME */
    else if ( strncmp(argv[i], "-colo", 5) == 0 ) {
      if ( strncmp(argv[i+1], "atom", 4) == 0 )
	mols->colorscheme = MOLS_ATOMIC;
      else if ( strncmp(argv[i+1], "mono", 4) == 0 )
	mols->colorscheme = MOLS_MONOCHROME;	
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -colorscheme option\" command; option should be atomic or monocolor\n",
		 argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -DRAWSTYLE */
    else if ( strncmp(argv[i], "-draw", 5) == 0 ) {
      if ( strcmp(argv[i+1], "solid") == 0 ) 
	mols->dispt.drawstyle = ISOSURF_SOLID;
      else if ( strcmp(argv[i+1], "wire") == 0 ) 
	mols->dispt.drawstyle = ISOSURF_WIRE;
      else if ( strcmp(argv[i+1], "dot") == 0 ) 
	mols->dispt.drawstyle = ISOSURF_DOT;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -drawstyle option\" command; option should be solid, wire or dot\n",
		 argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -TRANSPARENT */
    else if ( strncmp(argv[i], "-tran", 5) == 0 ) {
      int transp;
      if ( Tcl_GetInt(interp, argv[i+1], &transp) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      mols->dispt.transparent = transp;
    }
    /* -SHADEMODEL */
    else if ( strncmp(argv[i], "-shade", 6) == 0 ) {
      if ( strncmp(argv[i+1], "smooth", 6) == 0 ) 
	mols->dispt.shademodel = GL_SMOOTH;
      else if ( strncmp(argv[i+1], "flat", 4) == 0 ) 
	mols->dispt.shademodel = GL_FLAT;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -shademodel option\" command; option should be smooth or flat\n",
		 argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -TESSELLATION */
    else if ( strncmp(argv[i], "-tessel", 6) == 0 ) {
      if ( strncmp(argv[i+1], "cube", 4) == 0 ) 
	mols->tessellation_algorithm = ISOSURF_MARCHING_CUBES;
      else if ( strncmp(argv[i+1], "tetra", 5) == 0 ) 
	mols->tessellation_algorithm = ISOSURF_TETRAHEDRAL;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -tessellation option\" command; option should be cube or tetrahedron\n",
		 argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -NORMALS */
    else if ( strncmp(argv[i], "-norm", 5) == 0 ) {
      if ( strncmp(argv[i+1], "grad", 4) == 0 ) 
	mols->normals_algorithm = ISOSURF_NORMALS_GRADIENT;
      else if ( strncmp(argv[i+1], "trian", 5) == 0 ) 
	mols->normals_algorithm = ISOSURF_NORMALS_TRIANGLE;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -tessellation option\" command; option should be cube or tetrahedron\n",
		 argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -MONOCOLOR */
    else if ( strncmp(argv[i], "-mono", 5) == 0 ) {
      double rgb;
      int argcList;
      const char **argvList;
      Tcl_SplitList(interp, argv[i+1], &argcList, &argvList);      
      for (j=0; j<argcList; j++) {
	if ( Tcl_GetDouble(interp, argvList[j], &rgb) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	mols->monocolor[j]      = (float) rgb;
	mols->back_monocolor[j] = (float) rgb-0.3;
	if ( mols->back_monocolor[j] < 0.0 )
	  mols->back_monocolor[j] = 0.0;
      }
      mols->monocolor[3]=0.5;
      mols->back_monocolor[3]=0.5;
    }
    /* -SMOOTHSTEPS */
    else if ( strncmp(argv[i], "-smooths", 8) == 0 ) {
      int nstep;
      if ( Tcl_GetInt(interp, argv[i+1], &nstep) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      mols->smooth_nstep = nstep;
    }
    /* -SMOOTHWEIGHT */
    else if ( strncmp(argv[i], "-smoothw", 8) == 0 ) {
      double weight;
      if ( Tcl_GetDouble(interp, argv[i+1], &weight) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      mols->smooth_weight = (float) weight;
    }
    else {
      Tcl_SetResult(interp, "Usage: xc_molsurf <tolg> -ident indent ...", TCL_STATIC);
      return TCL_ERROR;
    }
  }
  /***************************************************************************/

  iso = FindIsoSurf( mols->isosurf_index );

  iso->smooth_nstep  = mols->smooth_nstep;
  iso->smooth_weight = mols->smooth_weight;

  if ( s_level != mols->isolevel || s_algorithm != mols->tessellation_algorithm || s_normals != mols->normals_algorithm )
    {
      /*fprintf(stderr,"iso=%ld\n", iso);*/
      MarchCubeNew(iso, mols->molgrid, mols->vec, mols->lowcoor, 
		   mols->i, mols->j, mols->k, mols->isolevel, 
		   mols->tessellation_algorithm, mols->dispt.shademodel, mols->normals_algorithm);
    }

  if ( s_nstep  != iso->smooth_nstep  || s_weight != iso->smooth_weight || s_level != mols->isolevel ||
       s_algorithm != mols->tessellation_algorithm || s_normals != mols->normals_algorithm ) 
    {

      xcSurfSmoothing( iso, mols->isolevel );

      /* now add color */

      if ( mols->colorscheme == MOLS_ATOMIC ) {
	MolSurfGridAtomicColor( iso, mols->d );
      } 
      else if ( mols->colorscheme == MOLS_DISTANCE ) {
	;
      } 
      else if ( mols->colorscheme == MOLS_MONOCHROME ) {
	;
      }
    }

  s_nstep     = iso->smooth_nstep;
  s_weight    = iso->smooth_weight;
  s_level     = mols->isolevel;
  s_algorithm = mols->tessellation_algorithm;
  s_normals   = mols->normals_algorithm;

  UpdateProjection();
  Togl_PostRedisplay(togl);

  return TCL_OK;
}


int 
XC_MolSurfCmd(ClientData clientData,Tcl_Interp *interp,
	      int argc, char *argv[])
{
  register int i, j;
  struct Togl *togl;
  /*float ***molgrid;*/
  /*float Vec[3][3];*/
  ISOSURFACE *iso;

  /* <togl> so far ignored !!! */
  /* find togl associated with toglName */
  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* -IDENTIFIER */
  if ( strncmp(argv[2], "-iden", 5) == 0 ) {
    int ident;
    if ( Tcl_GetInt(interp, argv[3], &ident) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing %s %s %s %s ...", argv[3], argv[0], argv[1], argv[2], argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    /***************************/
    mols = FindMolSurf( ident );
    /***************************/
    if (!mols) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"couldn't find MolSurface %d !!!\n",
	       ident);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  } else {
    Tcl_SetResult(interp, "Usage: xc_molsurf <tolg> -ident indent ...", TCL_STATIC);
    return TCL_ERROR;
  }

  for (i=4; i<argc; i+=2) {
    /* -TYPE */
    if ( strcmp(argv[i], "-type") == 0 ) {
      if ( strncmp(argv[i+1], "gaus", 4) == 0 )  mols->type = MOLS_GAUSSIAN;
      else if ( strcmp(argv[i+1], "exp") == 0 )  mols->type = MOLS_EXP;
      else if ( strncmp(argv[i+1], "unig", 4) == 0 ) 
	mols->type = MOLS_UNIGAUSS;
      else if ( strncmp(argv[i+1], "unie", 4) == 0 ) 
	mols->type = MOLS_UNIEXP;
      else if ( strncmp(argv[i+1], "distf", 5) == 0 )
	mols->type = MOLS_DISTFUNC;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -type option\" command; option should be gauss, exp, unigauss, uniexp or distfunc\n", argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }     
    /* -RADIUS */
    else if ( strncmp(argv[i], "-radi", 5) == 0 ) {
      if ( strcmp(argv[i+1], "cov") == 0 )  mols->rad = rcov;
      else if ( strcmp(argv[i+1], "VdW") == 0 ) mols->rad = rvdw;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -radius option\" command; option should be cov or VdW\n", argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -LEVEL */
    else if ( strncmp(argv[i], "-leve", 5) == 0 ) {
      double level;
      if ( Tcl_GetDouble(interp, argv[i+1], &level) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      mols->isolevel = (float) level;
    }
    /* -CUTOFF */
    else if ( strncmp(argv[i], "-cuto", 5) == 0 ) {
      double cutoff;
      if ( Tcl_GetDouble(interp, argv[i+1], &cutoff) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      mols->cutoff = (float) cutoff;
    }
    /* -COLORSHEME */
    else if ( strncmp(argv[i], "-colo", 5) == 0 ) {
      if ( strncmp(argv[i+1], "atom", 4) == 0 )
	mols->colorscheme = MOLS_ATOMIC;
      else if ( strncmp(argv[i+1], "mono", 4) == 0 )
	mols->colorscheme = MOLS_MONOCHROME;	
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -colorscheme option\" command; option should be atomic or monocolor\n", argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -DRAWSTYLE */
    else if ( strncmp(argv[i], "-draw", 5) == 0 ) {
      if ( strcmp(argv[i+1], "solid") == 0 ) 
	mols->dispt.drawstyle = ISOSURF_SOLID;
      else if ( strcmp(argv[i+1], "wire") == 0 ) 
	mols->dispt.drawstyle = ISOSURF_WIRE;
      else if ( strcmp(argv[i+1], "dot") == 0 ) 
	mols->dispt.drawstyle = ISOSURF_DOT;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -drawstyle option\" command; option should be solid, wire or dot\n", argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -TRANSPARENT */
    else if ( strncmp(argv[i], "-tran", 5) == 0 ) {
      int transp;
      if ( Tcl_GetInt(interp, argv[i+1], &transp) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      mols->dispt.transparent = transp;
    }
    /* -SHADEMODEL */
    else if ( strncmp(argv[i], "-shade", 6) == 0 ) {
      if ( strncmp(argv[i+1], "smooth", 6) == 0 ) 
	mols->dispt.shademodel = GL_SMOOTH;
      else if ( strncmp(argv[i+1], "flat", 4) == 0 ) 
	mols->dispt.shademodel = GL_FLAT;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -shademodel option\" command; option should be smooth or flat\n",
		 argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -TESSELLATION */
    else if ( strncmp(argv[i], "-tessel", 6) == 0 ) {
      if ( strncmp(argv[i+1], "cube", 4) == 0 ) 
	mols->tessellation_algorithm = ISOSURF_MARCHING_CUBES;
      else if ( strncmp(argv[i+1], "tetra", 5) == 0 ) 
	mols->tessellation_algorithm = ISOSURF_TETRAHEDRAL;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -tessellation option\" command; option should be cube or tetrahedron\n",
		 argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -NORMALS */
    else if ( strncmp(argv[i], "-norm", 5) == 0 ) {
      if ( strncmp(argv[i+1], "grad", 4) == 0 ) 
	mols->normals_algorithm = ISOSURF_NORMALS_GRADIENT;
      else if ( strncmp(argv[i+1], "trian", 5) == 0 ) 
	mols->normals_algorithm = ISOSURF_NORMALS_TRIANGLE;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -tessellation option\" command; option should be cube or tetrahedron\n",
		 argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -MONOCOLOR */
    else if ( strncmp(argv[i], "-mono", 5) == 0 ) {
      double rgb;
      int argcList;
      const char **argvList;
      Tcl_SplitList(interp, argv[i+1], &argcList, &argvList);      
      for (j=0; j<argcList; j++) {
	if ( Tcl_GetDouble(interp, argvList[j], &rgb) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	mols->monocolor[j]=(float)rgb;
      }
      mols->monocolor[3]=0.5;
    }
    /* -SURFACETYPE */
    else if ( strncmp(argv[i], "-surf", 5) == 0 ) {
      if ( strncmp(argv[i+1], "mols", 4) == 0 ) 
	mols->surftype = MOLS_MOLSURF;
      else if ( strncmp(argv[i+1], "gap", 3) == 0 ) 
	mols->surftype = MOLS_GAP;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"xc_molsurf %s -surftype option\" command; option should be molsurf or gap\n",
		 argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -RESOLUTION */
    else if ( strncmp(argv[i], "-res", 4) == 0 ) {
      double res;
      if ( Tcl_GetDouble(interp, argv[i+1], &res) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      mols->d = (float) res;
    }
    /* -SMOOTHSTEPS */
    else if ( strncmp(argv[i], "-smooths", 8) == 0 ) {
      int nstep;
      if ( Tcl_GetInt(interp, argv[i+1], &nstep) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      mols->smooth_nstep = nstep;
    }
    /* -SMOOTHWEIGHT */
    else if ( strncmp(argv[i], "-smoothw", 8) == 0 ) {
      double weight;
      if ( Tcl_GetDouble(interp, argv[i+1], &weight) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      mols->smooth_weight = (float) weight;
    }
    else {
      Tcl_SetResult(interp, "Usage: xc_molsurf <tolg> -ident indent options ...", TCL_STATIC);
      return TCL_ERROR;
    }
  }
  /***************************************************************************/

  if (mols->type == MOLS_GAUSSIAN) {
    MolSurfGauss_Init( mols );
    mols->MolSurfFunc = MolSurfGauss;
  } 
  else if (mols->type == MOLS_EXP) {
    MolSurfExp_Init( mols );
    mols->MolSurfFunc = MolSurfExp;
  }
  else if (mols->type == MOLS_UNIGAUSS) 
    mols->MolSurfFunc = MolSurfGaussPlain;
  else if (mols->type == MOLS_UNIEXP) 
    mols->MolSurfFunc = MolSurfExpPlain;

  /********************/
  /* xcr.dim == 0,1,2 */
  /********************/
  /* define size[], cutoff, lowcoor */
  if ( xcr.dim < 3 || mols->surftype == MOLS_MOLSURF ) {
    mols->lightm.two_side_iso[0] = 0.0;

    mols->lowcoor[0] = ss.minX	- (max.r + 1.05 * mols->cutoff);
    mols->lowcoor[1] = ss.minY	- (max.r + 1.05 * mols->cutoff);
    mols->lowcoor[2] = ss.minZ	- (max.r + 1.05 * mols->cutoff);

    mols->size[0] = (ss.maxX - ss.minX) + 2.0*(max.r + 1.05*mols->cutoff);
    mols->size[1] = (ss.maxY - ss.minY) + 2.0*(max.r + 1.05*mols->cutoff);
    mols->size[2] = (ss.maxZ - ss.minZ) + 2.0*(max.r + 1.05*mols->cutoff);

    /*
      mols->size[0] = (max.x+mx) + max.r + 1.05*mols->cutoff - mols->lowcoor[0];
      mols->size[1] = (max.y+my) + max.r + 1.05*mols->cutoff - mols->lowcoor[1];
      mols->size[2] = (max.z+mz) + max.r + 1.05*mols->cutoff - mols->lowcoor[2];
    */

    for(i=0; i<3; i++) {    
      for(j=0; j<3; j++) {
	mols->vec[i][j] = 0.0;
	mols->isoexpand.rep_vec[i][j] = 0.0;
      }
      mols->vec[i][i] = mols->size[i];
    }
    if ( xcr.dim == 1 && mols->surftype == MOLS_GAP ) {
      /* insert code */
      ;
    } else if ( xcr.dim == 2 && mols->surftype == MOLS_GAP ) {
      /* insert code */
      ;
    }
  } else if ( xcr.dim == 3 ) {
    double (*vecP)[4];
    mols->lightm.two_side_iso[0] = 0.0;

    /* insert code */
    mols->lowcoor[0] = 0.0;
    mols->lowcoor[1] = 0.0;
    mols->lowcoor[2] = 0.0;
    if ( xcr.celltype == XCR_PRIMCELL ) vecP = vec.prim;
    else vecP = vec.conv;

    mols->size[0] = distdv( vecP[0] );
    mols->size[1] = distdv( vecP[1] );
    mols->size[2] = distdv( vecP[2] );

    for(i=0; i<3; i++) {
      for(j=0; j<3; j++) {
	mols->vec[i][j] = vecP[i][j];
	mols->isoexpand.rep_vec[i][j] = vecP[i][j];
      }
      mols->isoexpand.irepvec[i] = xcr.nunit[i];
    }
  }
  mols->i = iroundf(mols->size[0] / mols->d) + 1;
  mols->j = iroundf(mols->size[1] / mols->d) + 1;
  mols->k = iroundf(mols->size[2] / mols->d) + 1;    

  mols->molgrid = xcMallocTensor3f(mols->i, mols->j, mols->k);

  MolSurfGrid(mols, mols->molgrid);

  if(mols->isosurf_index < 0) mols->isosurf_index = NewIsoSurf();
  iso = FindIsoSurf( mols->isosurf_index );
  VPf.surface[mols->isosurf_index] = 1;

  iso->smooth_nstep  = mols->smooth_nstep;
  iso->smooth_weight = mols->smooth_weight;

  MarchCubeNew(iso, mols->molgrid, mols->vec, mols->lowcoor, 
	       mols->i, mols->j, mols->k, mols->isolevel, 
	       mols->tessellation_algorithm, mols->dispt.shademodel, mols->normals_algorithm);

  /* now add color */

  if ( mols->colorscheme == MOLS_ATOMIC ) {
    MolSurfGridAtomicColor( iso, mols->d );
  } 
  else if ( mols->colorscheme == MOLS_DISTANCE ) {
    ;
  } 
  else if ( mols->colorscheme == MOLS_MONOCHROME ) {
    ;
  }

  UpdateProjection();
  Togl_PostRedisplay(togl);

  /*free( (FREE_ARG) mols->molgrid );*/

  return TCL_OK;
}


void
MolSurfGridAtomicColor(ISOSURFACE *iso, float mols_d)
{
  register int i, j;
  float dis, mindis, mindis1, mindis2;
  int ind, ind1, ind2;
  float f0, f1, f2, d, d_2, div_1;

  iso->color = 
    (float (*)[4]) realloc(iso->color, 
			   (size_t)iso->nvertex * sizeof(float [4]));

  d   = 0.3;
  d_2 = 0.5 / d ;

  for (i=0; i<iso->nvertex; i++) {
    mindis = mindis1 = mindis2 = +9.9e99;
    ind  = 1;
    ind1 = 1;
    ind2 = 1;

    /*------------------------------------------------------------------------
      this was an idea, #0/#1/#2, will it give good results ???
      -----------------------------------------------------------------------*/
    /*
      find index #0
    */
    for (j=1; j<=natoms; j++) {
      dis = dist6( iso->vertex[i].x, xat[j]+mx,
		   iso->vertex[i].y, yat[j]+my,
		   iso->vertex[i].z, zat[j]+mz );
      if ( mindis > dis ) {
	ind    = j;
	mindis = dis;
      }
    }

    /*
      find index #1
    */
    for (j=1; j<=natoms; j++) {
      if (j == ind ) continue;
      dis = dist6( iso->vertex[i].x, xat[j]+mx,
		   iso->vertex[i].y, yat[j]+my,
		   iso->vertex[i].z, zat[j]+mz );
      if ( mindis1 > dis ) {
	ind1    = j;
	mindis1 = dis;
      }
    }

    /*
      find index #2
    */
    for (j=1; j<=natoms; j++) {
      if (j == ind || j == ind1 ) continue;
      dis = dist6( iso->vertex[i].x, xat[j]+mx,
		   iso->vertex[i].y, yat[j]+my,
		   iso->vertex[i].z, zat[j]+mz );
      if ( mindis2 > dis ) {
	ind2    = j;
	mindis2 = dis;
      }
    }
    /*
      if mindisX are with some range -> mix the colors

      let this parameter be d !!!
    */
    f0 = 1.0;
    f1 = 0.0;
    f2 = 0.0;

    if ( mindis1 - mindis < d ) {
      f0 = 0.5 + d_2 * (mindis1 - mindis);
      f1 = 1 - f0;
    }
    if ( mindis2 - mindis < d ) {
      f2    = 0.5 - d_2 * (mindis2 - mindis);
      div_1 = 1 / (f0 + f1 + f2);
      f0 *= div_1;
      f1 *= div_1;
      f2 *= div_1;
    }
    for (j=0; j<3; j++) {
      iso->color[i][j] = 
	f0*atm.col[ind][j] + f1*atm.col[ind1][j]+f2*atm.col[ind2][j];
      /*iso->color[i][j] = atcol[nat[ind]][j];*/
      /* so far !!! */
      iso->color[i][3] = 0.4;
      /*
	fprintf(stderr,"col: %f %f %f %f\n", 
	iso->color[i][0], iso->color[i][1], iso->color[i][2], 
	iso->color[i][3]);
      */
      /*fprintf(stderr,"f0: %10.5f; f1: %10.5f; f2: %10.5f\n", f0, f1, f2);*/
    }

    /*
      for(i=0; i<iso->nvertex; i++)
      fprintf(stderr,"col: %f %f %f %f\n", 
      iso->color[i][0], iso->color[i][1], iso->color[i][2], 
      iso->color[i][3]);
    */
  }
}


/* *** molgrid must not be pre-allocated */
int
MolSurfGrid(MOL_SURF *mols, float ***molgrid) 
{
  register int ix, iy, iz, i;
  register float dens, dx, dy, dz, dis, mindis;
  float tmpXYZ[3];
  /*float rad;*/

  if (mols->surftype == MOLS_MOLSURF) {
    /* take care of mols.res.$ !!!!!!!!!!!!!!!!!!! */
    dx = mols->size[0] / (float) (mols->i - 1);
    dy = mols->size[1] / (float) (mols->j - 1);
    dz = mols->size[2] / (float) (mols->k - 1);

    /* rad is cutoff criteria */
    /*rad = max.r + 1.05 * mols->cutoff;*/

    /* now made a loop over whole grid; take a XY plane and travel along Z 
     * direction
     */
    fprintf(stderr,"Molecular Surface: Grid generation\n");
    fprintf(stderr,"  Grid resolution: %d x %d x %d\n", 
	    mols->i, mols->j, mols->k);
    fprintf(stderr,"      Grid bounds: lowcoor=[%6.3f,%6.3f,%6.3f], size=[%6.3f,%6.3f,%6.3f]\n",
	    mols->lowcoor[0], mols->lowcoor[1], mols->lowcoor[2],
	    mols->size[0], mols->size[1], mols->size[2]);
    xcErrDebug("Beginning Grid calculation");
    for (ix=0; ix<mols->i; ix++) {
      tmpXYZ[0] = mols->lowcoor[0] + (float) ix * dx; 

      for (iy=0; iy<mols->j; iy++) {
	tmpXYZ[1] = mols->lowcoor[1] + (float) iy * dy;

	for (iz=0; iz<mols->k; iz++) {
	  tmpXYZ[2] = mols->lowcoor[2] + (float) iz * dz;	

	  dens   = 0.0;
	  mindis = 999.9;
	  for (i=1; i<=natoms; i++) {	
	    dis = DIST3( tmpXYZ[0] - (xat[i]+mx),
			 tmpXYZ[1] - (yat[i]+my),
			 tmpXYZ[2] - (zat[i]+mz) );
	    if ( mols->type != MOLS_DISTFUNC ) {
	      dens += mols->MolSurfFunc(mols, nat[i], dis);
	    } else {
	      if ( dis < mindis ) {
		mindis = dis;
		dens   = dis;
	      }
	    }
	  }
	  molgrid[ix][iy][iz] = dens;
	}

      }
      xcErrDebug("step YZ calculation performed");
    } /* ix */
    xcErrDebug("Grid calculation is Done");
    /*
      mols->natm = 0;
      for (i=1; i<=natoms; i++) {	
      xl = tmpXYZ[0] - (xat[i]+mx);
      yl = tmpXYZ[1] - (yat[i]+my);
      dis = sqrt( xl * xl  +  yl * yl );
      if ( dis < rad ) {
      mols->atmindex[mols->natm] = i;
      mols->natm++;
      }
      }

      for (iz=0; iz<mols->k; iz++) {
      tmpXYZ[2] = mols->lowcoor[2] + (float) iz * dz;

      dens = 0.0;
      mindis = 2.0 * ss.maxX;
      for (i=0; i<mols->natm; i++) {
      dis = dist3f( tmpXYZ[0] - (xat[ mols->atmindex[i] ]+mx),
      tmpXYZ[1] - (yat[ mols->atmindex[i] ]+my),
      tmpXYZ[2] - (zat[ mols->atmindex[i] ]+mz) );
      if ( mols->type != MOLS_DISTFUNC ) {
      dens += mols->MolSurfFunc(mols, nat[ mols->atmindex[i] ], dis);
      } else {
      if ( dis < mindis ) {
      mindis = dis;
      dens   = dis;
      }
      }
      }
      molgrid[ix][iy][iz] = dens;
      } */

  } else {
    /******************************/
    /* mols->surftype == MOLS_GAP */
    /******************************/
    register int ii, jj, kk, idn, iup, jdn, jup, kdn, kup;
    float x[3], a, b, c;
    int size, *atn;
    double (*fcoor)[3], (*vP)[4], dii, djj, dkk;

    size  = xcr.natr;
    if ( xcr.celltype == XCR_PRIMCELL) {
      fcoor = xcr.prim_fcoor;      
      atn   = xcr.prim_nat;
      vP    = vec.prim;
    }
    else if ( xcr.celltype == XCR_CONVCELL ) {
      size *= xcr.ncell;      
      fcoor = xcr.conv_fcoor;
      atn   = xcr.conv_nat;
      vP    = vec.conv;
    }

    /* now made a loop over whole grid */
    for (ix=0; ix<mols->i; ix++) {
      dx = (float)ix / (float)(mols->i-1);
      for (iy=0; iy<mols->j; iy++) {
	dy = (float)iy / (float)(mols->j-1);
	for (iz=0; iz<mols->k; iz++) {
	  dz = (float)iz / (float)(mols->k-1);

	  dens = 0.0;
	  mindis = 9.9e+99;
	  for (i=0; i<size; i++) {
	    idn = jdn = kdn = -1;
	    iup = jup = kup = +1;
	    /*
	      if ( fcoor[i][0] < 0.3 ) {
	      idn=-1; iup=0;
	      } else if ( fcoor[i][0] > 0.7 ) {
	      idn=0; iup=1;
	      }
	      if ( fcoor[i][1] < 0.3 ) {
	      jdn=-1; jup=0;
	      } else if ( fcoor[i][1] > 0.7 ) {
	      jdn=0; jup=1;
	      }
	      if ( fcoor[i][2] < 0.3 ) {
	      kdn=-1; kup=0;
	      } else if ( fcoor[i][2] > 0.7 ) {
	      kdn=0; kup=1;
	      }
	    */

	    /* in cell 000 */
	    for (ii=idn; ii<=iup; ii++) {
	      dii = (double) ii;
	      for (jj=jdn; jj<=jup; jj++) {
		djj = (double) jj;
		for (kk=kdn; kk<=kup; kk++) {
		  dkk = (double) kk;
		  a = dx - (dii + fcoor[i][0]);
		  b = dy - (djj + fcoor[i][1]);
		  c = dz - (dkk + fcoor[i][2]);

		  /* convert tmpXYZ to Cartesian coordinates !!! */
		  x[0] = a*vP[0][0] + b*vP[1][0] + c*vP[2][0];
		  x[1] = a*vP[0][1] + b*vP[1][1] + c*vP[2][1];
		  x[2] = a*vP[0][2] + b*vP[1][2] + c*vP[2][2];
		  dis = distfv( x );
		  if ( mols->type != MOLS_DISTFUNC ) {
		    dens += mols->MolSurfFunc(mols, atn[i], dis);
		  } else {
		    if ( dis < mindis ) {
		      mindis = dis;
		      dens   = dis;
		    }
		  }
		}
	      }
	    }
	  }
	  molgrid[ix][iy][iz] = dens;
	} /* iz */
      }
    }
  }
  return XC_OK;
}


