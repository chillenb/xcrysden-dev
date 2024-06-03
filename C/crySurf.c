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
* Source: $XCRYSDEN_TOPDIR/C/crySurf.c
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
#include "vector.h"
#include "xcfunc.h"

#define NUM_OF_POINTS(old_n,degree) ( (old_n - 1)*(degree) + 1 )

extern void SetUnitCellCage( double vec[4][4], CellCage *cage );
extern NEW_WIN_CONTEXT *FindWinContextByTogl(struct Togl *togl);
extern MOL_SURF *FindMolSurf(int index);
extern int cryReadSurfARGS(ClientData clientData, Tcl_Interp *interp,
			   int argc, char *argv[], NEW_WIN_CONTEXT *wc, 
			   MOL_SURF *mols, struct READSURF_OPT *opt);
extern void MarchCubeNew(ISOSURFACE *iso,
			 float ***gridvalue,
			 float vec[][3],
			 float origin[],
			 int ix, int iy, int iz,
			 float isolevel, int algorithm, int shade_model, int normals_model);
extern void xcSurfSmoothing( ISOSURFACE *iso, float isolevel );

extern float ***cryGeneralGridRegPeriodInterpolator_f_LCASI3D(int n[3], int degree[3], float ***src);


static void Set_Orig_and_Size( NEW_WIN_CONTEXT *wc, MOL_SURF *m, 
			       struct Togl *togl );
/*
  -TYPE 
  -FS
  -RADIU 
  -LEVEL 
  -CUTOFF
  -COLORSHEME 
  -DRAWSTYLE 
  -TRANSPARENT 
  -SHADEMODEL 
  -MONOCOLOR 
  -FRONTMONOCOLOR 
  -BACKMONOCOLOR 
  -SURFACETYPE 
  -RESOLUTION 
  -SMOOTHSTEPS 
  -SMOOTHWEIGHT
  -FRONTFACE
  -REVERTNORMALS

  Sample READSURF_OPT struct:

  struct READSURF_OPT *ReadSurf_Options = {
  {"-type",           XC_ALLOWED},
  {"-fs",             XC_ALLOWED},
  {"-render",         XC_ALLOWED},
  {"-radiu", 	    XC_ALLOWED},
  {"-level", 	    XC_ALLOWED},
  {"-cutoff"	    XC_ALLOWED},
  {"-colorsheme",     XC_ALLOWED},
  {"-drawstyle",      XC_ALLOWED},
  {"-transparent",    XC_ALLOWED},
  {"-shademodel",     XC_ALLOWED},
  {"-monocolor",      XC_ALLOWED},
  {"-frontmonocolor", XC_ALLOWED},
  {"-backmonocolor",  XC_ALLOWED},
  {"-surfacetype",    XC_ALLOWED},
  {"-resolution",     XC_ALLOWED},
  {"-smoothsteps",    XC_ALLOWED},
  {"-smoothweight",   XC_ALLOWED},
  {"-frontface",      XC_ALLOWED},
  {"-revertnormals",  XC_ALLOWED},
  {NULL, NULL}
  }

  *** option: -fs { 
  -gridindex    <index> \
  -gridsubindex <subindex> \
  -bandindex    <index> \
  -celltype     para|bz \
  -cropbz       0|1  \    
  -displaycell  0|1 \		    ?  
  -celldisplaytype     wire|rod|solid  ?
  -interpolationdegree {n1 n2 n3}
  -wirecellcolor  {r g b a}
  -solidcellcolor {r g b a}
  }

*/

/*****************************************************************************
cry_surf <togl>   -ident        <identifier>
-type         gauss|exp|unigauss|uniexp|fermisurface \
-fs           {options}
-radius       cov|VdW \
-level        <isoleve> \
-cutoff       <cutoff> \
-colorscheme  atomic|monochrome \
-drawstyle    solid|wire|dot \
-transparent  0|1 \
-shademodel   smooth|flat \
-monocolor      {<r> <g> <b>} \
-frontmonocolor {<r> <g> <b>} \
-backmonocolor" {<r> <g> <b>} \
-surfacetype  molsurf|gap
-resolution   <x> \
-smoothsteps  <steps> \
-smoothweight <weight>
-frontface    CW|CCW
-revertnormals 0|1

ADD:

cry_surfreg <togl>

cry_surfconfig <tolg>   -ident <identifier> \
-render 0|1 \
-level <isolevel >
-colorscheme atomic|monochrome \
-drawstyle solid|wire|dot \
-transparent  0|1 \
-shademodel   smooth|flat
-monocolor    {<r> <g> <b>} \			  
-frontmonocolor {<r> <g> <b>} \
-backmonocolor" {<r> <g> <b>} \
-surfacetype  molsurf|gap \
-smoothsteps  <steps> \
-smoothweight <weight>
-frontface    CW|CCW
-revertnormals 0|1



preveri se pointerje pri iso. V zacetku jih postavi na NULL, 
potem pa v MarchCube, ce niso NULL jih najprej free();
*****************************************************************************/

int FindSurfaceOfWindow( NEW_WIN_CONTEXT *wc, MOL_SURF *mols )
{ 
  register int i;
  for(i=0; i<wc->VPf.nsurface; i++)
    if ( wc->VPf.surfaceInd[i] == mols->isosurf_index )
      return i;

  return -1;
}


int 
CRY_SurfRegCmd(ClientData clientData,Tcl_Interp *interp,
	       int argc, char *argv[])
{
  int ident;
  char *result = (char *) Tcl_Alloc((size_t) 10 * sizeof(char));
  struct Togl *togl;
  NEW_WIN_CONTEXT *wc;
  MOL_SURF *m;

  /* find togl associated with toglName */
  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  wc = FindWinContextByTogl( togl );

  ident=NewMolSurf();
  m = FindMolSurf( ident );
  if( m->isosurf_index < 0) m->isosurf_index = NewIsoSurf();

  sprintf(result,"%d",ident);
  Tcl_SetResult(interp, result, TCL_DYNAMIC);

  wc->VPf.surface    = 
    (int *) realloc( wc->VPf.surface, 
		     (size_t) (wc->VPf.nsurface+1) * sizeof(int) );
  wc->VPf.surfaceInd = 
    (int *) realloc( wc->VPf.surfaceInd, 
		     (size_t) (wc->VPf.nsurface+1) * sizeof(int) );
  wc->VPf.surfacePtr = 
    (void **) realloc( wc->VPf.surfacePtr, 
		       (size_t) (wc->VPf.nsurface+1) * sizeof(void*) );

  wc->VPf.surface[wc->VPf.nsurface]      = 0;
  wc->VPf.surfaceInd[wc->VPf.nsurface]   = m->isosurf_index;
  wc->VPf.surfacePtr[wc->VPf.nsurface++] = (void *) m;
  return TCL_OK;
}

int
CRY_SurfConfigCmd(ClientData clientData,Tcl_Interp *interp,
		  int argc, char *argv[])
{
  register int        i;
  short int           crop_bz;
  int                 s_nstep, interp_n[3]; 
  float               s_weight, s_level;
  struct Togl         *togl;
  FS_DATA             fs;
  NEW_WIN_CONTEXT     *wc;
  MOL_SURF            *mols;
  ISOSURFACE          *iso;
  struct READSURF_OPT opt[] = {
    {"-type",           XC_FORBIDDEN},
    {"-fs",             XC_ALLOWED},
    {"-render",         XC_ALLOWED},
    {"-radiu", 	        XC_FORBIDDEN},
    {"-level", 	        XC_ALLOWED},
    {"-cutoff",	        XC_FORBIDDEN},
    {"-colorsheme",     XC_FORBIDDEN},
    {"-drawstyle",      XC_ALLOWED},
    {"-transparent",    XC_ALLOWED},
    {"-shademodel",     XC_ALLOWED},
    {"-monocolor",      XC_ALLOWED},
    {"-frontmonocolor", XC_ALLOWED},
    {"-backmonocolor",  XC_ALLOWED},
    {"-surfacetype",    XC_FORBIDDEN},
    {"-resolution",     XC_FORBIDDEN},
    {"-smoothsteps",    XC_ALLOWED},
    {"-smoothweight",   XC_ALLOWED},
    {"-frontface",      XC_ALLOWED},
    {"-revertnormals",  XC_ALLOWED},
    {NULL, 0}
  };

  /* find togl associated with toglName */
  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  wc   = FindWinContextByTogl( togl );

  /* -IDENTIFIER */
  if ( strncmp(argv[2], "-iden", 5) == 0 ) {
    int ident;
    if ( Tcl_GetInt(interp, argv[3], &ident) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing %s %s %s %s ...", argv[2], argv[0], argv[1], argv[2], argv[3]);
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
      Tcl_SetResult(interp, "cry_surf must be called before cry_surfconfig!!!",TCL_STATIC);
      return TCL_ERROR;
    }
  } else {
    Tcl_SetResult(interp, "Usage: cry_surfconfig <tolg> -ident indent options ...",TCL_STATIC);
    return TCL_ERROR;
  }
  fs.gridindex    = mols->fs.gridindex;
  fs.gridsubindex = mols->fs.gridsubindex;
  fs.bandindex    = mols->fs.bandindex;
  fs.celltype     = mols->fs.celltype;
  fs.cropbz       = mols->fs.cropbz;

  s_level  = mols->isolevel;
  s_nstep  = mols->smooth_nstep;
  s_weight = mols->smooth_weight;
  crop_bz  = mols->fs.cropbz;

  interp_n[0]     = mols->interp_n[0];
  interp_n[1]     = mols->interp_n[1];
  interp_n[2]     = mols->interp_n[2];

  /* read COMM LINE ARGS ****************************************/
  if (cryReadSurfARGS(clientData, interp, argc, argv, wc, mols, opt) ==
      TCL_ERROR) {
    return TCL_ERROR;
  }
  /**************************************************************/

  iso = FindIsoSurf( mols->isosurf_index );
  iso->smooth_nstep  = mols->smooth_nstep;
  iso->smooth_weight = mols->smooth_weight;

  if ( fs.gridindex    != mols->fs.gridindex ||
       fs.gridsubindex != mols->fs.gridsubindex ||
       fs.bandindex    != mols->fs.bandindex ||
       fs.celltype     != mols->fs.celltype ) {
    Tcl_SetResult(interp, "can't change gridindex, gridsubindex, bandindex and celltype by cry_surfconfig command", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( interp_n[0] != mols->interp_n[0] ||
       interp_n[1] != mols->interp_n[1] ||
       interp_n[2] != mols->interp_n[2] ) {
    /* LCASI intepolation */
    int _n[3];
    /*
      _n[0] = mols->i - 1;
      _n[1] = mols->j - 1;
      _n[2] = mols->k - 1;
    */
    _n[0] = mols->i;
    _n[1] = mols->j;
    _n[2] = mols->k;
    if ( mols->interp_molgrid != (float***) NULL ) {
      xcFree_Tensor3f (mols->interp_molgrid);
      mols->interp_molgrid = (float***) NULL;
    }
    /*mols->interp_molgrid = cryRegPeriodInterpolator_f_LCASI3D(_n, mols->interp_n, mols->molgrid);*/
    if ( mols->interp_n[0] * mols->interp_n[1] * mols->interp_n[2] > 1 ) {

      /* mols->interp_molgrid = cryGeneralGridRegPeriodInterpolator_f_LCASI3D (_n, mols->interp_n, mols->molgrid); */

      mols->interp_molgrid = general_grid_fft_interpolator_tensor3f(_n, mols->interp_n, mols->molgrid);
      /* mols->interp_molgrid = general_grid_shankland_interpolator_tensor3f(_n, mols->interp_n, mols->molgrid); */

      MarchCubeNew(iso, mols->interp_molgrid, mols->vec, mols->lowcoor, 
		   NUM_OF_POINTS(mols->i, mols->interp_n[0]), 
		   NUM_OF_POINTS(mols->j, mols->interp_n[1]), 
		   NUM_OF_POINTS(mols->k, mols->interp_n[2]), mols->isolevel,
		   ISOSURF_MARCHING_CUBES, ISOSURF_SHADE_SMOOTH, ISOSURF_NORMALS_GRADIENT);
    } else {
      MarchCubeNew(iso, mols->molgrid, mols->vec, mols->lowcoor, 
		   mols->i, mols->j, mols->k, mols->isolevel,
		   ISOSURF_MARCHING_CUBES, ISOSURF_SHADE_SMOOTH, ISOSURF_NORMALS_GRADIENT);
    }
    /*xcSurfSmoothing( iso, mols->isolevel );*/
    if ( mols->fs.cropbz && mols->fs.celltype == XCR_BZ )
      CropBz_New( wc, iso );
  } else if ( interp_n[0] * interp_n[1] * interp_n[2] > 1 &&
	      s_level != mols->isolevel ) {
    MarchCubeNew(iso, mols->interp_molgrid, mols->vec, mols->lowcoor, 
		 NUM_OF_POINTS(mols->i, mols->interp_n[0]), 
		 NUM_OF_POINTS(mols->j, mols->interp_n[1]), 
		 NUM_OF_POINTS(mols->k, mols->interp_n[2]), mols->isolevel,
		 ISOSURF_MARCHING_CUBES, ISOSURF_SHADE_SMOOTH, ISOSURF_NORMALS_GRADIENT);
    /*xcSurfSmoothing( iso, mols->isolevel );*/
    if ( mols->fs.cropbz && mols->fs.celltype == XCR_BZ )
      CropBz_New( wc, iso );
  } else if ( s_level != mols->isolevel) {
    MarchCubeNew(iso, mols->molgrid, mols->vec, mols->lowcoor, 
		 mols->i, mols->j, mols->k, mols->isolevel,
		 ISOSURF_MARCHING_CUBES, ISOSURF_SHADE_SMOOTH, ISOSURF_NORMALS_GRADIENT);
    /*xcSurfSmoothing( iso, mols->isolevel );*/
    if ( mols->fs.cropbz && mols->fs.celltype == XCR_BZ )
      CropBz_New( wc, iso );
  }
  else if ( s_nstep  != iso->smooth_nstep  || s_weight != iso->smooth_weight ) {
    /*
      xcSurfSmoothing( iso, mols->isolevel );
      if ( mols->fs.cropbz && mols->fs.celltype == XCR_BZ )
      CropBz_New( wc, iso );
    */
    ;
  }
  else if ( crop_bz != mols->fs.cropbz && mols->fs.celltype == XCR_BZ ) {
    if ( !mols->fs.cropbz && mols->fs.celltype == XCR_BZ ) {
      /* previous cropping has froget about croped triangles -->
	 retriamgulate 
      */
      MarchCubeNew(iso, mols->molgrid, mols->vec, mols->lowcoor, 
		   mols->i, mols->j, mols->k, mols->isolevel,
		   ISOSURF_MARCHING_CUBES, ISOSURF_SHADE_SMOOTH, ISOSURF_NORMALS_GRADIENT);
    }
    /*xcSurfSmoothing( iso, mols->isolevel );*/
    if ( mols->fs.cropbz && mols->fs.celltype == XCR_BZ )
      CropBz_New( wc, iso );
  }        

  /* HERE GOES THE CROPPING algorithm */
  /* this is the safer, although time demanding */

  /*
    if ( mols->fs.cropbz && mols->fs.celltype == XCR_BZ )
    CropBz_New( wc, iso );
  */

  if ( !mols->fs.cropbz )
    for(i=0; i<iso->ntriangl; i++)
      iso->triangl_status[i] = 1;

  Togl_PostRedisplay(togl);

  return TCL_OK;
}


int 
CRY_SurfCmd(ClientData clientData,Tcl_Interp *interp,
	    int argc, char *argv[])
{
  register int        i;
  struct Togl         *togl; 
  NEW_WIN_CONTEXT     *wc;
  MOL_SURF            *mols;
  ISOSURFACE          *iso;
  struct READSURF_OPT opt[] = {
    {"-type",           XC_ALLOWED},
    {"-fs",             XC_ALLOWED},
    {"-render",         XC_ALLOWED},
    {"-radiu", 	        XC_FORBIDDEN},
    {"-level", 	        XC_ALLOWED},
    {"-cutoff",	        XC_FORBIDDEN},
    {"-colorsheme",     XC_FORBIDDEN},
    {"-drawstyle",      XC_ALLOWED},
    {"-transparent",    XC_ALLOWED},
    {"-shademodel",     XC_ALLOWED},
    {"-monocolor",      XC_ALLOWED},
    {"-frontmonocolor", XC_ALLOWED},
    {"-backmonocolor",  XC_ALLOWED},
    {"-surfacetype",    XC_FORBIDDEN},
    {"-resolution",     XC_FORBIDDEN},
    {"-smoothsteps",    XC_ALLOWED},
    {"-smoothweight",   XC_ALLOWED},
    {"-frontface",      XC_ALLOWED},
    {"-revertnormals",  XC_ALLOWED},
    {NULL, 0},
  };

  /* find togl associated with toglName */
  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  wc = FindWinContextByTogl( togl );

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
      snprintf(rss, sizeof(rss),"couldn't find Surface %d !!!\n",
	       ident);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  } else {
    char rss[1024];
    snprintf(rss, sizeof(rss),"Usage: %s <tolg> -ident indent ...", argv[0]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* read COMMAND LINE ARGS ***************************************/
  if (cryReadSurfARGS(clientData, interp, argc, argv, wc, mols, opt) ==
      TCL_ERROR) {
    return TCL_ERROR;
  }
  /*************************************************************/

  if (mols->type == MOLS_FERMISURFACE) {
    struct DATAGRID     *grid;
    CellCage            *cage;
    grid = FindDataGrid( mols->fs.gridindex );
    fsReadBand(grid, mols);    
    /* generate lattice-cage */
    cage = (CellCage *) realloc(wc->recprim_cage, sizeof(CellCage));   
    SetUnitCellCage ( vec.recprim, cage );
    wc->recprim_cage = cage;
  } 
  else {
    Tcl_SetResult(interp, "Sorry, so far just fermisurface type surfaces supported by cry_surf command", TCL_STATIC);
    return TCL_ERROR;
  }



  iso = FindIsoSurf( mols->isosurf_index );
  i = FindSurfaceOfWindow( wc, mols );
  if ( i < 0 ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "couldn't find surface #%d on window %s !!!\n",
	     mols->isosurf_index, argv[1]);  
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  wc->VPf.surface[i] = 1;

  iso->smooth_nstep  = mols->smooth_nstep;
  iso->smooth_weight = mols->smooth_weight;

  fprintf(stderr,"FS: interp_n[] == %d, %d, %d\n", 
	  mols->interp_n[0], mols->interp_n[1], mols->interp_n[2] );

  /* t.k. */
  if ( mols->interp_n[0] * mols->interp_n[1] * mols->interp_n[2] > 1 ) {
    int _n[3];
    _n[0] = mols->i;
    _n[1] = mols->j;
    _n[2] = mols->k;

    if ( mols->interp_molgrid != (float***) NULL ) {
      xcFree_Tensor3f (mols->interp_molgrid);
      mols->interp_molgrid = (float***) NULL;
    }

    /*mols->interp_molgrid = cryRegPeriodInterpolator_f_LCASI3D(_n, mols->interp_n, mols->molgrid);*/

    /* mols->interp_molgrid = cryGeneralGridRegPeriodInterpolator_f_LCASI3D (_n, mols->interp_n, mols->molgrid); */

    mols->interp_molgrid = general_grid_fft_interpolator_tensor3f(_n, mols->interp_n, mols->molgrid);

    /* mols->interp_molgrid = general_grid_shankland_interpolator_tensor3f(_n, mols->interp_n, mols->molgrid); */

    MarchCubeNew(iso, mols->interp_molgrid, mols->vec, mols->lowcoor, 
		 NUM_OF_POINTS(mols->i, mols->interp_n[0]), 
		 NUM_OF_POINTS(mols->j, mols->interp_n[1]), 
		 NUM_OF_POINTS(mols->k, mols->interp_n[2]), mols->isolevel,
		 ISOSURF_MARCHING_CUBES, ISOSURF_SHADE_SMOOTH, ISOSURF_NORMALS_GRADIENT);
  } else {
    /********/
    MarchCubeNew(iso, mols->molgrid, mols->vec, mols->lowcoor, 
		 mols->i, mols->j, mols->k, mols->isolevel,
		 ISOSURF_MARCHING_CUBES, ISOSURF_SHADE_SMOOTH, ISOSURF_NORMALS_GRADIENT);
  }

  /* HERE GOES THE CROPPING algorithm */
  /* this is the safer, although time demanding */

  if ( mols->fs.cropbz && mols->fs.celltype == XCR_BZ )
    CropBz_New( wc, iso );  

  if ( !mols->fs.cropbz )
    for(i=0; i<iso->ntriangl; i++)
      iso->triangl_status[i] = 1;

  wc->VPf.stropened = 1;
  Set_Orig_and_Size( wc, mols, togl );
  crySetProjection( wc, togl );

  Togl_PostRedisplay(togl);

  return TCL_OK;
}


static void
Set_Orig_and_Size( NEW_WIN_CONTEXT *wc, MOL_SURF *m, struct Togl *togl )
{
  int i;

  if ( m->fs.celltype == XCR_BZ ) {
    for (i=0; i<3; i++)
      wc->MVf.o_shift[i] = 0.0;

    wc->ss.minX = m->lowcoor[0];
    wc->ss.maxX = -wc->ss.minX;	  

    wc->ss.minY = m->lowcoor[1];
    wc->ss.maxY = -wc->ss.minY;	  

    wc->ss.minZ = m->lowcoor[2];
    wc->ss.maxZ = -wc->ss.minZ;      
  } else {
    /*******************/
    /* XCR_PARAPIPEDAL */
    /*******************/
    float size, orig[3] = {0.0, 0.0, 0.0};
    size = DetermineParapipedSize(m->vec[0], m->vec[1], m->vec[2], orig)/2.0;
    for (i=0; i<3; i++) {
      wc->MVf.o_shift[i] = -(m->vec[0][i] + m->vec[1][i] + m->vec[2][i])/2.0;
    }
    /*
      wc->ss.minX =  wc->MVf.o_shift[0];
      wc->ss.maxX = -wc->MVf.o_shift[0];

      wc->ss.minY =  wc->MVf.o_shift[1];
      wc->ss.maxY = -wc->MVf.o_shift[1];

      wc->ss.minZ =  wc->MVf.o_shift[2];
      wc->ss.maxZ = -wc->MVf.o_shift[2];
    */

    wc->ss.minX = -size;
    wc->ss.maxX =  size;

    wc->ss.minY = -size;
    wc->ss.maxY =  size;

    wc->ss.minZ = -size;
    wc->ss.maxZ =  size;
  }
}

