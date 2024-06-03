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
* Source: $XCRYSDEN_TOPDIR/C/xcIsoSurf.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <togl.h>
#include <tk.h>
#include <GL/gl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "xcGLparam.h"
#include "struct.h"
#include "isosurf.h"
#include "xcfunc.h"


extern Options3D is;

/* --- extern function prototypes ---*/

PLANEVERTEX ***plvertex = NULL;
/*****************************************************************************
 * display-atributes for isosurface & isoplane                             */
ISO_ATTRIB isoDisp = { ISOSURF_SOLID, GL_SMOOTH, ISOSURF_TRANSP_OFF, GL_TRUE };
ISO_ATTRIB isoplaneDisp[MAX_ISOOBJECTS];

/* when we perform xc_iso* commands, we must follow some order of commands, 
 * isostate variable trace that
 */
ISOSTATE isostate = { 0, {0, 0, 0, 0}, {0, 0, 0, 0}, 
		      {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, 
		      {
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			 0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			 0, 0, 0, 0, 0, 0, 0, 0}
		      }, 0, 0, ISO_NULL };


void  xcIsoError(void);
static void  InitIsoDataArr(void);
static void  InitIsoExpand(void);
static float SetIsosurf_VPf(void);

/*****************************************************************************
 * IN THIS FILE THE FOLLOWING COMMANDS ARE IMPLEMENTED:: 
 * ---------------------------------------------------
 * 
 * xc_iso init|finish|(end isostack|isosign)|grid|save| \
 *        minvalue|maxvalue|(isolevel <value>)|polygonise ?-algorithm cubes|tetrahedrons? ?-shademodel smoth|flat? ?-normals gradient|triangles?| \
 *        (smoothsteps <nstep>)|smoothweight <weight>)| \
 *        (get smoothsteps|smoothweight)|smoothing \
 *        (isoplane <type> <colorBase> <scalefunc> <what> ?<islide>?)| \
 *        interpolate <degree> |\
 *        (isoplaneconfig <object> -isoplanemin  <min> \
 *                                 -isoplanemax  <max> \
 *                                 -isolinecolor {monocolor r g b}|polycolor \
 *                                 -isolinewidth width \
 *                                 -isolinedash  nodash | negdash | fulldash \
 *                                 -isolinenlevels <nlevels> \
 *                                 -isoplanelighting   0/1)
 *
 * xc_isostack <nstack>|(<frame_in_stack3> <stack=[2-0]> \
 *                       <n_of_frames_in_stack[2-0]>)
 * 
 * xc_isosign <stack_level> <frame> <frame_sign>
 *
 * xc_isofiles bin_vertex_filename bin_filename filename0 \
 *                                              ?filename1? ?filename2? ...."
 *
 * xc_isopoints <object> 2D|3D points
 *
 * xc_isodata (<frame[i]_range>,i=ISODATA_MAXSTACK,0,-1)
 *
 * xc_isosurf <mesaWin> -drawstyle solid|wire|dot -shademodel smoth|flat \
 *                      -transparency on|off -isosurf none|xxx
 * xc_isoplane <mesaWin> <object> \
 *                   -planetype  colorplane|isoline|both \
 *                   -transparency on|off
 *                   -render now|after
 *
 * xc_isoexpand <mesaWin> <object> \
 *                        -repeattype default|convcell|primcell \
 *                        -shape      default|parapipedal|hexagonal \
 *                        -expand     whole|{nx ny nz}|... 
 *                        -render     now|after
 *****************************************************************************/


/*****************************************************************************
 * auxilary functions                                                        */
void
xcIsoError(void)
{
  int i, j;
  if ( isostate.gridvertex_malloc ) {
    xcFree_GRIDVERTEX( gridvertex );
    isostate.gridvertex_malloc = 0;
  }

  for(i=ISOOBJ_BASE; i<MAX_ISOOBJECTS; i++) {
    if ( isostate.plvertex_malloc[i] ) {
      xcFree_PLANEVERTEX( plvertex[i] );
      isostate.plvertex_malloc[i] = 0;
    }
    if ( isostate.isoline2D_malloc[i] ) {
      for (j=0; j<ISOLINE_MAXLEVEL; j++)
	xcFree_LINE( isoline2D[i].segment[j] );
      isostate.isoline2D_malloc[i] = 0;
    }
  }

  for (i=0; i<MAX_ISOSURFLEVELS; i++) {
    if ( isostate.vertex_malloc[i] ) {
      free( vertex[i] );
      isostate.vertex_malloc[i] = 0;
    }
    if ( isostate.triangl_malloc[i] ) {
      free( triangles[i] );
      isostate.triangl_malloc[i] = 0;
    }  
    if ( isostate.tri2verIN_malloc[i] ) {
      free( tri2verIN[i] );
      isostate.tri2verIN_malloc[i] = 0;
    }  
  }

  if ( isostate.bin_file_open ) {
    fclose(isodata.bin_fp);
    isostate.bin_file_open = 0;
  }
  if ( isostate.bin_vertex_file_open ) {
    fclose(isodata.bin_vertex_fp);
    isostate.bin_vertex_file_open = 0;
  }

  /* if an error occure, that means that we have error in Tcl/Tk
   * script. It is very dificcult to recover correctly, so we will initialize
   * all arrays in ISOSURF_DATA structure
   */
  InitIsoDataArr();

  for (i=0; i<MAX_ISOOBJECTS; i++) {
    /* THIS IS NEEDED TO RENDER ISOSURFACE */
    VPf.isosurf[i] = 0;
    /* THIS IS NEEDED TO RENDER COLORPLANE */
    VPf.colorplane[i] = 0;
    /* THIS IS NEEDED TO RENDER ISOLINE */
    VPf.isoline[i] = 0;
  }

}


/* initialize all arrays in isodata structure */
static void
InitIsoDataArr(void)
{
  int i, j, k;
  /*    static int error_malloced = 0; */

  for (i=0; i<ISODATA_MAXFILES; i++)
    for(j=0; j<ISODATA_MAXSTACK; j++)
      isodata.nframe[i][j] = 0;

  for (i=0; i<ISODATA_MAXSTACK; i++)
    for (j=0; j<ISODATA_MAXFRAME; j++)
      isodata.framesign[i][j] = 0;

  for (i=0; i<ISODATA_MAXFILES; i++)
    isodata.isofiles[i] = "";

  for (i=ISOOBJ_BASE; i<MAX_ISOOBJECTS; i++)
    for (j=0; j<3; j++) {
      isodata.vec[i][j].x = 0.0;
      isodata.vec[i][j].y = 0.0;
      isodata.vec[i][j].z = 0.0;
      for (k=0; k<4; k++)
	isodata.points[i][k][j] = 0.0;
    }

  /* DEFAULT for surface smoothing algorithm */
  isodata.smoothsteps  = 0;
  isodata.smoothweight = 0.5;

  /*    if (!error_malloced) { */
  /*      isodata.error = (char *) malloc ( sizeof(char) * 512 ); */
  /*      error_malloced = 1; */
  /*    } */
}


static void
InitIsoExpand(void)
{
  int i, j, k;

  for (k=0; k<MAX_ISOOBJECTS; k++) {
    for (i=0; i<4; i++)
      for (j=0; j<4; j++)
	isoexpand[k].rep_vec[i][j] = 0.0;

    for (i=0; i<3; i++)
      for (j=0; j<3; j++)
	isoexpand[k].transl[i][j] = 0;

    isoexpand[k].shape  = COM_ISOEXPAND_SHAPE_DEFAULT;
    isoexpand[k].expand = COM_ISOEXPAND_EXPAND_WHOLE;

    for (i=0; i<3; i++)
      isoexpand[k].irepvec[i] = 1;
  }
}


/*
  When rendering isosurface I have a 3D mesh of points, so I can render some 
  isoplanes as well; ISOOBJ_PLANE1, ISOOBJ_PLANE2, ISOOBJ_PLANE3 are 
  particular interesting -> they represent all three faces for some 
  paralepiped; check if it is one of such plane
*/
int
IsIsoPlane123( int type ) {
  if ( type == ISOOBJ_PLANE1 ||
       type == ISOOBJ_PLANE2 ||
       type == ISOOBJ_PLANE3  )
    return 1;
  else
    return 0;
}


/******************************************************************************
 * XC_IsoCmd --> inplementation of 'xc_iso' custom Tcl command 
 * --------------------- 
 * Usage: xc_iso init|finish|(end isostack|isosign)|grid|save| \
 *        minvalue|maxvalue|(isolevel <value>)|polygonise ?-algorithm cubes|tetrahedrons? ?-shademodel smoth|flat? ?-normals gradient|triangles?| \
 *        (smoothsteps <nstep>)|smoothweight <weight>)| \
 *        (get smoothsteps|smoothweight)|smoothing \
 *        (isoplane <type> <colorBase> <scalefunc> <what> ?<islide>?)| \
 *        interpolate <degree> |\
 *        (isoplaneconfig <object> -isoplanemin  <min> \
 *                                 -isoplanemax  <max> \
 *                                 -isolinecolor {monocolor r g b}|polycolor \
 *                                 -isolinewidth width \
 *                                 -isolinedash  nodash | negdash | fulldash \
 *                                 -isolinenlevels <nlevels> \
 *                                 -isoplanelighting   0/1)
 *
 */
int
XC_IsoCmd(ClientData clientData, Tcl_Interp *interp,
	  int argc, const char *argv[])
{ 

  /* printf("%d;   %s %s %s %s\n",argc,argv[0],argv[1],argv[2],argv[3]); */
  if ( (argc < 2 ) || 
       (argc == 2 && ( strcmp(argv[1],"init") != 0 && 
		       strcmp(argv[1],"finish") !=0 &&
		       strcmp(argv[1],"grid") != 0 &&
		       strcmp(argv[1],"minvalue") != 0 &&
		       strcmp(argv[1],"maxvalue") !=0 &&
		       strcmp(argv[1],"smoothing") != 0 &&
		       strcmp(argv[1],"polygonise") != 0 ) ) ||
       (argc == 3 && ( strcmp(argv[1],"end") != 0 &&
		       strcmp(argv[1],"isolevel") != 0 &&
		       strcmp(argv[1],"interpolate") != 0 &&
		       strcmp(argv[1],"get") != 0 &&
		       strcmp(argv[1],"smoothsteps") != 0 &&
		       strcmp(argv[1],"smoothweight") !=0 ) ) ||
       (argc == 4 && ( strcmp(argv[1],"isolevel") != 0 &&
		       strcmp(argv[1],"save") != 0 &&
		       strcmp(argv[1],"get") != 0 &&
		       strcmp(argv[1],"polygonise") != 0 ) ) ||
       (argc == 6 && ( strcmp(argv[1],"isoplane") != 0 &&
		       strcmp(argv[1],"isoplaneconfig") !=0 &&
		       strcmp(argv[1],"polygonise") != 0 ) ) ||
       (argc == 7 && ( strcmp(argv[1],"isoplane") != 0 &&
		       strcmp(argv[1],"isoplaneconfig") != 0 ) ) ||
       (argc > 7 && (strcmp(argv[1],"isoplaneconfig") != 0 && strcmp(argv[1],"polygonise") != 0) ) ||
       (argc != 2 && argc != 4 && argc != 6 && argc != 8 && strcmp(argv[1],"polygonise") == 0 ) )       
    {
      fprintf(stderr,"debug> argc=%d. argv[1]=%s\n", argc, argv[1]);
      fprintf(stderr,"debug> expr = %d\n", (argc != 2 && argc != 4 && argc != 6 && strcmp(argv[1],"polygonise")));
      Tcl_SetResult(interp, "Usage: xc_iso init|finish|(end isostack|isosign)|grid|save filename identifier|minvalue|maxvalue|(isolevel <value> ?<value>?)|polygonise ...|(smoothsteps <nstep>)|smoothweight <weight>)|(get smoothsteps|smoothweight)|(isoplane <type> <colorBase> <scalefunc> <what> ?<islide>?)|isoplaneconfig <object> -isoplanemin  <min>\n-isoplanemax  <max>\n-isolinecolor {monocolor r g b}\n-isolinewidth width\n-isoplanelighting 0/1| polycolor\nisolinedash  negdash | fulldash\n-isolinenlevels <nlevels>\n", TCL_STATIC);
      return TCL_ERROR;
    }

  /* ------------ XC_ISO INIT ---------------------------------------------- */
  if ( strcmp(argv[1],"init") == 0 ) {
    isostate.stateflag = ISO_INIT; /* here comes initialization */
    InitIsoDataArr();
    InitIsoExpand();
    /* malloc ***plvertex */
    if (!plvertex)
      plvertex = (PLANEVERTEX ***) malloc( (size_t) MAX_ISOOBJECTS *
					   sizeof(PLANEVERTEX **) );
    if (!plvertex) xcError("allocation failure for plvertex");
  } 

  /* ------------ XC_ISO FINISH -------------------------------------------- */
  else if ( strcmp(argv[1],"finish") == 0 ) {
    /* xcIsoError() will do precisely what we want */
    xcIsoError();
    isostate.stateflag = ISO_NULL;
    if ( dimType==XC_3D ) {
      if (is.stickmode && !is.ballmode) xcMakeProjection3D("sticks");
      if (is.ballmode) xcMakeProjection3D("balls");
      if (is.spacefillmode) xcMakeProjection3D("space");
    }
  }

  /* ------------ XC_ISO GRID -------------------------------------------- */
  else if ( strcmp(argv[1],"grid") == 0 ) {
    char *grid_string = (char *) Tcl_Alloc( sizeof(char) * 128);
    sprintf(grid_string, "%d %d %d", grd.nx, grd.ny, grd.nz);
    Tcl_SetResult(interp, grid_string, TCL_DYNAMIC);
  }

  /* ------------ XC_ISO SAVE ---------------------------------------------- */
  else if ( strcmp(argv[1],"save") == 0 ) {
    FILE *fp;
    /* save the calculated grid of points */
    if ( !(isostate.stateflag & 
	   (ISO_STACK | ISO_SIGN | ISO_FILES | ISO_POINTS | ISO_DATA)) )  {
      Tcl_SetResult(interp, "\"xc_isostack\", \"xc_isosign\", \"xc_isofiles\", \"xc_isopoints\" and \"xc_isodata\" should be called before \"xc_iso minvalue\" command", TCL_STATIC);
      return TCL_ERROR;
    }
    if ( (fp = fopen(argv[2], "w")) == NULL ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "couldn't open %s file for writing", argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return XC_ERROR;
    }
    if ( isodata.dim[ISOOBJ_BASE] == 2 ) 
      WriteDataGrid(fp, DATAGRID_2D, argv[3], ISOOBJ_BASE);
    if ( isodata.dim[ISOOBJ_BASE] == 3 )
      WriteDataGrid(fp, DATAGRID_3D, argv[3], ISOOBJ_BASE);      
    fclose(fp);
  }

  /* ------------ XC_ISO END ----------------------------------------------- */
  else if ( strcmp(argv[1],"end") == 0 ) {
    if ( strcmp(argv[2],"isostack") == 0 ) {
      /* check if all isodata.nframe[][1] are the same */
      int i;
      for (i=1; i<isodata.nframe[0][3]; i++) 
	if ( isodata.nframe[0][1] != isodata.nframe[i][1] ) {
	  Tcl_SetResult(interp, "stacks #1 were not specified correctly, because they don't have equal number of frames", TCL_STATIC);
	  return TCL_ERROR;
	}
      /* here is a place for assigning grd.nz */
      grd.nz    = isodata.nframe[0][1];
      newgrd.nz = grd.nz;
      if ( !(isostate.stateflag & ISO_STACK) ) isostate.stateflag |= ISO_STACK;
    }
    else if ( strcmp(argv[2],"isosign") == 0 ) {
      if ( !(isostate.stateflag & ISO_SIGN) ) isostate.stateflag |= ISO_SIGN;
    }
    else {
      Tcl_SetResult(interp, "Usage: xc_iso end isostack|isosign", TCL_STATIC);
      xcIsoError();
      return TCL_ERROR;
    }
  }

  /* ------------ XC_ISO MINVALUE ------------------------------------------ */
  else if ( strcmp(argv[1],"minvalue") == 0 ) {
    char *result = Tcl_Alloc( sizeof(char) * 128); /* maximum lenght of result 
						      string is 128 characters */
    double min_value = 99.9e+99;
    int ix, iy, iz;
    /* to execute "xc_iso minvalue" everything must be prepared */
    if ( !(isostate.stateflag & 
	   (ISO_STACK | ISO_SIGN | ISO_FILES | ISO_POINTS | ISO_DATA)) ) {
      Tcl_SetResult(interp, "\"xc_isostack\", \"xc_isosign\", \"xc_isofiles\", \"xc_isopoints\" and \"xc_isodata\" should be called before \"xc_iso minvalue\" command"     , TCL_STATIC);
      return TCL_ERROR;
    }

    /* BUG: below works just for 3D grids */
    if ( isodata.dim[ISOOBJ_BASE] == 3 ) {
      for (ix = 0; ix < newgrd.nx; ix++) 
	for (iy = 0; iy < newgrd.ny; iy++) 
	  for (iz = 0; iz < newgrd.nz; iz++)
	    if ( gridvertex[ix][iy][iz].val < min_value ) {
	      min_value = gridvertex[ix][iy][iz].val;
	    }
    } else if ( isodata.dim[ISOOBJ_BASE] == 2 ) {
      for (ix = 0; ix < newgrd.nx; ix++) 
	for (iy = 0; iy < newgrd.ny; iy++) 
	  if ( plvertex[ISOOBJ_BASE][ix][iy].val < min_value ) {
	    min_value = plvertex[ISOOBJ_BASE][ix][iy].val;
	  }
    }


    sprintf(result,"%f", min_value);
    Tcl_SetResult(interp, result, TCL_DYNAMIC);
  }

  /* ------------ XC_ISO MAXVALUE------------------------------------------- */
  else if ( strcmp(argv[1],"maxvalue") == 0 ) {
    char *result = Tcl_Alloc( sizeof(char) * 128); /* maximum lenght of result 
						      string is 128 characters */
    double max_value = -99.9e+99;
    int ix, iy, iz;
    /* to execute "xc_iso minvalue" everything must be prepared */
    if ( !(isostate.stateflag & 
	   (ISO_STACK | ISO_SIGN | ISO_FILES | ISO_POINTS | ISO_DATA)) ) {
      Tcl_SetResult(interp, "\"xc_isostack\", \"xc_isosign\", \"xc_isofiles\", \"xc_isopoints\" and \"xc_isodata\" should be called before \"xc_iso maxvalue\" command", TCL_STATIC);
      return TCL_ERROR;
    }
    /* BUG: below works just for 3D grids */
    if ( isodata.dim[ISOOBJ_BASE] == 3 ) {
      for (ix = 0; ix < newgrd.nx; ix++) 
	for (iy = 0; iy < newgrd.ny; iy++) 
	  for (iz = 0; iz < newgrd.nz; iz++)
	    if ( gridvertex[ix][iy][iz].val > max_value ) {
	      max_value = gridvertex[ix][iy][iz].val;
	    }
    } else if ( isodata.dim[ISOOBJ_BASE] == 2 ) {
      for (ix = 0; ix < newgrd.nx; ix++) 
	for (iy = 0; iy < newgrd.ny; iy++) 
	  if ( plvertex[ISOOBJ_BASE][ix][iy].val > max_value ) {
	    max_value = plvertex[ISOOBJ_BASE][ix][iy].val;
	  }
    }

    sprintf(result, "%f", max_value);
    Tcl_SetResult(interp, result, TCL_DYNAMIC);
  } 

  /* ------------ XC_ISO GET------------------------------------------- */
  else if ( strcmp(argv[1],"get") == 0 ) {
    char *result = Tcl_Alloc( sizeof(char) * 128); /* maximum lenght of result 
						      string is 128 characters */
    if ( strncmp(argv[2], "smooths", 7) == 0 ) 
      /* -- SMOOTHSTEPS -- */
      sprintf(result,"%d", isodata.smoothsteps);
    else if ( strncmp(argv[2], "smoothw", 7) == 0 )
      /* -- SMOOTHWEIGHT -- */
      sprintf(result,"%f", isodata.smoothweight);
    else {
      Tcl_SetResult(interp, "Usage: xc_iso get smoothsteps|smoothweight", TCL_STATIC);
      return TCL_ERROR;
    }
    Tcl_SetResult(interp, result, TCL_DYNAMIC);
  }

  /* ------------ XC_ISO ISOLEVEL ------------------------------------------ */
  else if ( strcmp(argv[1],"isolevel") == 0 ) {
    int i, ii;
    double isolevel;
    if ( argc < 3 && argc > 4 ) {
      Tcl_SetResult(interp, "Ugase: xc_iso isolevel <value> ?<value>?", TCL_STATIC);
      return TCL_ERROR;
    }

    isodata.nlevel = argc - 2;

    for (i=2; i<argc; i++) {
      ii = i - 2;
      if ( Tcl_GetDouble(interp, argv[i], &isolevel) == TCL_ERROR ) {
	if ( argc == 3 ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s %s", argv[2], argv[0], argv[1], argv[2]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	} else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s %s %s", argv[2], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }

      isodata.isolevel[ii] = (float) isolevel;
    }

    /* maybe isolevel is out of range */      
    if ( argc == 3 && ( isodata.isolevel[ii] < isodata.min || 
                        isodata.isolevel[ii] > isodata.max ) ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"specified isolevel %f is out of range, should be between [%f - %f]", isodata.isolevel[ii], isodata.min, isodata.max);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    else if ( argc == 4 ) {
      /*
        xc_iso isolevel +value -value
      
        at least one value has to be in range !!
      */
      int l0 = ( isodata.isolevel[0] < isodata.min || isodata.isolevel[0] > isodata.max );
      int l1 = ( isodata.isolevel[1] < isodata.min || isodata.isolevel[1] > isodata.max );
      if ( l0 && l1 ) {        
	char rss[1024];
	snprintf(rss, sizeof(rss),"specified isolevels %f,%f are out of range, should be between [%f - %f]",
                 isodata.isolevel[0], isodata.isolevel[1], isodata.min, isodata.max);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }

    if ( !(isostate.stateflag & ISO_ISOLEVEL) ) 
      isostate.stateflag |= ISO_ISOLEVEL;    
  }

  /* ------------ XC_ISO POLYGONISE ---------------------------------------- */
  else if ( strcmp(argv[1],"polygonise") == 0 ) {
    int i;
    int algorithm     = ISOSURF_TETRAHEDRAL;
    int shade_model   = ISOSURF_SHADE_SMOOTH;
    int normals_model = ISOSURF_NORMALS_GRADIENT;

    for (i=2; i<argc; i+=2)
      {
	/* argv[i] must be -algorithm || -shademodel */
	if ( strncmp(argv[i],"-algo",5) == 0 ) 
	  {
	    if ( strncmp(argv[i+1], "cube", 4) == 0 ) 
	      { 
		algorithm = ISOSURF_MARCHING_CUBES; 
	      }
	    else if ( strncmp(argv[i+1], "tetr", 4) == 0 ) 
	      { 
		algorithm = ISOSURF_TETRAHEDRAL; 
	      }
	    else 
	      {
		char rss[1024];
		snprintf(rss, sizeof(rss),"wrong usega of xc_iso polygonise, must be: polygonise ?-algorithm cubes|tetrahedrons?\n");
		Tcl_SetResult(interp, rss, TCL_VOLATILE);
		return TCL_ERROR;
	      }
	  }
	else if ( strncmp(argv[i],"-shade",6) == 0 ) 
	  {
	    if ( strcmp(argv[i+1],"smooth") == 0)     shade_model = ISOSURF_SHADE_SMOOTH;
	    else if ( strcmp(argv[i+1],"flat") == 0 ) shade_model = ISOSURF_SHADE_FLAT;
	    else {
	      char rss[1024];
	      snprintf(rss, sizeof(rss),"unknown value \"%s\" for -shademodel option, must be either \"smooth\" or \"flat\", while executing %s %s",
		       argv[i+1], argv[0], argv[1]);
	      Tcl_SetResult(interp, rss, TCL_VOLATILE);
	      return TCL_ERROR;
	    }
	  }
	else if ( strncmp(argv[i], "-norm", 5) == 0)
	  {
	    if ( strncmp(argv[i+1],"grad",4) == 0)         normals_model = ISOSURF_NORMALS_GRADIENT;
	    else if ( strncmp(argv[i+1],"triang",6) == 0 ) normals_model = ISOSURF_NORMALS_TRIANGLE;
	    else {
	      char rss[1024];
	      snprintf(rss, sizeof(rss),"unknown value \"%s\" for -normals option, must be either \"gradient\" or \"trianglet\", while executing %s %s", argv[i+1], argv[0], argv[1]);
	      Tcl_SetResult(interp, rss, TCL_VOLATILE);
	      return TCL_ERROR;
	    }
	  }
      }

    /* to execute "xc_iso polygonise" everything must be prepared */
    if ( !(isostate.stateflag & 
	   (ISO_STACK | ISO_SIGN | ISO_FILES | ISO_POINTS | 
	    ISO_DATA | ISO_ISOLEVEL)) ) {
      Tcl_SetResult(interp, "\"xc_isostack\", \"xc_isosign\", \"xc_isofiles\", \"xc_isopoints\", \"xc_isodata\" and \xc_iso isolevel <level>\" should be called before \"xc_iso polygonise\" command", TCL_STATIC);
      return TCL_ERROR;
    }
    /* for polygonisation isodata.dim[ISOOBJ_BASE] must be 3 */
    if ( isodata.dim[ISOOBJ_BASE] != 3 ) {
      Tcl_SetResult(interp, "\"xc_isopoints <object> 3D\" must be set, before calling \"xc_iso polygonise\"", TCL_STATIC);
      return TCL_ERROR;
    }

    for (i=0; i<isodata.nlevel; i++) 
      {	
	/* polygonise */
	MarchingCubes(isodata.isolevel[i], i, isodata.smoothsteps,
                      isodata.smoothweight, algorithm, shade_model, normals_model);
      }

    if ( !(isostate.stateflag & ISO_POLYGONISE) ) 
      isostate.stateflag |= ISO_POLYGONISE;
  }

  /* ------------ XC_ISO SMOOTHSTEP ---------------------------------------- */
  else if ( strcmp(argv[1],"smoothsteps") == 0 ) {
    /* it just remember the # of smoothing steps */
    if ( Tcl_GetInt(interp, argv[2], &isodata.smoothsteps) == TCL_ERROR ) {
      isodata.smoothsteps=0;
      Tcl_SetResult(interp, "read character instead of integer, while executing \"xc_iso smoothsteps <nstep>\" command", TCL_STATIC);
      return TCL_ERROR;
    }
  }

  /* ------------ XC_ISO SMOOTHWEIGHT--------------------------------------- */
  else if ( strcmp(argv[1],"smoothweight") == 0 ) {
    double weight;
    /* it just remember the # of smoothing steps */
    if ( Tcl_GetDouble(interp, argv[2], &weight) == TCL_ERROR ) {
      isodata.smoothweight=0.5;
      Tcl_SetResult(interp, "read character instead of double, while executing \"xc_iso smoothweight <weight>\" command", TCL_STATIC);
      return TCL_ERROR;
    }      
    isodata.smoothweight = (float) weight;
  }

  /* ------------ XC_ISO SMOOTHING --------------------------------------- */
  else if ( strcmp(argv[1],"smoothing") == 0 ) {
    int i;
    /* it is not OK as that OK as that */
    /* to execute "xc_iso poligonise" everything must be prepared */
    if ( !(isostate.stateflag & 
	   (ISO_STACK | ISO_SIGN | ISO_FILES | ISO_POINTS | 
	    ISO_DATA | ISO_ISOLEVEL | ISO_POLYGONISE)) ) {
      Tcl_SetResult(interp, "\"xc_isostack\", \"xc_isosign\", \"xc_isofiles\", \"xc_isopoints\", \"xc_isodata\" \"xc_iso isolevel <level>\" and \"xc_iso polygonise\" should be called before \"xc_iso smoothing\" command", TCL_STATIC);
      return TCL_ERROR;
    }
    /* for polygonisation isodata.dim[ISOOBJ_BASE] must be 3 */
    if ( isodata.dim[ISOOBJ_BASE] != 3 ) {
      Tcl_SetResult(interp, "\"xc_isopoints <object> 3D\" must be set, before calling \"xc_iso smoothing\"", TCL_STATIC);
      return TCL_ERROR;
    }

    for (i=0; i<isodata.nlevel; i++)
      /* insert here */
      WrapperSurfSmoothing(i);
  }

  /* ------------ XC_ISO ISOPLANE ---------------------------------------- */
  else if ( strcmp(argv[1], "isoplane") == 0 ) {
    int type; 
    int cb = COLORBASE_MONO;
    int fn = SCALE_FUNC_LIN;

    if ( argc < 6 || argc > 7) {
      Tcl_SetResult(interp, "Usage: xc_iso isoplane <type> <colorBase> <scalefunc> <what> ?<islide>?", TCL_STATIC);
      return TCL_ERROR;
    }
    /*
     * <what> can be: none|colorplane|isoline|both
     */

    /* to execute "xc_iso isoplane" everything must be prepared */
    if ( !(isostate.stateflag & 
	   (ISO_STACK | ISO_SIGN | ISO_FILES | ISO_POINTS | ISO_DATA )) ) {
      Tcl_SetResult(interp, "\"xc_isostack\", \"xc_isosign\", \"xc_isofiles\", \"xc_isopoints\" and \"xc_isodata\" should be called before \"xc_iso isoplane\" command", TCL_STATIC);
      return TCL_ERROR;
    }

    /*
     * get the TYPE of isoplane
     */
    if ( Tcl_GetInt(interp, argv[2], &type ) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" while executing %s %s %s %s %s", argv[2], argv[0], argv[1], argv[2], argv[3], argv[4]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( type < ISOOBJ_BASE || type >= MAX_ISOOBJECTS ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"object type \"%d\" for xc_iso isoplane is out of range, must be within [%d,%d]", type, ISOOBJ_BASE,MAX_ISOOBJECTS-1);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    /* if type == ISOOBJ_BASE --> 
       we deal with isoplane as base iso-object, so:
       isodata.dim[ISOOBJ_BASE] must be 2

       if type >  ISOOBJ_PLANE1|2|3 --> 
       we want to draw isoplane as a secondary isoobject, so:
       isodata.dim[ISOOBJ_BASE] must be 3
    */
    if ( type == ISOOBJ_BASE && isodata.dim[ISOOBJ_BASE] != 2 ) {
      Tcl_SetResult(interp, "\"xc_isopoints 0 2D\" must be set, before calling \"xc_iso isoplane 0\"", TCL_STATIC);
      return TCL_ERROR;
    } else if ( IsIsoPlane123( type ) && isodata.dim[ISOOBJ_BASE] != 3 ) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "\"xc_isopoints 0 3D\" must be set, before calling \"xc_iso isoplane %d\"", type);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    /*
     * get colorbase
     */
    if ( Tcl_GetInt(interp, argv[3], &cb ) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" while executing %s %s %s %s %s", argv[3], argv[0], argv[1], argv[2], argv[3], argv[4]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( cb < COLORBASE_FIRST || cb > COLORBASE_LAST ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"colorbase %d is out of range, must be within [%d,%d]", cb, COLORBASE_FIRST, COLORBASE_LAST);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    /*
     * get scalefunction
     */
    if ( Tcl_GetInt(interp, argv[4], &fn ) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" while executing %s %s %s %s %s", argv[4], argv[0], argv[1], argv[2], argv[3], argv[4]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( fn < SCALE_FUNC_FIRST || fn > SCALE_FUNC_LAST ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"scale-function %d is out of range, must be within [%d,%d]", fn, SCALE_FUNC_FIRST, SCALE_FUNC_LAST);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    /*
     * get what
     */
    if ( strcmp(argv[5], "none") == 0 ) {
      int flag = ISO_NULL;
      if ( type == ISOOBJ_BASE ) {
	flag |= ISO_COLORPLANE | ISO_ISOLINE;
	VPf.colorplane[ISOOBJ_BASE] = GL_FALSE;
	VPf.isoline[ISOOBJ_BASE]    = GL_FALSE;
      } else if ( type == ISOOBJ_PLANE1 ) {
	flag |= ISO_COLORPLANE1 | ISO_ISOLINE1;
	VPf.colorplane[ISOOBJ_PLANE1] = GL_FALSE;
	VPf.isoline[ISOOBJ_PLANE1]    = GL_FALSE;
      }
      else if ( type == ISOOBJ_PLANE2 ) {
	flag |= ISO_COLORPLANE2 | ISO_ISOLINE2;
	VPf.colorplane[ISOOBJ_PLANE2] = GL_FALSE;
	VPf.isoline[ISOOBJ_PLANE2]    = GL_FALSE;
      }
      else if ( type == ISOOBJ_PLANE3 ) {
	flag |= ISO_COLORPLANE3 | ISO_ISOLINE3;
	VPf.colorplane[ISOOBJ_PLANE3] = GL_FALSE;
	VPf.isoline[ISOOBJ_PLANE3]    = GL_FALSE;
      }
      xcDeleteBitFlags( &isostate.stateflag, flag ); 
    }
    else if ( strcmp(argv[5], "colorplane") == 0 || 
	      strcmp(argv[5], "isoline") == 0 ||
	      strcmp(argv[5], "both") == 0 ) {
      int i;
      int islide;
      int nslide = newgrd.nz;
      int a      = newgrd.nx;
      int b      = newgrd.ny; /* nslide, a, b will be reassigned if PLANE123 */

      /*
       * do we have ISOOBJ_PLANE1|2|3
       */
      if ( IsIsoPlane123( type ) ) {
	/*
	  PLANE1 == vec0 x vec1;
	  PLANE2 == vec2 x vec0;
	  PLANE3 == vec1 x vec2;
	*/
	if ( type == ISOOBJ_PLANE2 ) {
	  a      = newgrd.nz;
	  b      = newgrd.nx;
	  nslide = newgrd.ny;
	}
	if ( type == ISOOBJ_PLANE3 ) {
	  a      = newgrd.ny;
	  b      = newgrd.nz;
	  nslide = newgrd.nx;
	}
	/* first get the number of slide; islide = [1..newgrd.nx|ny|nz] */
	if ( Tcl_GetInt(interp, argv[6], &islide ) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" while executing %s %s %s %s %s ...", argv[6], argv[0], argv[1], argv[2], argv[3], argv[4]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	/* 
	 * check if islide is within [1..nslide]
	 */
	if ( islide < 1 || islide > nslide ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"islide #%d is out of range, should be within [%d,%d]", islide, 1, nslide);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}

	/*
	 * malloc plvertex[type]
	 */
	if ( !isostate.plvertex_malloc[type] ) {
	  plvertex[type] = xcMallocPLANEVERTEX( a, b );
	  isostate.plvertex_malloc[type] = 1;
	}
	if ( (strcmp( argv[5], "isoline") == 0 ||
	      strcmp( argv[5], "both") == 0) && 
	     !isostate.isoline2D_malloc[type] ) {
	  for(i=0; i<ISOLINE_MAXLEVEL; i++) {
	    isostate.max_n_isoline2D[type][i] = a;
	    isoline2D[type].segment[i] = 
	      xcMallocLINE( isostate.max_n_isoline2D[type][i] );
	  }
	  isostate.isoline2D_malloc[type] = 1;
	}    
	/*
	 * now read data from gridvertex
	 */
	ReadPlvertex123( type, islide );
      }

      if ( strcmp( argv[5], "colorplane" ) == 0 ||
	   strcmp(argv[5], "both") == 0) {
	/*
	 * update state flag
	 */
	if ( type == ISOOBJ_BASE && !(isostate.stateflag & ISO_COLORPLANE) )
	  isostate.stateflag |= ISO_COLORPLANE;
	if ( type == ISOOBJ_PLANE1 && !(isostate.stateflag & ISO_COLORPLANE1) )
	  isostate.stateflag |= ISO_COLORPLANE1;
	if ( type == ISOOBJ_PLANE2 && !(isostate.stateflag & ISO_COLORPLANE2) )
	  isostate.stateflag |= ISO_COLORPLANE2;
	if ( type == ISOOBJ_PLANE3 && !(isostate.stateflag & ISO_COLORPLANE3) )
	  isostate.stateflag |= ISO_COLORPLANE3;

	/* make a colorplane */
	if ( !isostate.plvertex_malloc[type] ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "bug in xcIsoSurf.c; plvertex[%d] was not malloced; before ColorPlane( type, cb, fn, a, b )", type);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	ColorPlane( type, cb, fn, a, b );
      } 
      if ( strncmp(argv[5], "isoline", 7) == 0 ||
	   strcmp(argv[5], "both") == 0 ) {
	/*
	 * update state flag
	 */
	if ( type == ISOOBJ_BASE && !(isostate.stateflag & ISO_ISOLINE) )
	  isostate.stateflag |= ISO_ISOLINE;
	if ( type == ISOOBJ_PLANE1 && !(isostate.stateflag & ISO_ISOLINE1) )
	  isostate.stateflag |= ISO_ISOLINE1;
	if ( type == ISOOBJ_PLANE2 && !(isostate.stateflag & ISO_ISOLINE2) )
	  isostate.stateflag |= ISO_ISOLINE2;
	if ( type == ISOOBJ_PLANE3 && !(isostate.stateflag & ISO_ISOLINE3) )
	  isostate.stateflag |= ISO_ISOLINE3;

	/* make isolines */
	if ( !isostate.isoline2D_malloc[type] ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "bug in xcIsoSurf.c; isoline2D[%d] was not malloced; before Isoline2D( ptype, cb, fn, a, b )", type);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	IsoLine2D( type, cb, fn, a, b );
      }    
    } 
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss), "unknown <what> \"%s\" for xc_iso isoplane, must be one of colorplane, isoline or both, while executing %s %s %s %s %s %s", argv[5], argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }   
  }
  /* ------------ XC_ISO INTERPOLATE --------------------------------------- */
  else if ( strcmp(argv[1],"interpolate") == 0 ) {
    int i, degree;
    if ( argc != 3 ) {
      Tcl_SetResult(interp, "Usage: xc_iso interpolate <degree>", TCL_STATIC);
      return TCL_ERROR;
    }
    /* to execute "xc_iso interpolate" everything must be prepared */
    if ( !(isostate.stateflag & 
	   (ISO_STACK | ISO_SIGN | ISO_FILES | ISO_POINTS | ISO_DATA )) ) {
      Tcl_SetResult(interp, "\"xc_isostack\", \"xc_isosign\", \"xc_isofiles\", \"xc_isopoints\" and \"xc_isodata\" should be called before \"xc_iso interpolate\" command", TCL_STATIC);
      return TCL_ERROR;
    }
    if ( Tcl_GetInt(interp, argv[2], &degree ) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" while executing %s %s %s", argv[3], argv[0], argv[1], argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    /*
     * so far we can interpolate just ISOOBJ_BASE, and this is needed
     */
    for (i=ISOOBJ_PLANE1; i<=ISOOBJ_PLANE3; i++) {
      if ( isostate.plvertex_malloc[i] ) {
	xcFree_PLANEVERTEX( plvertex[i] );
	isostate.plvertex_malloc[i] = 0;
      }
      /* 
	 we do not need to do this, seens we use xcReallocate2xBigger routines:

	 if ( isostate.isoline2D_malloc[i] ) {
	 xcFree_LINE( isoline2D[i].segment );
	 isostate.isoline2D_malloc[i] = 0;
	 }
      */
    }

    /* now we should delete all flags about PLANE123 */
    xcDeleteBitFlags( &isostate.stateflag, 
		      (ISO_COLORPLANE1 | ISO_COLORPLANE2 | ISO_COLORPLANE3 | 
		       ISO_ISOLINE1 | ISO_ISOLINE2 | ISO_ISOLINE3 ) );

    isoInterpolate( degree );
  }

  /* ------------ XC_ISO ISOPLANECONFIG------------------------------------- */
  else if ( strcmp(argv[1],"isoplaneconfig") == 0 ) {
    int i, j, obj;
    static int first_time[MAX_ISOOBJECTS] = {1, 1, 1, 1};
    /* Usage:
     *  xc_iso isoplaneconfig <object> \
     *                            -isoplanemin  <min> \
     *                            -isoplanemax  <max> \
     *                            -isolinecolor {monocolor r g b} | polycolor \
     *                            -isolinewidth width \
     *                            -isolinedash  nodash | negdash | fulldash \
     *                            -isolinenlevels <nlevels> \
     *                            -isoplanelighting   0/1)
     *
     * Defaults: -isolinecolor   == monocolor 1.0, 1.0, 1.0, 1.0
     *           -isolinedash    == negdash
     *           -isolinenlevels == 15
     *           -isolinewidth   == 2
     */
    if ( argc % 2 != 1) {
      Tcl_SetResult(interp, "Usage: xc_iso isoplaneconfig <object> -isoplanemin  <min>\n-isoplanemax  <max>\n-isolinecolor monocolor\n-isolinewidth width\n-isoplanelighting 0/1 | polycolor\n-isolinedash  negdash | fulldash\n-isolinenlevels <nlevels>)", TCL_STATIC);
      return TCL_ERROR;
    }
    /* get object */
    if ( Tcl_GetInt(interp, argv[2], &obj ) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"expected integer, but got %s", argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( obj < ISOOBJ_BASE || obj >= MAX_ISOOBJECTS ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"<object> type \"%d\" is out of range, should be between [%d,%d]", obj, ISOOBJ_BASE, MAX_ISOOBJECTS-1);      
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    /* set DEFAULTS */
    if (first_time[obj]) {
      isoline2D[obj].linecolor = ISOLINE_MONOCOLOR;
      for (j=0; j<4; j++)
	isoline2D[obj].monocolor[j] = 1.0;
      isoline2D[obj].linedash = ISOLINE_NEGDASH;
      isoline2D[obj].nlevel = 15;
      first_time[obj] = 0;
    }
    for (i=3; i<argc; i+=2) {
      if ( strcmp(argv[i], "-isolinecolor") == 0 ) {
	if ( strcmp(argv[i+1], "polycolor") == 0 ) {
	  /* colors of isolines according to color base */
	  isoline2D[obj].linecolor = ISOLINE_POLYCOLOR;
	} else {
	  /* {monocolor r, g, b ?a?} */
	  int   argcList;
	  double c[4];
	  const char  **argvList;
	  Tcl_SplitList(interp, argv[i+1], &argcList, &argvList);
	  /* monocolor */
	  if ( strcmp(argvList[0], "monocolor") != 0 || 
	       argcList < 3 || argcList > 4 ) {
	    char rss[1024];
	    snprintf(rss, sizeof(rss), "expected \"monocolor r, g, b, ?a?\", but got \"%s\"", argv[i+1]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    Tcl_Free((char *) argvList);
	    return TCL_ERROR;
	  }
	  for (j=1; j<argcList; j++)
	    /* r,g,b */
	    if ( Tcl_GetDouble(interp, argvList[j], &c[j-1] ) == TCL_ERROR ) {
	      char rss[1024];
	      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" while executing xc_iso isoplane configure -isolinecolor {monocolor r, g, b, ?a?}", 
		       argvList[j]);
	      Tcl_SetResult(interp, rss, TCL_VOLATILE);
	      Tcl_Free((char *) argvList);
	      return TCL_ERROR;
	    }
	  /* now assing monocolor */
	  isoline2D[obj].monocolor[0] = (float) c[0];
	  isoline2D[obj].monocolor[1] = (float) c[1];
	  isoline2D[obj].monocolor[2] = (float) c[2];
	  isoline2D[obj].monocolor[3] = 1.0;
	  if ( argcList == 5 ) isoline2D[obj].monocolor[3] = (float) c[3];
	  /* colors of isolines according to color base */
	  isoline2D[obj].linecolor = ISOLINE_MONOCOLOR;
	  Tcl_Free((char *) argvList);
	}
      }
      else if ( strcmp(argv[i], "-isolinewidth") == 0 ) {
	double width;
	if ( Tcl_GetDouble(interp, argv[i+1], &width) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s ... %s %s ...", argv[i+1], argv[0], argv[1], argv[i], argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	isoline2D[obj].linewidth = (float) width;
      }
      else if ( strcmp(argv[i], "-isolinedash") == 0 ) {
	if ( strcmp(argv[i+1], "nodash") == 0 ) {
	  isoline2D[obj].linedash = ISOLINE_NODASH;	  
	} else if ( strcmp(argv[i+1], "negdash") == 0 ) {
	  isoline2D[obj].linedash = ISOLINE_NEGDASH;
	} else if ( strcmp(argv[i+1], "fulldash") == 0 ) {
	  isoline2D[obj].linedash = ISOLINE_FULLDASH;
	} else {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "unknown isolinedash %s, should be negdash or fulldash", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }
      else if ( strcmp(argv[i], "-isolinenlevels") == 0 ) {
	int nlevel;
	if ( Tcl_GetInt( interp, argv[i+1], &nlevel ) == TCL_ERROR ) {
	  return TCL_ERROR;
	}
	isoline2D[obj].nlevel = nlevel;
      }
      else if ( strcmp(argv[i], "-isoplanemin") == 0 ) {
	double min;
	if ( Tcl_GetDouble( interp, argv[i+1], &min ) == TCL_ERROR ) {
	  return TCL_ERROR;
	}
	isodata.min_allowed[obj] = min;
      }
      else if ( strcmp(argv[i], "-isoplanemax") == 0 ) {
	double max;
	if ( Tcl_GetDouble( interp, argv[i+1], &max ) == TCL_ERROR ) {
	  return TCL_ERROR;
	}
	isodata.max_allowed[obj] = max;
      }
      else if ( strcmp(argv[i], "-isoplanelighting") == 0 ) {
	int lighting;
	if ( Tcl_GetInt( interp, argv[i+1], &lighting ) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "wanted integer, but got %s", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}	
	isoplaneDisp[obj].lighting = lighting;
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss), "unknown option \"%s\", must be one of -isolinecolor, -isolinedash, -isolinenlevels, -isoplanemin, -isoplanemax", argv[i]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
  }

  return TCL_OK;
}


/******************************************************************************
 * XC_IsostackCmd --> inplementation of 'xc_isostack' custom Tcl command 
 * --------------------- 
 * Usage:
 *       FIRST  specify:  xc_isostack <nstack> 
 *       if (nstack < 4) {
 *          specify:  xc_isostack 0 <stack=[2-0]> <n_of_frames_in_stack[2-0]>
 *       } else {
 *          specify:  xc_isostack <frame_in_stack3> <stack=[3-0]> 
 *                                <n_of_frames_in_stack[3,0]>
 * 
 *          SPECIAL CASE of this command is:
 *   xc_isostack 0 3 <n_of_frames_in_stack3> .... this specify number of frames
 *                                                in stack3
 *       }
 *
 * or
 *
 *       xc_isostack clean    --- clean stacks
 *
 * ----------------------------------------------------------------------------
 *            stack order:
 *                         order 0   ..... open-shell  case: nframe(0) = 2
 *                                         close-shell case: nframe(0) = 1
 *
 *                         order 1   ..... if 3D case, then this is number of
 *                                         points in Z directions
 *
 *                         order 2   ..... used for differentiall maps,
 *                                         for example: stack1 = PSCF,
 *                                                      stack2 = PATO
 *
 *                         order 3   ..... several files may be merged 
 *                                         together (example: comparison of
 *                                         methods); nframes in order 3 is
 *                                         number of files merged toghether;
 *                                         Each order-3 stack can have a 
 *                                         different number of frames in stack
 *                                         order 2 & order 0, but must have 
 *                                         a matching stacks of order 1 
 */
int
XC_IsostackCmd(ClientData clientData, Tcl_Interp *interp,
	       int argc, const char *argv[])
{ 
  int i, stack, nframe;
  static int frame3;
  static int is_nstack = 0;

  if ( argc != 2 && argc != 4 ) {
    Tcl_SetResult(interp, "Usage: xc_iso isoplane <type> <colorBase> <scalefunc> <what> ?<islide>?", TCL_STATIC);
    return TCL_ERROR;
  } 

  /* first consider CASE: 
   *                xc_isostack <n_of_stacks>
   */
  if ( argc == 2 && strcmp(argv[1],"clean") == 0 ) {
    InitIsoDataArr();
    is_nstack = 0;
  }
  else if ( argc == 2 ) {
    if ( Tcl_GetInt(interp, argv[1], &(isodata.nstack)) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" as first argument in \"xc_isostack\" command", argv[1]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    if ( isodata.nstack > ISODATA_MAXSTACK ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"to many stacks in \"xc_isostack\" command, maximum number of stacks is %d", ISODATA_MAXSTACK);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      isodata.nstack = 0;
      return TCL_ERROR;
    }
    is_nstack = 1;
    /* if isodata.nstack < 4, that means we have just one frame of order 3 */
    if ( isodata.nstack < 4 ) isodata.nframe[0][3] = 1;
  } else if ( argc == 4 ) {
    /* check if xc_isostack <nstack> was specified */
    if ( !is_nstack ) {
      Tcl_SetResult(interp,"before using xc_isostack <frame_in_stack3> <stack> <n_of_frames_in_stack>, xc_isostack <nstack> must be specified", TCL_STATIC);
      return TCL_ERROR;
    }
    /* CASE: xc_isostack <frame_in_stack3> <stack> <n_of_frames_in_stack> */
    if ( Tcl_GetInt(interp, argv[1], &frame3) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" when executing %s %s %s", argv[0], argv[1], argv[2], argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }    
    /* if isodata.nstack < 4 ==> if frame3 > 0 ==> ERROR */
    if ( isodata.nstack < 4 && frame3 > 0 ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"because only %d stacks are defined, only 0 can be specified as first argument to \"xc_isostack\" command", isodata.nstack);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    if ( Tcl_GetInt(interp, argv[2], &stack) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" as second argument to \"xc_isostack\" command", argv[2]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    /* stack can only be in range [0-2] for isodata.nstack < 4 and
     * in range [0-3] for isodata.nstack = 4 
     */
    if ( stack < 0 || ( (isodata.nstack <  4 && stack > 2) ||
			(isodata.nstack >= 4 && stack > 3) ) ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"stack too large - %d, number must not be greater than %d; while executing %s %s %s %s", 
	       stack, isodata.nstack-1, argv[0], argv[1], argv[2], argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }

    if ( Tcl_GetInt(interp, argv[3], &nframe) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" as the third argument to \"xc_isostack\" command", argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
    /* is stack 3 defined */
    if ( isodata.nstack < 4 ) {
      isodata.nframe[0][stack] = nframe;
    } else {
      if ( stack == 3 ) 
	for (i=0; i<nframe; i++)
	  isodata.nframe[i][3] = nframe;
      /* maybe we have stack-3 defined, but number of frames in stack 3 wasn't
       * yet specified */
      if (isodata.nframe[0][3] == 0) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"please specify number of frames in stack-3 first, before executing %s %s %s %s", argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }

      isodata.nframe[frame3][stack] = nframe;
    }
  }

  /* everything is OK */
  return TCL_OK;
}


/******************************************************************************
 * XC_IsosignCmd --> inplementation of 'xc_isosign' custom Tcl command
 * --------------------- 
 * Usage:
 *       xc_isosign <stack_level> <frame> <frame_sign>
 *   
 *            frame      .... frame's number or 'all' 
 *                             
 *            frame_sign ....-1|0|+1 => if -1 then data are multiplied by -1
 *                                      if 0 by 0 , if 1 by 1
 * DEFAULT sign is 0
 */
int
XC_IsosignCmd(ClientData clientData, Tcl_Interp *interp,
	      int argc, const char *argv[])
{
  int stacklevel, framenum, sign;

  if (argc != 4) {
    Tcl_SetResult(interp,"Usage: xc_isostack <stack_order> <frame> <frame's_sign>", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( !(isostate.stateflag & ISO_STACK) ) {
    Tcl_SetResult(interp, "\"xc_isostack\" command must be called before calling \"xc_isosign\" command", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( Tcl_GetInt(interp, argv[1], &(stacklevel)) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing \"%s %s %s\"", argv[0],argv[1],argv[2],argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  /* stacklevel shouldn't be greater than MAXSTACK-1 */
  if ( stacklevel > ISODATA_MAXSTACK - 1 ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"stack of order %d was specified, but maximim order of stack allowed is %d; while executing \"%s %s %s %s\"", 
	     stacklevel, ISODATA_MAXSTACK-1, argv[0],argv[1],argv[2],argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( strcmp(argv[2],"all") == 0 ) {
    ;
  } else {
    if ( Tcl_GetInt(interp, argv[2], &(framenum)) == TCL_ERROR ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing \"%s %s %s\"", argv[0],argv[1],argv[2],argv[3]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  if ( Tcl_GetInt(interp, argv[3], &(sign)) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing \"%s %s %s\"", argv[0],argv[1],argv[2],argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* sign must be -1 | 0 | 1 */
  if ( sign != -1 && sign != 0 && sign != 1 ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"sign must be one of -1, 0, 1; but got \"%d\", while executing \"%s %s %s %s\"", sign, argv[0],argv[1],argv[2],argv[3]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* UPDATING SIGN */
  if ( strcmp(argv[2],"all") != 0 ) {
    isodata.framesign[stacklevel][framenum] = sign;
  } else {
    int i;
    for (i=0; i<ISODATA_MAXFRAME; i++)
      isodata.framesign[stacklevel][i] = sign;
  }

  return TCL_OK;
}



/******************************************************************************
 * XC_IsofilesCmd --> inplementation of 'xc_isofiles' custom Tcl command
 * --------------------- 
 * Usage:
 *       xc_isofiles <binary_vertex_filename> <binary_filename> <filename0> 
 *                   ?<filename1>? ?<filename2>? ....
 * --------------------------------------------------------------------------
 * xc_isofiles require:
 *                      pre-calling of: xc_isostack
 *                                      xc_isosign
 */
int
XC_IsofilesCmd(ClientData clientData, Tcl_Interp *interp,
	       int argc, const char *argv[])
{ 
  int i;
  const char *version;

  if ( argc < 4 ) {
    Tcl_SetResult(interp, "Usage: xc_isofiles <binary_vertex_filename> <binary_filename> <filename0> ?<filename1>? ?<filename2>? ....", TCL_STATIC);
    return TCL_ERROR;
  } 

  if ( !(isostate.stateflag & (ISO_STACK | ISO_SIGN)) ) {
    Tcl_SetResult(interp, "\"xc_isostack\" & \"xc_isosign\" must be called before calling \"xc_isofiles\" command", TCL_STATIC);
    return TCL_ERROR;
  }

  /* open binary filenames */
  if (isostate.bin_vertex_file_open) {
    fclose( isodata.bin_vertex_fp );
    isostate.bin_vertex_file_open = 0;
  }
  if ( (isodata.bin_vertex_fp = fopen(argv[1],"w+")) == NULL ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"can't open file \"%s\"",argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    xcIsoError();
    return TCL_ERROR;
  }
  isostate.bin_vertex_file_open = 1;

  if ( isostate.bin_file_open ) {
    fclose(isodata.bin_fp);
    isostate.bin_file_open = 0;
  }
  if ( (isodata.bin_fp = fopen(argv[2],"w")) == NULL ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"can't open file \"%s\"",argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    xcIsoError();
    return TCL_ERROR;
  }
  isostate.bin_file_open        = 1;

  /*
    what is the CRYSTAL version we are using ???
    Load the value from system(c95_version) variable
  */
  fprintf(stderr, "c95_version-1\n"); fflush(stderr);
  version = Tcl_GetVar2(interp, "system", "c95_version", TCL_GLOBAL_ONLY);

  if ( strncmp("95",version,2) == 0 ) {
    crystal_version = 95;
  } 
  else if ( strncmp("98",version,2) == 0 ) {
    crystal_version = 98;
  }
  else if ( strncmp("03",version,2) == 0 ) {
    crystal_version = 3;
  } 
  else if ( strncmp("06",version,2) == 0 ) {
    crystal_version = 6;
  } 
  else if ( strncmp("09",version,2) == 0 ) {
    crystal_version = 9;
  } 
  else if ( strncmp("14",version,2) == 0 ) {
    crystal_version = 14;
  } 
  else {
    crystal_version = 0;
  }
  fprintf(stderr, "c95_version-2\n"); fflush(stderr);

  /* now checks if files really exist */
  for (i=3; i<argc; i++) {
    if ( (isodata.fp = fopen(argv[i],"r")) == NULL ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"can't open file \"%s\"",argv[i]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      xcIsoError();
      return TCL_ERROR;
    }
    if (!ReadBlock0(RB0_FIND, 0, 0, 0, i-3)) {
      char rss[1024];
      snprintf(rss, sizeof(rss), "%s", isodata.error);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      xcIsoError();
      return TCL_ERROR;
    }
    fclose(isodata.fp);
  }

  fclose(isodata.bin_fp);
  isodata.bin_fp = fopen(argv[2],"r"); /* don't need to put 
					* isostate.bin_file_open to 1, because
					* it is already 1 (look up 20 lines
					*/

  /* everything is OK */
  if ( !(isostate.stateflag & ISO_FILES) ) isostate.stateflag |= ISO_FILES;  
  return TCL_OK;
}



/******************************************************************************
 * XC_IsopointsCmd --> inplementation of 'xc_isopoints' custom Tcl command
 * --------------------- 
 * Usage:   xc_isopoints <object> 2D|3D <points>
 *
 *                       <object> .... ISOOBJ_BASE
 *                                .... ISOOBJ_PLANE1
 *                                .... ISOOBJ_PLANE2
 *                                .... ISOOBJ_PLANE3
 */
int
XC_IsopointsCmd(ClientData clientData, Tcl_Interp *interp,
		int argc, const char *argv[])
{
  int obj, nobj;
  int i, j, ij;
  double value;
  /*float scalar, Zvec[3] = { 1.0, 1.0, 10.0 };*/

  if ( argc < 3 ) {
    Tcl_SetResult(interp, "Usage: xc_isopoints 2D|3D <object> <points>", TCL_STATIC);
    return TCL_ERROR;
  } 

  /*
   * get the object type
   */
  if ( Tcl_GetInt(interp, argv[1], &obj) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" while executing %s %s %s ...", argv[1], argv[0], argv[1], argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( obj < 0 || obj >= MAX_ISOOBJECTS ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"object type \"%d\" is out of range", obj);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( strcmp(argv[2],"2D") == 0 ) { 
    isodata.dim[obj] = 2;
    if ( argc != 12 ) {
      Tcl_SetResult(interp, "invalid number of arguments in \"xc_isopoints <objectType> 2D\" command", TCL_STATIC);
      return TCL_ERROR;
    }
    /* delete the following bit flags */
    xcDeleteBitFlags( &isostate.stateflag, (ISO_DATA | ISO_POLYGONISE) );
    /* initialize IsoExpand struct */
    InitIsoExpand();
  } else if ( strcmp(argv[2],"3D") == 0 ) {
    if ( obj > ISOOBJ_BASE ) {
      Tcl_SetResult(interp, "only BASE isoobjects can be of 3D type", TCL_STATIC);
      return TCL_ERROR;
    }
    isodata.dim[obj] = 3;
    if ( argc != 15 ) {
      Tcl_SetResult(interp, "invalid number of arguments in \"xc_isopoints <objectType> 3D\" command", TCL_STATIC);
      return TCL_ERROR;
    }
    /* delete the following bit flags */
    xcDeleteBitFlags( &isostate.stateflag, (ISO_DATA | ISO_COLORPLANE) );
    /* initialize IsoExpand struct */
    InitIsoExpand();
  } else {
    char rss[1024];
    snprintf(rss, sizeof(rss),"unknown option \"%s\", must be one of 2D, 3D", 
	     argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  for (i=0; i<=isodata.dim[obj]; i++)
    for (j=0; j<3; j++) {
      ij = 3 + i * 3 + j;
      if ( Tcl_GetDouble(interp, argv[ij], &value) == 
	   TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted double, but got \"%s\" as %d-th argument while executing \"xc_isopoints\" command", argv[ij], ij);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      isodata.points[obj][i][j] = (float) value;
      if ( isodata.dim[obj] == 3 && obj == ISOOBJ_BASE ) {
	int ip;
	/* get also the points for PLANE1|2|3

	   PLANE1 == vec0 x vec1; -+
	   PLANE2 == vec2 x vec0;  |-> origin is point0
	   PLANE3 == vec1 x vec2; -+
	*/
	for(ip=ISOOBJ_PLANE1; ip<=ISOOBJ_PLANE3; ip++) {
	  isodata.points[ip][i][j] = (float) value;
	  isodata.dim[ip] = 2;
	}
      }
    }
  /* now made vectors; order of points submited is:
   *           origin, Xpoint, Ypoint, Zpoint 
   */
  isodata.vec[obj][0].x = isodata.points[obj][1][0] - isodata.points[obj][0][0];
  isodata.vec[obj][0].y = isodata.points[obj][1][1] - isodata.points[obj][0][1];
  isodata.vec[obj][0].z = isodata.points[obj][1][2] - isodata.points[obj][0][2];

  isodata.vec[obj][1].x = isodata.points[obj][2][0] - isodata.points[obj][0][0];
  isodata.vec[obj][1].y = isodata.points[obj][2][1] - isodata.points[obj][0][1];
  isodata.vec[obj][1].z = isodata.points[obj][2][2] - isodata.points[obj][0][2];

  if ( isodata.dim[obj] == 3 ) {
    isodata.vec[obj][2].x = isodata.points[obj][3][0] - isodata.points[obj][0][0];
    isodata.vec[obj][2].y = isodata.points[obj][3][1] - isodata.points[obj][0][1];
    isodata.vec[obj][2].z = isodata.points[obj][3][2] - isodata.points[obj][0][2];

    if ( obj == ISOOBJ_BASE ) {
      float det = XYZ_det3x3(isodata.vec[obj][0], isodata.vec[obj][1], isodata.vec[obj][2]);
      if (det > 0) {
	fprintf(stderr,"isosurface spanning vectors orientation is right\n");
	isodata.cell_orientation = XC_RIGHT;
      } else {
	fprintf(stderr,"isosurface spanning vectors orientation is left\n");
	isodata.cell_orientation = XC_LEFT;
      }	
    }
  } else {
    isodata.vec[obj][2].x = 0.0;
    isodata.vec[obj][2].y = 0.0;
    isodata.vec[obj][2].z = 0.0;
  }

  nobj = obj;
  if ( isodata.dim[obj] == 3 && obj == ISOOBJ_BASE ) {
    /* get also the vectors for PLANE1|2|3

       PLANE1 == vec0 x vec1; -+
       PLANE2 == vec2 x vec0;  |-> origin je point0
       PLANE3 == vec1 x vec2; -+
    */    
    isodata.vec[ISOOBJ_PLANE1][0] = isodata.vec[obj][0];
    isodata.vec[ISOOBJ_PLANE1][1] = isodata.vec[obj][1];
    isodata.vec[ISOOBJ_PLANE1][2] = isodata.vec[obj][2];

    isodata.vec[ISOOBJ_PLANE2][0] = isodata.vec[obj][2];
    isodata.vec[ISOOBJ_PLANE2][1] = isodata.vec[obj][0];
    isodata.vec[ISOOBJ_PLANE2][2] = isodata.vec[obj][1];

    isodata.vec[ISOOBJ_PLANE3][0] = isodata.vec[obj][1];
    isodata.vec[ISOOBJ_PLANE3][1] = isodata.vec[obj][2];
    isodata.vec[ISOOBJ_PLANE3][2] = isodata.vec[obj][0];

    nobj = (ISOOBJ_PLANE3 > obj) ? ISOOBJ_PLANE3 : obj;
  }

  /* get the normals for isoplanes */
  for (i = obj; i<= nobj; i++) {
    isodata.colnml[i] = 
      VertexNormal(isodata.vec[i][0], isodata.vec[i][1]);
    /* normal must face toward positive Z; that means positive volume with
     * respect to (0,0,1), but take instead rather vector (1.0, 1.0, 10.0), to
     * be secure if z-comp of normal is zero 
     */
    /* scalar = Zvec[0] * isodata.colnml[i].x + \
       Zvec[1] * isodata.colnml[i].y + Zvec[2] * isodata.colnml[i].z;*/
    /* 
       if ( !IsIsoPlane123( i ) ) {
       if ( (scalar < 0.0 && fabsf(scalar) > 1.e-6) || 
       (isodata.colnml[i].z < 0.0 && fabsf(scalar) < 1.e-6) ) {
       isodata.colnml[i].x *= -1.0;
       isodata.colnml[i].y *= -1.0;
       isodata.colnml[i].z *= -1.0;
       }
    */

    /* 
       for PLANE123, it seems that i mixed up the vector products 
       (instead on counterclockwise I made clockwise, so correct that 
    */
    /*if ( (scalar > 0.0 && fabsf(scalar) > 1.e-6) || 
      (isodata.colnml[i].z > 0.0 && fabsf(scalar) < 1.e-6) ) {
      isodata.colnml[i].x *= -1.0;
      isodata.colnml[i].y *= -1.0;
      isodata.colnml[i].z *= -1.0;
      }
    */
    normalizepvf( &(isodata.colnml[i].x), &(isodata.colnml[i].y), 
		  &(isodata.colnml[i].z) );
  }

  /* everything is OK */
  if ( !(isostate.stateflag & ISO_POINTS) ) isostate.stateflag |= ISO_POINTS;  

  return TCL_OK;
}



/******************************************************************************
 * XC_IsodataCmd --> implementation of 'xc_isodata' custom Tcl command
 * -----------------------
 * Usage: xc_isodata (<frame[i]_range>,i=ISODATA_MAXSTACK,0,-1)
 *                                     ^^^^^^^^^^^^^^^^^^^^^^^
 *        or
 *
 *        xc_isodata make
 *        ------------------------------------------------------------
 *
 *        syntax of frame[i]_range -->   x     .... x-th frame
 *                                       x-y   .... from x-th to y-th frames 
 *
 * --------------------------------------------------------------------------
 * xc_isofiles require:
 *                      pre-calling of: xc_isostack
 *                                      xc_isosign
 *                                      xc_isofiles
 *
 * To specify all iso-data, somethimes it's required to call xc_isodata 
 * several times, thatwhy it is needed to know when we are done by 
 * xc_isodata. "xc_isodata make" command is implemented for specifying 
 * the end of xc_isodata.
 * When we call xc_isodata the first time initializations is performed before 
 * before anything happen. By specifying "xc_isodata make" we also set 
 * initialization flag to "uninitialized", so next call of xc_isodata will
 * perform initialization at a first place.
 */
int 
XC_IsodataCmd(ClientData clientData, Tcl_Interp *interp,
	      int argc, const char *argv[])
{
  int i, j, k, first;
  static int which_time = 0;
  static int framenum[2 * ISODATA_MAXFRAME][ISODATA_MAXSTACK][2];
  char fint[3], sint[3];

  /* initialization for framenum[][]; just when we enter the first time */
  if (which_time == 0) {
    if ( !(isostate.stateflag & 
	   (ISO_STACK | ISO_SIGN | ISO_FILES| ISO_POINTS)) ) {
      Tcl_SetResult(interp, "\"xc_isostack\", \"xc_isosign\", \"xc_isofiles\" and \"xc_isopoints\" should be called before \"xc_isodata\" command", TCL_STATIC);
      xcIsoError();
      return TCL_ERROR;
    }
    /*
     * malloc gridvertex or plvertex[ISOOBJ_BASE] !!!
     */
    if ( isodata.dim[ISOOBJ_BASE] == 3 && isostate.gridvertex_malloc == 0 ) {
      gridvertex = xcMallocGRIDVERTEX( grd.nx, grd.ny, grd.nz );
      isostate.gridvertex_malloc = 1;
    } 
    else if ( isodata.dim[ISOOBJ_BASE] == 2 && 
	      isostate.plvertex_malloc[ISOOBJ_BASE] == 0 ) {
      plvertex[ISOOBJ_BASE] = xcMallocPLANEVERTEX( grd.nx, grd.ny );
      isostate.plvertex_malloc[ISOOBJ_BASE] = 1;
      for (i=0; i<ISOLINE_MAXLEVEL; i++) {
	isostate.max_n_isoline2D[ISOOBJ_BASE][i] = grd.nx;
	isoline2D[ISOOBJ_BASE].segment[i] = 
	  xcMallocLINE( isostate.max_n_isoline2D[ISOOBJ_BASE][i] ); 
      }
      isostate.isoline2D_malloc[ISOOBJ_BASE] = 1;
    }

    for (k=0; k < ISODATA_MAXFILES; k++)
      for (i=0; i < ISODATA_MAXSTACK; i++) 
	for (j=0; j < 2; j++)
	  framenum[k][i][j] = -1;
    /* initialize ReadBlock0 function */
    ReadBlock0(RB0_INIT, 0, 0, 0, 0);
  } 

  /***************************************************************************
   * if isodata.dim[ISOOBJ_BASE] == 3 then 
   *                            malloc space for "vertex" & "triangles"        
   *
if ( isodata.dim[ISOOBJ_BASE] == 3 ) {
isostate.max_n_triangl = isostate.max_n_vertex = 
MAX_ISOSURFLEVELS * 2 * newgrd.nx * newgrd.ny * newgrd.nz;
if ( !isostate.vertex_malloc ) {
vertex_address = 
malloc( (size_t) sizeof(VERTEX) * isostate.max_n_vertex );
isostate.vertex_malloc = 1;
if (!vertex_address) xcError("allocation error for vertex_adress");
}
if ( !isostate.triangl_malloc ) {
triangl_address = 
malloc( (size_t) sizeof(TRIANGLE) * isostate.max_n_triangl );
isostate.triangl_malloc = 1;
if (!vertex_address) xcError("allocation error for triangl_adress");
}
}
  */

  /***************************************************************************/

  /* is usage of xc_isodata command correct !?? */
  if ( argc == 2 && strcmp(argv[1],"make") == 0 ) {   
    /* framenum[][0][0] must be entered in sorted fashion, check that */
    for (i=0; i<which_time-1; i++)
      if ( framenum[i+1][0][0] < framenum[i][0][0] ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"sequence of \"xc_isodata (<frame[i]_range>,i=%d,0,-1)\" was not used properly, because stack-0 frames must be entered in sorted fashion",ISODATA_MAXSTACK);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	xcIsoError();
	return TCL_ERROR;
      }

    /* fwrite origin & vectors to "bin_vertex_fp" */
    /*
      fwrite(isodata.points[ISOOBJ_BASE][0], 
      sizeof(double), 3, isodata.bin_vertex_fp);
      fwrite(isodata.vec[ISOOBJ_BASE][0], sizeof(XYZ), 3, 
      isodata.bin_vertex_fp);
    */

    if ( !ReadIsoData(framenum, which_time) ) {
      char rss[1024];
      snprintf(rss, sizeof(rss),"%s", isodata.error);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      xcIsoError();
      return TCL_ERROR;
    }
    /* now assign the isodata.min */
    isodata.mintol = fabs(isodata.max * 1.0e-6);
    if ( fabs(isodata.min) > fabs(isodata.max) )
      isodata.mintol = fabs(isodata.min * 1.0e-6);

    for (k=0; k < ISODATA_MAXFILES; k++)
      for (i=0; i < ISODATA_MAXSTACK; i++) 
	for (j=0; j < 2; j++)
	  framenum[k][i][j] = -1;
    which_time = 0; /* reset the which_time, because, next time we enter 
		     * in this routine, we will start from begining */
    return TCL_OK;
  }
  else if ( argc == 2 && strcmp(argv[1],"make") != 0 ) {
    /* wrong usage of xc_isodata command */
    char rss[1024];
    snprintf(rss, sizeof(rss),"wrong usage of \"xc_isodata\" command, should be: \nxc_isodata (<frame[i]_range>,i=%d,0,-1)\n  or\nxc_isodata make",ISODATA_MAXSTACK);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    xcIsoError();
    return TCL_ERROR;
  }
  else if ( argc != ISODATA_MAXSTACK + 1 ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wrong # of arguments, while executing \"xc_isodata\" command, should be \"xc_isodata (<frame[i]_range>,i=%d,0,-1)\"",ISODATA_MAXSTACK);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    xcIsoError();
    return TCL_ERROR;
  } 

  /* we come here only if NOT "xc_isodata make" */
  for (i=1; i < argc; i++) {
    int f3, f31;
    if ( strstr(argv[i],"-") == NULL ) {
      /* number, not range, was specified */
      if ( Tcl_GetInt(interp, argv[i], &(framenum[which_time][i-1][0])) == 
	   TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing \"xc_isodata (<frame[i]_range>,i=%d,0,-1)\" command", argv[i], ISODATA_MAXSTACK);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      /* maybe number is not in range, check that */
      if ( i == 1 ) {
	f3 = framenum[which_time][0][0];
	f31 = f3;
      }
      for (j=f3; j<=f31; j++)
	if ( framenum[which_time][i-1][0] > isodata.nframe[j][ISODATA_MAXSTACK - i] ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"specified frame - %d of stack %d is out of range, frame should be lower than %d", 
		   framenum[which_time][i-1][0], i-1, 
		   isodata.nframe[j][ISODATA_MAXSTACK - i] + 1);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      framenum[which_time][i-1][1] = -1;
    } else {
      /* range was specified */
      first = strcspn(argv[i],"-");
      strncpy(fint,argv[i],first);
      /*fint[first] = (char) NULL;*/
      fint[first] = ' ';
      strcpy(sint,&argv[i][first+1]);
      if ( Tcl_GetInt(interp, fint, &(framenum[which_time][i-1][0])) == 
	   TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing \"xc_isodata (<frame[i]_range>,i=%d,0,-1)\" command", argv[i], ISODATA_MAXSTACK);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      if ( Tcl_GetInt(interp, sint, &(framenum[which_time][i-1][1])) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing \"xc_isodata (<frame[i]_range>,i=%d,0,-1)\" command", argv[i], ISODATA_MAXSTACK);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      /* framenum[][][1] must be > framenum[][][0]; check that */
      if ( framenum[which_time][i-1][0] >= framenum[which_time][i-1][1] ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"range in \"xc_isodata\" command not specified properly, should be a-b, but got %s",argv[i]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      /* maybe frames specified are out of range, check that */
      if ( i == 1 ) {
	f3  = framenum[which_time][0][0];
	f31 = framenum[which_time][0][1];
      }
      for (j=f3; j<=f31; j++)
	if ( framenum[which_time][i-1][1] > isodata.nframe[j][ISODATA_MAXSTACK - i] ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"specified range of frames - [%d-%d] of stack %d is out of range, frame should be lower than %d", 
		   framenum[which_time][i-1][0], framenum[which_time][i-1][1], 
		   i-1, isodata.nframe[j][ISODATA_MAXSTACK - i]+1);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}      
    }
  }

  which_time++; /* which time trace the number of times we get in this func */
  if ( !(isostate.stateflag & ISO_DATA) ) isostate.stateflag |= ISO_DATA;    
  return TCL_OK;
}



/******************************************************************************
 * XC_IsosurfCmd --> implementation of 'xc_isosurf' custom Tcl command
 * -----------------------
 * Usage: xc_isosurf <toglName> \
 *                   -drawstyle  wire|solid  \
 *                   -shademodel smooth|flat \
 *                   -transparency on|off \
 *                   -render now|after \
 *                   -isosurf none| xxx
 *                   
 *
 * Note: default flag for options: 
 *                                -drawstyle    = "wire"
 *                                -shademodel   = "smooth"
 *                                -trensparency = "off"
 *                                -render       = after
 *                                -isosurf      = xxx
 *
 * if -isosurf == none; then VPf.isosurf is set to GL_FALSE && the isosurface
 *                      is not rendered
 */                     
int 
XC_IsosurfCmd(ClientData clientData, Tcl_Interp *interp,
	      int argc, const char *argv[])
{
  struct Togl *togl;
  int i;
  int render_isosurf = 1;
  logical render = 0;

  if ( (argc % 2) != 0 ) {
    Tcl_SetResult(interp, "Usage: xc_isosurf <toglName> -drawstyle wire|solid -shademodel smooth|flat -transparency on|off -render now|after -isosurf none|xxx", TCL_STATIC);
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

  /* it is not allowed to be in XC_2D mode when rendering isosurface;
     if in XC_2D mode return silently */
  if (dimType == XC_2D) return TCL_OK;

  /* isodata.dim must be 3; if not -> ERROR */
  if ( isodata.dim[ISOOBJ_BASE] != 3 ) {
    Tcl_SetResult(interp, "\"xc_isopoints 3D ...\" command must be used before evoking \"xc_isoplane\" command", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( !(isostate.stateflag & 
	 (ISO_STACK | ISO_SIGN | ISO_FILES | ISO_POINTS | 
	  ISO_DATA | ISO_POLYGONISE)) ) {
    Tcl_SetResult(interp, "\"xc_isostack\", \"xc_isosign\", \"xc_isofiles\", \"xc_isopoints\", \"xc_isodata\" and \"xc_iso poligonise\" should be called before \"xc_isosurf\" command", TCL_STATIC);
    return TCL_ERROR;
  }

  /****************/
  /* read options */
  for (i=2; i<argc; i+=2) {
    if ( strcmp(argv[i],"-drawstyle") == 0 ) {
      if ( strcmp(argv[i+1],"wire") == 0 ) isoDisp.drawstyle = ISOSURF_WIRE;
      else if ( strcmp(argv[i+1],"solid") == 0 ) isoDisp.drawstyle = ISOSURF_SOLID;
      else if ( strncmp(argv[i+1],"dot",3) == 0 ) isoDisp.drawstyle = ISOSURF_DOT;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown value \"%s\" for -drawstyle option, must be either \"wire\" or \"solid\" or \"dot\", while executing %s %s",
		 argv[i+1], argv[0], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else if ( strcmp(argv[i],"-shademodel") == 0 ) {
      if ( strcmp(argv[i+1],"smooth") == 0) isoDisp.shademodel = GL_SMOOTH;
      else if ( strcmp(argv[i+1],"flat") == 0 ) isoDisp.shademodel = GL_FLAT;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown value \"%s\" for -shademodel option, must be either \"smooth\" or \"flat\", while executing %s %s",
		 argv[i+1], argv[0], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else if ( strcmp(argv[i],"-transparency") == 0 ) {
      if ( strcmp(argv[i+1],"off") == 0 ) 
	isoDisp.transparent = ISOSURF_TRANSP_OFF;
      else if ( strcmp(argv[i+1],"on") == 0 ) 
	isoDisp.transparent = ISOSURF_TRANSP_ON;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown value \"%s\" for -transparency option, must be either \"off\" or \"on\", while executing %s %s",
		 argv[i+1], argv[0], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else if ( strcmp(argv[i],"-render") == 0 ) {
      if ( strcmp(argv[i+1],"now") == 0 ) 
	render = 1;
      else if ( strcmp(argv[i+1],"after") == 0 ) 
	render = 0;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown value \"%s\" for -render option, must be either \"now\" or \"after\", while executing %s %s",
		 argv[i+1], argv[0], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else if ( strcmp(argv[i],"-isosurf") == 0 ) {
      if ( strcmp(argv[i+1],"none") == 0 ) {
	render_isosurf = 0;
      }
    } else {
      char rss[1024];
      snprintf(rss, sizeof(rss),"unknown option \"%s\", must one of -drawstyle, -shademodel, -transparency, -render and -isosurf, while executing %s %s",
	       argv[i+1], argv[0], argv[1]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  VPf.isosurf[ISOOBJ_BASE]    = GL_TRUE;
  if (!render_isosurf) {
    VPf.isosurf[ISOOBJ_BASE]  = GL_FALSE;
  }
  /* turn off all VPf factors, that contadict with VPf.isosurf */
  VPf.colorplane[ISOOBJ_BASE] = GL_FALSE;
  VPf.isoline[ISOOBJ_BASE]    = GL_FALSE;

  /* because isosurface may be greater than structure, 
     take care of projection */
  if (is.stickmode && !is.ballmode) xcMakeProjection3D("sticks");
  if (is.ballmode) xcMakeProjection3D("balls");
  if (is.spacefillmode) xcMakeProjection3D("space");

  /* call (*xcDisplay)() if -render option was "now" */
  MVf.isosize = SetIsosurf_VPf();
  if ( render ) {
    Togl_PostRedisplay(togl);
  }

  return TCL_OK;
}


/******************************************************************************
 * XC_IsoplaneCmd --> implementation of 'xc_isoplane' custom Tcl command
 * -----------------------
 * Usage: xc_isoplane <toglName> <object> \
 *                   -planetype  colorplane|isoline|both \
 *                   -transparency on|off
 *                   -render now|after
 *
 * Note: default flag for options: 
 *                                -planetype    = "colorplane"
 *                                -transparency = "off"
 *                                -render       = "after"
 * Note #2: isoline is not implemented yet
 */
int 
XC_IsoplaneCmd(ClientData clientData, Tcl_Interp *interp,
	       int argc, const char *argv[])
{
  struct Togl *togl;
  int i, obj, planeflag, lineflag;
  logical render = 0;

  if ( (argc % 2) != 1 ) {
    Tcl_SetResult(interp, "Usage: xc_isoplane <toglName> <object> -planetype colorplane|isoline -transparency on|off -render now|after", TCL_STATIC);
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

  /* it is not allowed to be in XC_2D mode when rendering isoplane;
     if in XC_2D mode return silently */
  if (dimType == XC_2D) return TCL_OK;

  /*
   * get object type
   */
  if ( Tcl_GetInt(interp, argv[2], &obj ) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" while executing %s %s %s ...", argv[2], argv[0], argv[1], argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( obj < 0 || obj >= MAX_ISOOBJECTS ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"object type \"%d\" is out of range", obj);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* isodata.dim[ISOOBJ_BASE] must be 2; if not -> ERROR */
  if ( obj == ISOOBJ_BASE && isodata.dim[ISOOBJ_BASE] != 2 ) {
    Tcl_SetResult(interp, "\"xc_isopoints 0 2D ...\" command must be used before evoking \"xc_isoplane\" command", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( !(isostate.stateflag & 
	 (ISO_STACK | ISO_SIGN | ISO_FILES | ISO_POINTS | ISO_DATA)) ) {
    Tcl_SetResult(interp, "\"xc_isostack\", \"xc_isosign\", \"xc_isofiles\", \"xc_isopoints\", and \"xc_isodata\" should be called before \"xc_isosurf\" command", TCL_STATIC);
    return TCL_ERROR;
  }

  /*
   * set default options
   */
  isoplaneDisp[obj].transparent = ISOSURF_TRANSP_OFF;
  VPf.colorplane[obj] = GL_TRUE;
  VPf.isoline[obj]    = GL_FALSE;

  planeflag = ISO_COLORPLANE;
  lineflag  = ISO_ISOLINE;
  if ( obj == ISOOBJ_PLANE1 ) {
    planeflag = ISO_COLORPLANE1;
    lineflag  = ISO_ISOLINE1;
  } else if ( obj == ISOOBJ_PLANE2 ) {
    planeflag = ISO_COLORPLANE2;
    lineflag  = ISO_ISOLINE2;
  } else if ( obj == ISOOBJ_PLANE3 ) {
    planeflag = ISO_COLORPLANE3;
    lineflag  = ISO_ISOLINE3;
  } 

  /****************/
  /* read options */
  for (i=3; i<argc; i+=2) {
    if ( strcmp(argv[i],"-planetype") == 0 ) {
      if ( strcmp(argv[i+1],"colorplane") == 0 ) {
	/* check for ISO_COLORPLANE flag */
	if ( !(isostate.stateflag & planeflag) ) {
	  Tcl_SetResult(interp, "\"xc_iso isoplane ... colorplane\" should be called before \"xc_isoplane -planetype colorplane\" command", TCL_STATIC);
	  return TCL_ERROR;
	}	
	VPf.colorplane[obj] = GL_TRUE;
	VPf.isoline[obj]    = GL_FALSE;
      }
      else if ( strcmp(argv[i+1],"isoline") == 0 ) {
	if ( !(isostate.stateflag & lineflag) ) {
	  Tcl_SetResult(interp, "\"xc_iso isoplane ... isoline\" should be called before \"xc_isoplane -planetype isoline\" command", TCL_STATIC);
	  return TCL_ERROR;
	}
	VPf.isoline[obj]    = GL_TRUE;
	VPf.colorplane[obj] = GL_FALSE;
      } else if ( strcmp(argv[i+1], "both") == 0 ) {
	if ( !(isostate.stateflag & (planeflag | lineflag)) ) {
	  Tcl_SetResult(interp, "\"xc_iso isoplane ... isoline\" and \"xc_iso isoplane ... colorplane\"should be called before \"xc_isoplane -planetype isoline\" command", TCL_STATIC);
	  return TCL_ERROR;
	}
	VPf.isoline[obj]    = GL_TRUE;
	VPf.colorplane[obj] = GL_TRUE;
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"invalid argument %s for -planetype option, must be one of \"colorplane\", \"isoline\" or \"both\", while executing xc_isoplane %s ...",
		 argv[i+1],argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else if ( strcmp(argv[i],"-transparency") == 0 ) {
      if ( strcmp(argv[i+1],"off") == 0 ||
	   strcmp(argv[i+1],"0") == 0) {
	isoplaneDisp[obj].transparent = ISOSURF_TRANSP_OFF;
      }
      else if ( strcmp(argv[i+1],"on") == 0 ||
		strcmp(argv[i+1],"1")  == 0 ) {
	isoplaneDisp[obj].transparent = ISOSURF_TRANSP_ON;
      }
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"invalid argument %s for -transparency option, must be one of \"on\" or \"off\", while executing xc_isoplane %s ...",
		 argv[i+1],argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else if ( strcmp(argv[i],"-render") == 0 ) {
      if ( strcmp(argv[i+1],"now") == 0 ) 
	render = 1;
      else if ( strcmp(argv[i+1],"after") == 0 ) 
	render = 0;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown value \"%s\" for -render option, must be either \"now\" or \"after\", while executing %s %s",
		 argv[i+1], argv[0], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss),"unknown option %s, must be one of \"-planetype\", \"-transparency\" or \"-render\", while executing xc_isoplane %s ...",argv[i],argv[1]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }    
  }

  if ( obj == ISOOBJ_BASE ) /* colorplane is BASE object, so disable isosurf */
    VPf.isosurf[ISOOBJ_BASE] = GL_FALSE;

  /* call (*xcDisplay)() if -render option was "now" */
  MVf.isosize = SetIsosurf_VPf(); 
  if ( render ) {
    Togl_PostRedisplay(togl);
  }

  return TCL_OK;
}


/******************************************************************************
 * XC_IsoexpandCmd --> implementation of 'xc_isoexpand' custom Tcl command
 * -----------------------
 * Usage:
 * xc_isoexpand <toglName> <object> \
 *                         -repeattype default|convcell|primcell \
 *                         -shape      default|parapipedal|hexagonal \
 *                         -expand     whole|{nx ny nz}|... 
 *                         -render           now|after
 * default option values:
 *                        -repeattype       default
 *                        -shape            default
 *                        -expand           whole
 *                        -render           after
 * XC_IsoexpandCmd first check if ISO_POLYGONISE or ISO_COLORPLANE or
 * ISO_ISOLINE is turned on. So there is no need in specifing this option 
 * separatly. If someday it will be possible to render simultaneously then
 * this option will be needed.
 */
int 
XC_IsoexpandCmd(ClientData clientData, Tcl_Interp *interp,
		int argc, const char *argv[])
{
  struct Togl *togl;
  int obj, colorflag, lineflag;
  int i, j;
  int render = 0;
  unsigned int flag = ISO_NULL, flag1 = ISO_NULL;
  int repeat = COM_ISOEXPAND_REPEAT_DEFAULT;
  GetComOption expand = { COM_ISOEXPAND_EXPAND_WHOLE };

  if ( (argc % 2) != 1 ) {
    Tcl_SetResult(interp, "Usage: xc_isoexpand <toglName> <object> -repeattype default|convcell|primcell -shape default|parapipedal| hexagonal -expand whole|{nx ny nz}", TCL_STATIC);
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

  if ( Tcl_GetInt(interp, argv[2], &obj) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\" while executing %s %s %s ...", argv[2], argv[0], argv[1], argv[2]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( obj < 0 || obj >= MAX_ISOOBJECTS ) {
    char rss[1024];
    snprintf(rss, sizeof(rss),"object type \"%d\" is out of range", obj);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }


  colorflag = ISO_COLORPLANE;
  lineflag  = ISO_ISOLINE; 
  switch( obj ) 
    {
    case ISOOBJ_PLANE1:
      colorflag = ISO_COLORPLANE1;
      lineflag  = ISO_ISOLINE1; 
      break;
    case ISOOBJ_PLANE2:
      colorflag = ISO_COLORPLANE2;
      lineflag  = ISO_ISOLINE2; 
      break;
    case ISOOBJ_PLANE3:
      colorflag = ISO_COLORPLANE3;
      lineflag  = ISO_ISOLINE3; 
      break;
    }

  /* maybe xc_isoexpand was called two early */
  if ( isodata.dim[obj] == 2 ) {
    flag  = ISO_STACK | ISO_SIGN | ISO_FILES | ISO_POINTS | 
      ISO_DATA | colorflag;
    flag1 = ISO_STACK | ISO_SIGN | ISO_FILES | ISO_POINTS | 
      ISO_DATA | lineflag;
  }
  if ( isodata.dim[obj] == 3 )
    flag = ISO_STACK | ISO_SIGN | ISO_FILES | ISO_POINTS | 
      ISO_DATA | ISO_POLYGONISE;

  if ( !(isostate.stateflag & flag) && !(isostate.stateflag & flag1) ) {
    Tcl_SetResult(interp, "\"xc_isostack\", \"xc_isosign\", \"xc_isofiles\", \"xc_isopoints\", \"xc_isodata\" and \"xc_iso colorplane\" or \"xc_iso polygonise\" should be called before \"xc_isoexpand\" command", TCL_STATIC);
    return TCL_ERROR;
  }  

  /* if we are in XC_2D mode return silently */
  if ( dimType == XC_2D ) return TCL_OK;

  /* load default value for shape */
  isoexpand[obj].shape = COM_ISOEXPAND_SHAPE_DEFAULT;
  for (j=0; j<3; j++)
    isoexpand[obj].irepvec[j] = xcr.nunit[j];

  /****************/
  /* read options */
  for (i=3; i<argc; i+=2) {
    if ( strncmp(argv[i],"-repeat",7) == 0 ) {
      if ( strcmp(argv[i+1],"default") == 0 ) 
	repeat = COM_ISOEXPAND_REPEAT_DEFAULT;
      else if ( strcmp(argv[i+1],"convcell") == 0 ) 
	repeat = COM_ISOEXPAND_REPEAT_CONVCELL;
      else if ( strcmp(argv[i+1],"primcell" ) == 0 ) 
	repeat = COM_ISOEXPAND_REPEAT_PRIMCELL;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown value \"%s\" for -repeattype option, must be one of \"default\", \"convcell\", \"primcell\", while executing %s %s ...", argv[i+1], argv[0], argv[1]); 
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else if ( strcmp(argv[i],"-shape") == 0 ) {
      if ( strcmp(argv[i+1],"default") == 0 ) 
	isoexpand[obj].shape = COM_ISOEXPAND_SHAPE_DEFAULT;
      else if ( strncmp(argv[i+1],"para",4) == 0 ) 
	isoexpand[obj].shape = COM_ISOEXPAND_SHAPE_PARAPIPEDAL;      
      else if ( strncmp(argv[i+1],"hexa",4) == 0 )
	isoexpand[obj].shape = COM_ISOEXPAND_SHAPE_HEXAGONAL;      
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown value \"%s\" for -shape option, must be one of \"default\", \"parapipedal\", \"hexagonal\", while executing %s %s ...", argv[i+1], argv[0], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else if ( strcmp(argv[i],"-expand") == 0 ) {
      if ( strcmp(argv[i+1],"whole") == 0 ) {
	isoexpand[obj].expand = COM_ISOEXPAND_EXPAND_WHOLE;
	for (j=0; j<3; j++)
	  isoexpand[obj].irepvec[j] = xcr.nunit[j];
      }
      else {
	int type = 3;
	/* take care of list:
	 * if isodata.dim[obj] == 2 --> {nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, na, nb, nc} 
	 * if isodata.dim[obj] == 3 --> {nx, ny, nz}
	 */
	if ( isodata.dim[obj] == 2 ) type = XC_GET_INT12;
	if ( isodata.dim[obj] == 3 ) type = XC_GET_INT3;
	if ( !xcSplitList( type, interp, argv + i + 1, &expand) )
	  return TCL_ERROR;			
	if ( isodata.dim[obj] == 2 ) {
	  /* usage:
	     -expand {<vec1> <vec2> <vec3> n1 n2 n3}
	     <vecx> == (nx, ny, nz)
	  */
	  for (j=0; j<3; j++) {	    
	    isoexpand[obj].transl[0][j] = iroundf(expand.vec[j]);
	    isoexpand[obj].transl[1][j] = iroundf(expand.vec[j+3]);
	    isoexpand[obj].transl[2][j] = iroundf(expand.vec[j+6]);
	  }
	  isoexpand[obj].irepvec[0] = iroundf(expand.vec[9]);
	  isoexpand[obj].irepvec[1] = iroundf(expand.vec[10]);
	  isoexpand[obj].irepvec[2] = iroundf(expand.vec[11]);
	}
	else if ( isodata.dim[obj] == 3 ) {
	  for (j=0; j<3; j++) 
	    isoexpand[obj].irepvec[j] = iroundf(expand.vec[j]);
	}
	isoexpand[obj].expand = COM_ISOEXPAND_EXPAND_LIST;
      }
    }
    else if ( strcmp(argv[i], "-render") == 0 ) {
      if ( strcmp(argv[i+1], "now") == 0 ) {
	render = 1;
      } else if ( strcmp(argv[i+1], "after") == 0 ) {
	render = 0;
      } else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown -render value \"%s\", should be now or after", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    else {
      char rss[1024];
      snprintf(rss, sizeof(rss),"unknown option \"%s\", must be one of \"-repeat\",\"-shape\", \"-expand\", \"-render\", while executing %s %s ...", 
	       argv[i], argv[0], argv[1]); 
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  if ( (repeat == COM_ISOEXPAND_REPEAT_DEFAULT && 
	xcr.celltype == XCR_PRIMCELL) ||
       repeat == COM_ISOEXPAND_REPEAT_PRIMCELL ) {
    if ( isodata.dim[obj] == 3 ) 
      xcMat44Copyd( isoexpand[obj].rep_vec, vec.prim, 3, 3 );
    else if ( isodata.dim[obj] == 2 ) {
      if ( isoexpand[obj].expand == COM_ISOEXPAND_EXPAND_LIST ) 
	for (i=0; i<3; i++) {
	  isoexpand[obj].rep_vec[0][i] = 
	    vec.prim[0][i] * (double) isoexpand[obj].transl[0][0] +
	    vec.prim[1][i] * (double) isoexpand[obj].transl[0][1] +
	    vec.prim[2][i] * (double) isoexpand[obj].transl[0][2];
	  isoexpand[obj].rep_vec[1][i] =
	    vec.prim[0][i] * (double) isoexpand[obj].transl[1][0] +
	    vec.prim[1][i] * (double) isoexpand[obj].transl[1][1] +
	    vec.prim[2][i] * (double) isoexpand[obj].transl[1][2];
	  isoexpand[obj].rep_vec[2][i] = 0.0;
	}
      else {
	for (i=0; i<3; i++)
	  for (j=0; j<3; j++)
	    isoexpand[obj].rep_vec[i][j] = vec.prim[i][j];
      }
    }
  }
  if ( (repeat == COM_ISOEXPAND_REPEAT_DEFAULT && 
	xcr.celltype == XCR_CONVCELL ) ||
       repeat == COM_ISOEXPAND_REPEAT_CONVCELL ) {
    if ( isodata.dim[obj] == 3 )
      xcMat44Copyd( isoexpand[obj].rep_vec, vec.conv, 3, 3 ); 
    else if ( isodata.dim[obj] == 2 ) {
      if ( isoexpand[obj].expand == COM_ISOEXPAND_EXPAND_LIST )
	for(i=0; i<3; i++) {
	  isoexpand[obj].rep_vec[0][i] = 
	    vec.conv[0][i] * (double) isoexpand[obj].transl[0][0] +
	    vec.conv[1][i] * (double) isoexpand[obj].transl[0][1] +
	    vec.conv[2][i] * (double) isoexpand[obj].transl[0][2];
	  isoexpand[obj].rep_vec[1][i] =
	    vec.conv[0][i] * (double) isoexpand[obj].transl[1][0] +
	    vec.conv[1][i] * (double) isoexpand[obj].transl[1][1] +
	    vec.conv[2][i] * (double) isoexpand[obj].transl[1][2];
	  isoexpand[obj].rep_vec[2][i] = 
	    vec.conv[0][i] * (double) isoexpand[obj].transl[2][0] +
	    vec.conv[1][i] * (double) isoexpand[obj].transl[2][1] +
	    vec.conv[2][i] * (double) isoexpand[obj].transl[2][2];
	}
      else {
	for (i=0; i<3; i++)
	  for (j=0; j<3; j++)
	    isoexpand[obj].rep_vec[i][j] = vec.conv[i][j];
      }
    }
  }

  /* now according to GENGEOM for slabs, polymer & molecule I get
   * (nx,ny,0), (nx,0,0) & (0,0,0) respectively, but the lowest
   * possible vector here is (1,1,1) -> correct that
   */
  for (j=0; j<3; j++)
    if ( isoexpand[obj].irepvec[j] == 0 ) isoexpand[obj].irepvec[j] = 1;

  /* isosurface's space can be greater that structure's space, so we must 
   * update MVf.structsize to allow all isosurface to be rendered
   */
  MVf.isosize = SetIsosurf_VPf();
  if (is.stickmode && !is.ballmode) xcMakeProjection3D("sticks");
  if (is.ballmode) xcMakeProjection3D("balls");
  if (is.spacefillmode) xcMakeProjection3D("space");

  if ( !(isostate.stateflag & ISO_EXPAND) ) isostate.stateflag |= ISO_EXPAND;  

  if (render) Togl_PostRedisplay(togl);

  return TCL_OK;
}


static float
SetIsosurf_VPf(void)
{
  register int i;
  float isosize, max_isosize = 0.0;
  float  orig[3], vec0[3], vec1[3], vec2[3];

  for (i=ISOOBJ_BASE; i<MAX_ISOOBJECTS; i++) {
    orig[0] = -isodata.points[i][0][0] + mx;
    orig[1] = -isodata.points[i][0][1] + my;
    orig[2] = -isodata.points[i][0][2] + mz;

    vec0[0]=(float)(isoexpand[i].irepvec[0]-1) * isoexpand[i].rep_vec[0][0] +
      isodata.vec[i][0].x;
    vec0[1]=(float)(isoexpand[i].irepvec[0]-1) * isoexpand[i].rep_vec[0][1] +
      isodata.vec[i][0].y;
    vec0[2]=(float)(isoexpand[i].irepvec[0]-1) * isoexpand[i].rep_vec[0][2] +
      isodata.vec[i][0].z;

    vec1[0]=(float)(isoexpand[i].irepvec[1]-1) * isoexpand[i].rep_vec[1][0] +
      isodata.vec[i][1].x;
    vec1[1]=(float)(isoexpand[i].irepvec[1]-1) * isoexpand[i].rep_vec[1][1] +
      isodata.vec[i][1].y;
    vec1[2]=(float)(isoexpand[i].irepvec[1]-1) * isoexpand[i].rep_vec[1][2] +
      isodata.vec[i][1].z;

    vec2[0]=(float)(isoexpand[i].irepvec[2]-1) * isoexpand[i].rep_vec[2][0] +
      isodata.vec[i][2].x;
    vec2[1]=(float)(isoexpand[i].irepvec[2]-1) * isoexpand[i].rep_vec[2][1] +
      isodata.vec[i][2].y;
    vec2[2]=(float)(isoexpand[i].irepvec[2]-1) * isoexpand[i].rep_vec[2][2] +
      isodata.vec[i][2].z;

    isosize = DetermineParapipedSize( vec0, vec1, vec2, orig );
    if ( isosize > max_isosize ) max_isosize = isosize;
  }

  return max_isosize;
}
