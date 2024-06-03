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
* Source: $XCRYSDEN_TOPDIR/C/crySurfArgs.c
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


int
VerifyReadSurfOpt(const char *option, size_t n, 
		  struct READSURF_OPT *opts, char *argv[])
{
  char *opt;
  int i;

  i=0;
  opt = opts[i].option;
  while( opt ) {
    if ( strncmp( option, opt, n ) == 0 ) {
      if ( opts[i].flag == XC_FORBIDDEN ) {
	fprintf(stderr, "option %s not allowed !!!\nError while executing %s %s", option, argv[0], argv[1]);
	return XC_ERROR;
      } else {
	return XC_OK;
      }
    }
    opt = opts[++i].option;
  }

  fprintf(stderr, "option %s not known !!!\nError while executing %s %s", 
	  option, argv[0], argv[1]);
  return XC_ERROR;
}


int
cryReadSurfARGS(ClientData clientData,Tcl_Interp *interp,
		int argc, char *argv[], NEW_WIN_CONTEXT *wc, MOL_SURF *mols,
		struct READSURF_OPT *opt) 
{
  register int i, j;

  for (i=4; i<argc; i+=2) {
    /* -TYPE */
    if ( strcmp(argv[i], "-type") == 0 ) {
      if ( !VerifyReadSurfOpt("-type", 5, opt, argv) ) return XC_ERROR;

      if ( strncmp(argv[i+1], "gaus", 4) == 0 )  mols->type = MOLS_GAUSSIAN;
      else if ( strcmp(argv[i+1], "exp") == 0 )  mols->type = MOLS_EXP;
      else if ( strncmp(argv[i+1], "unig", 4) == 0 ) 
	mols->type = MOLS_UNIGAUSS;
      else if ( strncmp(argv[i+1], "unie", 4) == 0 ) 
	mols->type = MOLS_UNIEXP;
      else if ( strncmp(argv[i+1], "distf", 5) == 0 )
	mols->type = MOLS_DISTFUNC;
      else if ( strncmp(argv[i+1], "fs", 2) == 0 || 
		strncmp(argv[i+1], "fermis", 6) == 0 ) 
	mols->type = MOLS_FERMISURFACE;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"%s %s -type option\" command; option should be gauss, exp, unigauss, uniexp, distfuncor fermisurface\n", argv[i+1], argv[0], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }     
    /* -RENDER */
    else if ( strncmp(argv[i], "-rend", 5) == 0 ) {
      int index, rend;
      if ( !VerifyReadSurfOpt("-rend", 3, opt, argv) ) return TCL_ERROR;

      if ( Tcl_GetInt(interp, argv[i+1], &rend) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted integer, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }	
      index = FindSurfaceOfWindow( wc, mols );
      if ( index < 0 ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), "couldn't find surface #%d on window %s !!!\n",
		 mols->isosurf_index, argv[1]);  
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      wc->VPf.surface[index] = 1;

    }     
    /* -FS */
    else if ( strncmp(argv[i], "-fs", 3) == 0 ) {
      int argcList;
      const char **argvList;
      if ( !VerifyReadSurfOpt("-fs", 3, opt, argv) ) return TCL_ERROR;

      Tcl_SplitList(interp, argv[i+1], &argcList, &argvList);
      if ( !ReadFS_ARGS(interp, argcList, argvList, wc, mols) ) 
	return TCL_ERROR;
    }
    /* -RADIUS */
    else if ( strncmp(argv[i], "-radi", 5) == 0 ) {
      if ( !VerifyReadSurfOpt("-radi", 5, opt, argv) ) return TCL_ERROR;

      if ( strcmp(argv[i+1], "cov") == 0 )  mols->rad = rcov;
      else if ( strcmp(argv[i+1], "VdW") == 0 ) mols->rad = rvdw;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"%s %s -radius option\" command; option should be cov or VdW\n",
		 argv[i+1], argv[0], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -LEVEL */
    else if ( strncmp(argv[i], "-leve", 5) == 0 ) {
      double level;
      if ( !VerifyReadSurfOpt("-leve", 5, opt, argv) ) return TCL_ERROR;

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
      if ( !VerifyReadSurfOpt("-cuto", 5, opt, argv) ) return TCL_ERROR;

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
      if ( !VerifyReadSurfOpt("-colo", 5, opt, argv) ) return TCL_ERROR;

      if ( strncmp(argv[i+1], "atom", 4) == 0 )
	mols->colorscheme = MOLS_ATOMIC;
      else if ( strncmp(argv[i+1], "mono", 4) == 0 )
	mols->colorscheme = MOLS_MONOCHROME;	
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"cry_molsurf %s -colorscheme option\" command; option should be atomic or monocolor\n",	argv[i+1], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -DRAWSTYLE */
    else if ( strncmp(argv[i], "-draw", 5) == 0 ) {
      if ( !VerifyReadSurfOpt("-draw", 5, opt, argv) ) return TCL_ERROR;

      if ( strcmp(argv[i+1], "solid") == 0 ) 
	mols->dispt.drawstyle = ISOSURF_SOLID;
      else if ( strcmp(argv[i+1], "wire") == 0 ) 
	mols->dispt.drawstyle = ISOSURF_WIRE;
      else if ( strcmp(argv[i+1], "dot") == 0 ) 
	mols->dispt.drawstyle = ISOSURF_DOT;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"%s %s -drawstyle option\" command; option should be solid, wire or dot\n",
		 argv[i+1], argv[0], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -TRANSPARENT */
    else if ( strncmp(argv[i], "-tran", 5) == 0 ) {
      int transp;
      if ( !VerifyReadSurfOpt("-tran", 5, opt, argv) ) return TCL_ERROR;

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
      if ( !VerifyReadSurfOpt("-shade", 6, opt, argv) ) return TCL_ERROR;

      if ( strncmp(argv[i+1], "smooth", 6) == 0 ) 
	mols->dispt.shademodel = GL_SMOOTH;
      else if ( strncmp(argv[i+1], "flat", 4) == 0 ) 
	mols->dispt.shademodel = GL_FLAT;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"%s %s -shademodel option\" command; option should be smooth or flat\n",
		 argv[i+1], argv[0], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }
    /* -MONOCOLOR */
    else if ( strncmp(argv[i], "-mono", 5) == 0 ) {
      double rgba;
      int argcList;
      const char **argvList;
      if ( !VerifyReadSurfOpt("-mono", 5, opt, argv) ) return TCL_ERROR;

      Tcl_SplitList(interp, argv[i+1], &argcList, &argvList);      
      for (j=0; j<argcList; j++) {
	if ( Tcl_GetDouble(interp, argvList[j], &rgba) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	mols->monocolor[j]      = (float) rgba;
	mols->back_monocolor[j] = 1.0 - (float) rgba;
      }
      if ( j < 4) {
	mols->monocolor[3]      = 0.5;
	mols->back_monocolor[3] = 0.5;
      }
    }

    /* -FRONTMONOCOLOR */
    else if ( strncmp(argv[i], "-frontmonocolor", 10) == 0 ) {
      double rgba;
      int argcList;
      const char **argvList;
      if ( !VerifyReadSurfOpt("-frontmonocolor", 10, opt, argv) ) return TCL_ERROR;

      Tcl_SplitList(interp, argv[i+1], &argcList, &argvList);      
      for (j=0; j<argcList; j++) {
	if ( Tcl_GetDouble(interp, argvList[j], &rgba) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"wanted double, but got \"%s\"", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	mols->monocolor[j] = (float) rgba;
      }
      if ( j < 4 ) mols->monocolor[3] = 0.5;
    }

    /* -BACKMONOCOLOR */
    else if ( strncmp(argv[i], "-backmonocolor", 9) == 0 ) {
      double rgba;
      int    argcList;
      const char   **argvList;
      if ( !VerifyReadSurfOpt("-backmonocolor", 9, opt, argv) ) return TCL_ERROR;

      Tcl_SplitList(interp, argv[i+1], &argcList, &argvList);      
      for (j=0; j<argcList; j++) {
	if ( Tcl_GetDouble(interp, argvList[j], &rgba) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss),"wanted double, but got \"%s\"", argv[i+1]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
	mols->back_monocolor[j] = (float) rgba;
      }
      if ( j < 4 ) mols->back_monocolor[3] = 0.5;
    }

    /* -SURFACETYPE */
    else if ( strncmp(argv[i], "-surf", 5) == 0 ) {
      if ( !VerifyReadSurfOpt("-surf", 5, opt, argv) ) return TCL_ERROR;

      if ( strncmp(argv[i+1], "mols", 4) == 0 ) 
	mols->surftype = MOLS_MOLSURF;
      else if ( strncmp(argv[i+1], "gap", 3) == 0 ) 
	mols->surftype = MOLS_GAP;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown option \"%s\", while executing \"%s %s -surftype option\" command; option should be molsurf or gap\n",
		 argv[i+1], argv[0], argv[1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
    }

    /* -RESOLUTION */
    else if ( strncmp(argv[i], "-res", 4) == 0 ) {
      double res;
      if ( !VerifyReadSurfOpt("-res", 5, opt, argv) ) return TCL_ERROR;

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
      if ( !VerifyReadSurfOpt("-smooths", 8, opt, argv) ) return TCL_ERROR;

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
      if ( !VerifyReadSurfOpt("-smoothw", 8, opt, argv) ) return TCL_ERROR;

      if ( Tcl_GetDouble(interp, argv[i+1], &weight) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted double, but got \"%s\", while executing %s %s %s %s ...", argv[i+1], argv[0], argv[1], argv[2], argv[3]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);

	return TCL_ERROR;
      }
      mols->smooth_weight = (float) weight;
    }

    /* -FRONTFACE */
    else if ( strncmp(argv[i], "-frontface", 7) == 0 ) {      
      if ( strcmp(argv[i+1], "CW") == 0 ) 
	mols->frontface = GL_CW;
      else if ( strcmp(argv[i+1], "CCW") == 0 ) 
	mols->frontface = GL_CCW;
    }

    /* -REVERTNORMALS */
    else if ( strncmp(argv[i], "-revertnormals", 8) == 0 ) {
      int revertnormals;
      if ( Tcl_GetInt(interp, argv[i+1], &revertnormals) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted int, but got \"%s\"", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      mols->revertnormals = revertnormals;
    }

    else {
      char rss[1024];
      snprintf(rss, sizeof(rss),
	       "Usage: %s <tolg> -ident indent options ...", argv[0]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return TCL_ERROR;
    }
  }
  return TCL_OK;
}


int 
ReadFS_ARGS(Tcl_Interp *interp,	int argc, const char *argv[],
	    NEW_WIN_CONTEXT *wc, MOL_SURF *mols)
{
  register int    i;
  int             gridindex = -1, gridsubindex = -1, bandindex = -1;
  struct DATAGRID *grid;

  /* -fs { 
     -gridindex <index> \
     -gridsubindex <subindex> \
     -bandindex <index> \
     -celltype    para|bz \
     -cropbz    0|1  \    
     -displaycell 0|1 \		    ?  
     -celldisplaytype wire|rod|solid  ?
     }
  */
  for (i=0; i<argc; i+=2) {
    /* -GRIDINDEX */
    if ( strncmp(argv[i], "-gridi", 6) == 0 ) {
      if ( Tcl_GetInt(interp, argv[i+1], &gridindex) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted integer for -gridindex but got %s, while parsing -fs option !!!\n", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return XC_ERROR;
      } 
    }
    /* -GRIDSUBINDEX */
    else if ( strncmp(argv[i], "-grids", 6) == 0 ) {
      if ( Tcl_GetInt(interp, argv[i+1], &gridsubindex) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted integer for -gridsubindex but got %s, while parsing -fs option !!!\n", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return XC_ERROR;
      } 
    }
    /* -BANDINDEX */
    else if ( strncmp(argv[i], "-bandi", 6) == 0 ) {
      if ( Tcl_GetInt(interp, argv[i+1], &bandindex) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted integer for -bandindex but got %s, while parsing -fs option !!!\n", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return XC_ERROR;
      }
    }
    /* -CELLTYPE */
    else if ( strncmp(argv[i], "-cellt", 6) == 0 ) {
      if ( strncmp(argv[i+1], "para", 4) == 0 )
	mols->fs.celltype = XCR_PARAPIPEDAL;
      else if ( strncmp(argv[i+1], "bz", 2) == 0 )
	mols->fs.celltype = XCR_BZ;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"invalid cell type %s, must be para or bz !!!\n", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return XC_ERROR;
      }
    }
    /* -CROPBZ */
    else if ( strncmp(argv[i], "-cropb", 6) == 0 ) {
      int cropbz;
      if ( Tcl_GetInt(interp, argv[i+1], &cropbz) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted 0 or 1 for -cropbz but got %s, while parsing -fs option !!!\n", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return XC_ERROR;
      }
      mols->fs.cropbz = (short int) cropbz;
    }
    /* -DISPLAYCELL */
    else if ( strncmp(argv[i], "-displayc", 9) == 0 ) {
      if ( Tcl_GetInt(interp, argv[i+1], &wc->VPf.dispLattice) == TCL_ERROR ) {
	char rss[1024];
	snprintf(rss, sizeof(rss),"wanted 0 or 1 for -displaycell but got %s, while parsing -fs option !!!\n", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return XC_ERROR;
      }
    }
    /* -CELLDISPLAYTYPE */
    else if ( strncmp(argv[i], "-celld", 6) == 0 ) {
      if ( strncmp(argv[i+1], "w", 1) == 0 )
	wc->VPf.dispLatType = CELL_WIRE;
      else if ( strncmp(argv[i+1], "r", 1) == 0 )
	wc->VPf.dispLatType = CELL_ROD;
      else if ( strncmp(argv[i+1], "solidw", 6 ) == 0 )
	wc->VPf.dispLatType = CELL_SOLID_AND_WIRE;
      else if ( strncmp(argv[i+1], "solidr", 6 ) == 0 )
	wc->VPf.dispLatType = CELL_SOLID_AND_ROD;
      else if ( strncmp(argv[i+1], "s", 1) == 0 )
	wc->VPf.dispLatType = CELL_SOLID;
      else {
	char rss[1024];
	snprintf(rss, sizeof(rss),"unknown celldisplaytype %s, must be one of wire, rod, solid, solidwire or solidrod !!!\n", argv[i+1]);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return XC_ERROR;
      }
    }
    /* -INTERPOLATIONDEGREE */
    else if ( strncmp(argv[i], "-interpolationdegree", 15) == 0 ) {
      const char **argvList;
      int  argcList;
      int  ia, n[3];

      Tcl_SplitList(interp, argv[i+1], &argcList, &argvList);
      if ( argcList != 3 ) {
	char rss[1024];
	snprintf(rss, sizeof(rss), 
		 "wanted three elements but got %d elements", argcList);
	Tcl_SetResult(interp, rss, TCL_VOLATILE);
	return TCL_ERROR;
      }
      for (ia=0; ia<3; ia++) {
	if ( Tcl_GetInt(interp, argvList[ia], &(n[ia])) == TCL_ERROR ) {
	  char rss[1024];
	  snprintf(rss, sizeof(rss), "wanted integer but got %s", argvList[ia]);
	  Tcl_SetResult(interp, rss, TCL_VOLATILE);
	  return TCL_ERROR;
	}
      }

      mols->interp_n[0] = n[0];
      mols->interp_n[1] = n[1];
      mols->interp_n[2] = n[2];
    }

    /* -WIRECELLCOLOR */
    else if ( strncmp(argv[i], "-wirecellcolor", 10) == 0 ) {
      GetGlParam color;
      if ( !xcSplitList( XC_GET_RGBA, interp, argv + i + 1, &color ) )
	return TCL_ERROR;
      COPY_V (4, mols->fs.wirecellcolor, color.vec);
    }

    /* -SOLIDCELLCOLOR */
    else if ( strncmp(argv[i], "-solidcellcolor", 11) == 0 ) {
      GetGlParam color;
      if ( !xcSplitList( XC_GET_RGBA, interp, argv + i + 1, &color ) )
	return TCL_ERROR;
      COPY_V (4, mols->fs.solidcellcolor, color.vec);
    }

    else {
      char rss[1024];
      snprintf(rss, sizeof(rss), "unknown option %s, must be one of -gridindex, -bandindex, -celltype, -cropbz, -displaycell, -celldisplaytype, -wirecellcolor or -solidcellcolor", argv[i]);
      Tcl_SetResult(interp, rss, TCL_VOLATILE);
      return XC_ERROR;
    }
  }

  /**************************************************/
  /* mandatory options are -gridindex and -bandindex */
  /**************************************************/
  if ( bandindex < 0 || gridindex < 0 || gridsubindex < 0 ) {
    Tcl_SetResult(interp, "-gridindex, -gridsubindex and -bandindex options must be specified within -fs { options } !!!\n", TCL_STATIC);
    return XC_ERROR;
  }

  if ( (grid = FindDataGrid( gridindex )) == NULL ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), "can't find grid %d !!!\n", gridindex);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return XC_ERROR;
  } else mols->fs.gridindex = gridindex;

  if ( gridsubindex < grid->n_of_subgrids)
    mols->fs.gridsubindex = gridsubindex;
  else {
    char rss[1024];
    snprintf(rss, sizeof(rss), "to large gridsubindex %d, must be lower than %d !!!\n",
	     gridsubindex, grid->n_of_subgrids);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return XC_ERROR;
  }
  if ( bandindex < grid->nband ) 
    mols->fs.bandindex = bandindex;
  else {
    char rss[1024];
    snprintf(rss, sizeof(rss), "to large bandindex %d, must be lower than %d !!!\n",
	     gridsubindex, grid->nband);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return XC_ERROR;
  }
  return XC_OK;
}

