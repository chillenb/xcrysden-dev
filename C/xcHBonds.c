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
* Source: $XCRYSDEN_TOPDIR/C/xcHBonds.c
* ------                                                                    *
* Copyright (c) 1996-2003 by Anton Kokalj                                   *
*****************************************************************************

*/

#include <togl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "vector.h"
#include "xcfunc.h"
#include "memory.h"

extern H_Bond hbonds;

static int _Hbond_MakeList(Tcl_Interp *interp, const char *argvI1, int **list);
static void _Hbond_PrintList(Tcl_Interp *interp, int *list);

static void Usage(Tcl_Interp *interp) {
  Tcl_SetResult(interp, "Usage: xc_hbonds <toglName> OPTIONS\nwhere options are:\n    on|off\nor\n    H_like_list  {H_like_list atomic numbers}\n    O_like_list  {O_like_list atomic numbers}\n    length_min  <min-length>\nor\nlength_max  <max-length>", TCL_STATIC);
}

/*****************************************************************************
 * xc_hbonds <toglName> on
 * or                   off
 * or                   get what
 * or                   set what <what_values> ?what <value? ...
 *
 *                      what:        what_value:
 * ---------------------------------------------------------------------------
 *                      -H_like_list  {H_like_list atomic numbers} 
 *                      -O_like_list  {O_like_list atomic numbers}
 *                      -length_min   double-number
 *                      -length_max   double-number
 *                      -angle_min    angle-min
 *
 *****************************************************************************/
int 
XC_HBondsCmd(ClientData clientData,Tcl_Interp *interp,
	     int argc, const char *argv[])
{
  struct Togl *togl;
  char         c;
  int          i, len, inum;
  double       num;

  if ( argc < 3 ) {
    Usage(interp);
    return TCL_ERROR;
  }

  /* 
   * find togl associated with toglName 
   */
  if ( Togl_GetToglFromName(interp, argv[1], &togl) == TCL_ERROR ) {
    char rss[1024];
    snprintf(rss, sizeof(rss), 
	     "couldn't find %s togl widget", argv[1]);
    Tcl_SetResult(interp, rss, TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* 
   * parse options 
   */

  c=argv[2][0];
  len=strlen(argv[2]);

  if ((c == 'o') && (strncmp(argv[2], "off", len) == 0) && (len >= 2)) {
    /* ON */
    VPf.Hbond = GL_FALSE;
  }
  else if ((c == 'o') && (strncmp(argv[2], "on", len) == 0) && (len >= 2)) {
    /* OFF */
    VPf.Hbond = GL_TRUE;
  } else if ((c == 'g') && (strncmp(argv[2], "get", len) == 0) && (len >= 1)) {
    /* 
       GET 
    */
    char string[20];
    if ( strncmp(argv[3], "-H_", 3) == 0 ) 
      {
	/* -H_like_list */
	_Hbond_PrintList(interp, hbonds.H_like_list);
      } 
    else if ( strncmp(argv[3], "-O_", 3) == 0 ) 
      {
	/* -O_like_list */
	_Hbond_PrintList(interp, hbonds.O_like_list);
      } 
    else if ( strcmp(argv[3], "-length_min") == 0 ) 
      {
	sprintf(string,"%15.10f ", hbonds.length_min);
	Tcl_AppendResult(interp, string, (char*)NULL);
      } 
    else if ( strcmp(argv[3], "-length_max") == 0 ) 
      {
	sprintf(string,"%15.10f ", hbonds.length_max);
	Tcl_AppendResult(interp, string, (char*)NULL);
      }
    else if ( strcmp(argv[3], "-angle_min") == 0 ) 
      {
	sprintf(string,"%15.10f ", hbonds.angle_min);
	Tcl_AppendResult(interp, string, (char*)NULL);
      }
    else if ( strcmp(argv[3], "-line_width") == 0 ) 
      {
	sprintf(string,"%15.10f ", hbonds.line_width);
	Tcl_AppendResult(interp, string, (char*)NULL);
      }
    else if ( strcmp(argv[3], "-line_pattern") == 0 ) 
      {
	sprintf(string,"0x%04x", hbonds.line_pattern);
	Tcl_AppendResult(interp, string, (char*)NULL);
      }
    else if ( strcmp(argv[3], "-line_patternsize") == 0 ) 
      {
	sprintf(string,"%d", hbonds.line_patternsize);
	Tcl_AppendResult(interp, string, (char*)NULL);
      }
    else if ( strcmp(argv[3], "-color") == 0 ) 
      {
	sprintf(string,"%5.3f %5.3f %5.3f", 
		hbonds.color[0], hbonds.color[1], hbonds.color[2]);
	Tcl_AppendResult(interp, string, (char*)NULL);
      }
    else 
      {
	Usage(interp);
	return TCL_ERROR;
      }
  } else if ((c == 's') && (strncmp(argv[2], "set", len) == 0) && (len >= 1)) {
    /* 
       SET 
    */
    for (i=3; i<argc; i+=2) {
      if ( strncmp(argv[i], "-H_", 3) == 0 ) 
	{
	  /* -H_like_list */
	  if ( _Hbond_MakeList(interp, argv[i+1], &(hbonds.H_like_list)) == TCL_ERROR ) {
	    char rss[1024];
	    snprintf(rss, sizeof(rss), "error parsing H_like_list: %s\n", argv[i+1]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	} 
      else if ( strncmp(argv[i], "-O_", 3) == 0 ) 
	{
	  /* -O_like_list */
	  if ( _Hbond_MakeList(interp, argv[i+1], &(hbonds.O_like_list)) == TCL_ERROR ) {
	    char rss[1024];
	    snprintf(rss, sizeof(rss), "error parsing O_like_list: %s\n", argv[i+1]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	} 
      else if ( strcmp(argv[i], "-length_min") == 0 ) 
	{
	  /* -length_min */
	  if ( Tcl_GetDouble(interp, argv[i+1], &num) == TCL_ERROR ) {
	    char rss[1024];
	    snprintf(rss, sizeof(rss),"wanted double, but got %s\n", argv[i+1]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	  hbonds.length_min = num;			   
	} 
      else if ( strcmp(argv[i], "-length_max") == 0 ) 
	{
	  /* -length_max */
	  if ( Tcl_GetDouble(interp, argv[i+1], &num) == TCL_ERROR ) {
	    char rss[1024];
	    snprintf(rss, sizeof(rss),"wanted double, but got %s\n", argv[i+1]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	  hbonds.length_max = num;			   	
	}
      else if ( strcmp(argv[i], "-angle_min") == 0 ) 
	{
	  /* -angle_min */
	  if ( Tcl_GetDouble(interp, argv[i+1], &num) == TCL_ERROR ) {
	    char rss[1024];
	    snprintf(rss, sizeof(rss),"wanted double, but got %s\n", argv[i+1]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	  hbonds.angle_min = num;			   
	} 
      else if ( strcmp(argv[i], "-line_width") == 0 ) 
	{
	  /* -line_width */
	  if ( Tcl_GetDouble(interp, argv[i+1], &num) == TCL_ERROR ) {
	    char rss[1024];
	    snprintf(rss, sizeof(rss),"wanted double, but got %s\n", argv[i+1]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	  hbonds.line_width = (float) num;			   
	} 
      else if ( strcmp(argv[i], "-line_pattern") == 0 ) 
	{
	  /* -line_pattern */
	  if ( Tcl_GetInt(interp, argv[i+1], &inum) == TCL_ERROR ) {
	    char rss[1024];
	    snprintf(rss, sizeof(rss),"wanted integer, but got %s\n", argv[i+1]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	  hbonds.line_pattern = (unsigned short) inum;			   
	}
      else if ( strcmp(argv[i], "-line_patternsize") == 0 ) 
	{
	  /* -line_patternsize */
	  if ( Tcl_GetInt(interp, argv[i+1], &inum) == TCL_ERROR ) {
	    char rss[1024];
	    snprintf(rss, sizeof(rss),"wanted integer, but got %s\n", argv[i+1]);
	    Tcl_SetResult(interp, rss, TCL_VOLATILE);
	    return TCL_ERROR;
	  }
	  hbonds.line_patternsize = inum;			   
	}
      else if ( strcmp(argv[i], "-color") == 0 ) 
	{
	  /* -color */
	  GetGlParam color;
	  if ( !xcSplitList( XC_GET_RGB, interp, argv + i + 1, &color ) )
	    return TCL_ERROR;
	  COPY_V (3, hbonds.color, color.vec);
	}
      else 
	{
	  Usage(interp);
	  return TCL_ERROR;
	}
    }
    /*
      make Hbonds
    */
    make_H_Bonds();
  } else {
    /* 
       !!! wrong usage !!!
    */
    Usage(interp);
    return TCL_ERROR;
  }

  return TCL_OK;
}


static int _Hbond_MakeList(Tcl_Interp *interp, const char *argvI1, int **list) {  
  int argcList, j, intnum;
  const char **argvList; 
  int size = 0;
  int *listPtr;

  /* count the list-size */
  listPtr = *list; 
  while (*listPtr > 0) {
    listPtr++; size++;      
  }

  /* split the list of argvI1 */
  Tcl_SplitList(interp, argvI1, &argcList, &argvList);

  /* if size < argcList+1 ==> reallocate the list */
  if ( size < argcList+1 ) 
    *list = xcRealloc(*list, (size_t) (argcList+1)*sizeof(int));

  /* make the list */
  listPtr = *list;
  for (j=0; j<argcList; j++) 
    {
      if ( Tcl_GetInt(interp, argvList[j], &intnum) == TCL_ERROR ) return TCL_ERROR;
      listPtr[j] = intnum;
    }
  listPtr[argcList] = 0; /* mark the end-of-list with ZERO "0" */

  return TCL_OK;
}


static void _Hbond_PrintList(Tcl_Interp *interp, int *list) {
  char string[5];  
  while (*list > 0) {
    sprintf(string,"%3d ",*list);
    Tcl_AppendResult(interp, string, (char*)NULL);
    list++;
  }
}
