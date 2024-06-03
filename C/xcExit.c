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
 * Source: $XCRYSDEN_TOPDIR/C/xcExit.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/


/* XC_ExitCmd - implementation of xc_exit custom Tcl command
 *  
 * we need a custom EXIT command, because if we have a structure opened, 
 * we must first release a Mesa context, delete lists, free variables and only
 * then we cam exit
 */
int
XC_ExitCmd(ClientData clientData,Tcl_Interp *interp,
		  int argc, char *argv[])
{

  if (argc != 1) {
   Tcl_SetResult(interp, "Usage: xc_exit", TCL_STATIC);
    return TCL_ERROR;
  }
 
  printf("XC_ExitCmd; VPf.stropened=%d\n",VPf.stropened);
  fflush(stdout);
  if (VPf.stropened) {
    /* maybe there are some 3D lists; destroy them */
    /* xcMaybeDestroyLists(); */
    xcDisplayFunc( (void *)NULL );  
    /* destiry mesa context */
    xcDestroyMesaContext();  
    /* now FREE all variables */
    FreeAllVariables();
  }
  
  /* now we can exit normally*/
  Tcl_Exit(0);
  return TCL_OK;
}
