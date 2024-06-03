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
 * Source: $XCRYSDEN_TOPDIR/C/signal.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

void xcFPEHandler(void);

#if defined _LINUXALPHA || defined _DEC_CUSTOM_FPE_HANDLER

#include <stdio.h>
#include <signal.h>

typedef void (*sighandler_t)(int);
static sighandler_t FPEHandler(int signum);

void 
xcFPEHandler(void) {
  sighandler_t FPEHandler(int signum);
  signal(SIGFPE, (sighandler_t) FPEHandler);
}


static
sighandler_t FPEHandler(int signum) {
  static int first_time = 1;

  if (first_time) {
    fprintf(stderr, 
	    "Warning: Floating Point Exception occurred - Ignoring !!!\n", NULL);
    first_time = 0;
  }
    
  return SIG_IGN;
}
#else

void 
xcFPEHandler(void) {
  ;
}

#endif
