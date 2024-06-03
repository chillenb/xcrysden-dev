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
 * Source: $XCRYSDEN_TOPDIR/C/xcdebug.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/
#include <stdio.h>
/* xcdebug.c */
void xcdebug(const char *text);
void xcErrDebug(const char *text);
void breakpoint(const char *text);

void 
xcdebug(const char *text)
{
  printf("%s\n",text);
  fflush(stdout);
}

void xcErrDebug(const char *text) {
  fprintf(stderr,"%s\n",text);
  fflush(stderr);
}

void breakpoint(const char *text) { xcdebug( text ); fflush(stdout); }
