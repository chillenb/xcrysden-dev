#include <stdlib.h>
#include <GL/gl.h>

#define MINTOL   1.0e-15
#define ABS(x)   ( (x)>0 ? (x) : -(x) )

extern GLuint xcGenLists( GLsizei i );

/*
  The three routines below (findList1, compareParams1, and
  makeModelPtr1) are slightly modified versions of libaux findList,
  compareParams, and makeModelPtr routines (from "OpenGL Programming
  Guide", file: shapes.c). Below is the corresponding copyright.

  * (c) Copyright 1993, Silicon Graphics, Inc.
  * ALL RIGHTS RESERVED 
  * Permission to use, copy, modify, and distribute this software for 
  * any purpose and without fee is hereby granted, provided that the above
  * copyright notice appear in all copies and that both the copyright notice
  * and this permission notice appear in supporting documentation, and that 
  * the name of Silicon Graphics, Inc. not be used in advertising
  * or publicity pertaining to distribution of the software without specific,
  * written prior permission. 
  *
  * THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
  * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
  * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
  * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
  * GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
  * SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
  * KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
  * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
  * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
  * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
  * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
  * POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
  * 
  * US Government Users Restricted Rights 
  * Use, duplication, or disclosure by the Government is subject to
  * restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
  * (c)(1)(ii) of the Rights in Technical Data and Computer Software
  * clause at DFARS 252.227-7013 and/or in similar or successor
  * clauses in the FAR or the DOD or NASA FAR Supplement.
  * Unpublished-- rights reserved under the copyright laws of the
  * United States.  Contractor/manufacturer is Silicon Graphics,
  * Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
  *
  * OpenGL(TM) is a trademark of Silicon Graphics, Inc.

*/

GLuint findList1 (int lindex, GLdouble *paramArray, int size); 
static int compareParams1 (GLdouble *oneArray, GLdouble *twoArray, int size); 
GLuint makeModelPtr1 (int lindex, GLdouble *sizeArray, int count); 

/*****************************************************************************/
/* THIS IS FOR LISTS                                                         */
/*	structure for each geometric object	 */
typedef struct model { 
  GLuint       list;	        /*  display list to render object   */ 
  struct model *ptr;	        /*  pointer to next object	*/ 
  int          numParam;	/*  # of parameters		*/ 
  GLdouble     *params;	        /*  array with parameters	*/ 
} MODEL, *MODELPTR; 

/*	array of linked lists--used to keep track of display lists  
 *	for each different type of geometric object. 
 */
static MODELPTR lists[25] = { 
  NULL, NULL, NULL, NULL, NULL, 
  NULL, NULL, NULL, NULL, NULL, 
  NULL, NULL, NULL, NULL, NULL, 
  NULL, NULL, NULL, NULL, NULL, 
  NULL, NULL, NULL, NULL, NULL 
}; 


/****************NEXT THREE FUNCTIONS ARE LISTS MANAGEMENT *******************/
/*	linked lists--display lists for each different  
 *	type of geometric objects.  The linked list is  
 *	searched, until an object of the requested 
 *	size is found.  If no geometric object of that size 
 *	has been previously made, a new one is created. 
 */ 
GLuint  
findList1 (int lindex, GLdouble *paramArray, int size)  
{ 
    MODELPTR endList; 

    endList = lists[lindex]; 
    while (endList != NULL) { 
	if (compareParams1 (endList->params, paramArray, size)) 
	    return (endList->list);  
	endList = endList->ptr; 
    } 
    /*  if not found, return 0 and calling routine should 
     *  make a new list	 
     */ 
    return (0); 
} 


static int  
compareParams1 (GLdouble *oneArray, GLdouble *twoArray, int size)  
{ 
    int i; 
    int matches = 1; 
    GLdouble diff;

    for (i = 0; (i < size) && matches; i++) { 
      diff = *oneArray++ - *twoArray++;
      if ( ABS(diff) > MINTOL ) 
	matches = 0; 
    }
    return (matches); 
} 


GLuint  
makeModelPtr1 (int lindex, GLdouble sizeArray[], int count) 
{ 
  int i;
  MODELPTR newModel; 
  GLdouble *size;

  size = (GLdouble *) malloc( sizeof(GLdouble) * count );    
  for (i=0; i<count; i++) 
    size[i] = sizeArray[i];
  
  newModel = (MODELPTR) malloc (sizeof (MODEL)); 
  newModel->list     = xcGenLists (1); 
  newModel->numParam = count; 
  newModel->params   = size; 
  newModel->ptr      = lists[lindex]; 
  lists[lindex]      = newModel; 
  
  return (newModel->list); 
}
