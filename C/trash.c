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
 * Source: $XCRYSDEN_TOPDIR/C/trash.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

/* ========================================================================= */
/*              HERE ARE ROUTINES THAT ARE NO MORE NEEDED                    */
/* ========================================================================= */

/* preparing for 2D; */
/* going from "AtomBond3D  coor3D" -> "AtomBond  coor";
 * actually we are not going from coor3D -> coor, but just
 * update coor on the basis of MajorMat !!!!!
 */
void 
From3Dto2D(void)
{
  int i, atbn = natoms + nbonds;
  double x,y,z;

  for(i=1; i<=atbn; i++)
    {
      UpdateCoor( &(coor + i)->x1, &(coor + i)->y1, &(coor + i)->z1 );
      if ( (coor + i)->flag == BOND ) {
	UpdateCoor( &(coor + i)->x2, &(coor + i)->y2, &(coor + i)->z2 );
	/* if ATOM and BOND are on the same plane, ATOM must be drawn 
	 first so we will add very little to z-BOND */
	*(zorient + i)  = 0.5 * ((coor + i)->z1 + (coor + i)->z2) + 0.0001;
      }
    }  
}



void
From2Dto3D(void)
{
  int i;
  int atbn;

  printf("From2Dto3D!!!!\n",NULL);
  fflush(stdout);

  /* delete previous made 3D lists, if I do not do this than the rest of this
   * function is without any sense, because later I use 3D Lists for 
   * rendering, and if I do not use updated 3DLists ...
   */
  xcMaybeDelete3DLists();
  atbn = natoms + nbonds;
  for(i=1; i <= atbn; i++)
    {
      (coor3D + i)->flag = (coor + i)->flag;
      (coor3D + i)->nat = (coor + i)->nat;
      /* (coor3D + i)->list is not copied !!!! */
      (coor3D + i)->x1 = (coor + i)->x1;
      (coor3D + i)->y1 = (coor + i)->y1;
      (coor3D + i)->z1 = (coor + i)->z1;
      if ( (coor3D + i)->flag == BOND )
	GetCylinderPar(i, (coor + i)->x2 - (coor + i)->x1,
		          (coor + i)->y2 - (coor + i)->y1,
		          (coor + i)->z2 - (coor + i)->z1 );
    }
  /* now assign NewMat */
  LoadIdentity( MajorMat );
  MajorMatToVec();
}


void
xcClearScreen(void)
{
  /* make a black polygon on whole window -> so structure will disapear */
  glClear(GL_COLOR_BUFFER_BIT);
  glColor3f( 0.0, 0.0, 0.0);
  glBegin(GL_POLYGON);
    glVertex2d(-MVf.structsize, MVf.structsize);
    glVertex2d(MVf.structsize, MVf.structsize);
    glVertex2d(MVf.structsize, -MVf.structsize);
    glVertex2d(-MVf.structsize, -MVf.structsize);
  glEnd();
  glFlush();
  Togl_SwapBuffers(togl);
}


