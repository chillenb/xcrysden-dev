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
 * Source: $XCRYSDEN_TOPDIR/C/xcDisplayFunc2.c
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#include <togl.h>
#include <stdio.h>
#include <string.h>
#include "struct.h"
#include "xcfunc.h"

extern GLuint         PointList;
extern Options3D      is;
extern RasterFontSize rf;
extern OrthoProj      ort; 
extern PERSPECTIVE    persp;
extern AtomicLabel    *atomLabel, globalAtomLabel;
extern short          *do_not_display_atomlabel;
extern char           *element[];
extern GLuint         fontOffset;
extern realTimeMove makeMovie;

/* --- function prototypes --- */
void xcAssignDisplayFunc(const char *dispMode);
void RewriteCoor(GLenum type);
void xcUpdateCoor( double mat[][4], double *x, double *y, double *z);

void xcWireFrameL2D(struct Togl *togl);
void xcWireFrameF2D(struct Togl *togl);
void xcWireFrameFL2D(struct Togl *togl);

void xcPointLineL2D(struct Togl *togl);
void xcPointLineF2D(struct Togl *togl);
void xcPointLineFL2D(struct Togl *togl);

void xcBallStick1L2D(struct Togl *togl);
void xcBallStick1F2D(struct Togl *togl);
void xcBallStick1FL2D(struct Togl *togl);

void xcBallStick2L2D(struct Togl *togl);
void xcBallStick2F2D(struct Togl *togl);
void xcBallStick2FL2D(struct Togl *togl);

void xcTurnLabelsOn(void);
void xcTurnLabelsOff(void);
void xcTurnFramesOn(void);
void xcTurnFramesOff(void);
void WhichDispFunc(void);
void UpdateDispFunc(void);
void UpdateProjection(void);

/* ---- extern function prototypes ---- */
extern void xcBall2D(int atom);
extern void xcDisplayFunc( void (*Func)(struct Togl *togl) );
extern void xcWireFrame2D(struct Togl *togl);
extern void xcPointLine2D(struct Togl *togl);
extern void xcPipeBall2D(struct Togl *togl);
extern void xcBallStick12D(struct Togl *togl);
extern void xcBallStick22D(struct Togl *togl);
extern void xcSpaceFill2D(struct Togl *togl);
extern void (*xcDisplay)(struct Togl *togl);
extern void xcDisplayXYZ(void);
/*extern void xcMakeBall3DLists(void);*/
/*extern void xcMakeSpaceFill3DLists(void);*/
/*extern void xcMakeStick3DLists(void);*/
/*  extern void xcMakeFrame3DLists(void); */
/*  extern void xcMakeBallLabel3DList(void); */
/*  extern void xcMakeSpaceLabel3DList(void); */
/* extern void xcMakeBallLists(void); */
extern int MakeSticks1(int i,
		       GLdouble *x1, GLdouble *y1, GLdouble *z1, 
		       GLdouble *x2, GLdouble *y2, GLdouble *z2, 
		       GLdouble *x3, GLdouble *y3, GLdouble *z3, 
		       GLdouble *x4, GLdouble *y4, GLdouble *z4);

extern int MakeSticks2(int i, int col, int flag);
extern void xcMakeProjection2D(const char *mode);


void
xcAssignDisplayFunc(const char *dispMode)
{
  /* 2D or 3D */
  if ( dimType == XC_2D ) {
    if ( strcmp(dispMode,"WF") == 0 ) { /* WireFrame */

      /* only bonds are taken into account */
      
      if ( VPf.labelsOn == 0 && VPf.framesOn == 0 ) {
	RewriteCoor(XC_BOND);
	xcDisplayFunc(xcWireFrame2D);
	/* BECAUSE "tmp_nobjects" may differ from last time we must 
	   Z-orientate, to push tmp_nobjects which will be display to front */
	hpsort_index1(tmp_nobjects, zorient, iwksp); 
      }

      /* bonds + atoms (labels) are taken into account */
      if ( VPf.labelsOn == 1 && VPf.framesOn == 0 ) {
	/* Rewritecoor xat -> coor */
	RewriteCoor(XC_BOND | XC_LABEL);
	xcDisplayFunc(xcWireFrameL2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }

      /* bonds + frames are taken into account */
      if ( VPf.labelsOn == 0 && VPf.framesOn == 1 ) {
	/* Rewritecoor xat -> coor */
	RewriteCoor(XC_BOND | XC_FRAME);
	xcDisplayFunc(xcWireFrameF2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }

      /* bonds + atoms (labels) + frames are taken into account */
      if ( VPf.labelsOn == 1 && VPf.framesOn == 1 ) {
	/* Rewritecoor xat -> coor */
	RewriteCoor(XC_BOND | XC_LABEL | XC_FRAME);
	xcDisplayFunc(xcWireFrameFL2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }
    }

    else if ( strcmp(dispMode,"PL") == 0 ) { /* PointLine */

      /* only atoms and bonds are taken into account */
      if ( VPf.labelsOn == 0 && VPf.framesOn == 0 ) {
	RewriteCoor(XC_BOND | XC_ATOM); 
	xcDisplayFunc(xcPointLine2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }

      /* bonds + atoms + labels are taken into account */
      if ( VPf.labelsOn == 1 && VPf.framesOn == 0 ) {
	/* Rewritecoor xat -> coor */
	RewriteCoor(XC_BOND | XC_ATOM); /* XC_LABEL and XC_ATOM are synomyms */
	xcDisplayFunc(xcPointLineL2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }

      /* bonds + atoms + frames are taken into account */
      if ( VPf.labelsOn == 0 && VPf.framesOn == 1 ) {
	/* Rewritecoor xat -> coor */
	RewriteCoor(XC_BOND | XC_ATOM | XC_FRAME);
	xcDisplayFunc(xcPointLineF2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }

      /* bonds + atoms + labels + frames are taken into account */
      if ( VPf.labelsOn == 1 && VPf.framesOn == 1 ) {
	/* Rewritecoor xat -> coor */
	RewriteCoor(XC_BOND | XC_ATOM | XC_FRAME);
	xcDisplayFunc(xcPointLineFL2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }
    }

    else if ( strcmp(dispMode,"PB") == 0 ) { /* PointLine */
      GLuint flag = XC_BOND | XC_ATOM;
      if ( VPf.framesOn == 1 ) flag |= XC_FRAME;
      RewriteCoor (flag); 
      xcDisplayFunc(xcPipeBall2D);
      hpsort_index1(tmp_nobjects, zorient, iwksp);
    }

    else if ( strcmp(dispMode,"BS1") == 0 ) { /* PointLine */

      /* only atoms and bonds are taken into account */
      if ( VPf.labelsOn == 0 && VPf.framesOn == 0 ) {
	RewriteCoor(XC_BOND | XC_ATOM);
	xcDisplayFunc(xcBallStick12D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }

      /* bonds + atoms + labels are taken into account */
      if ( VPf.labelsOn == 1 && VPf.framesOn == 0 ) {
	/* Rewritecoor xat -> coor */
	RewriteCoor(XC_BOND | XC_ATOM);
	xcDisplayFunc(xcBallStick1L2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }

      /* bonds + atoms + frames are taken into account */
      if ( VPf.labelsOn == 0 && VPf.framesOn == 1 ) {
	/* Rewritecoor xat -> coor */
	RewriteCoor(XC_BOND | XC_ATOM | XC_FRAME);
	xcDisplayFunc(xcBallStick1F2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }

      /* bonds + atoms + labels + frames are taken into account */
      if ( VPf.labelsOn == 1 && VPf.framesOn == 1 ) {
	/* Rewritecoor xat -> coor */
	RewriteCoor(XC_BOND | XC_ATOM | XC_FRAME);
	xcDisplayFunc(xcBallStick1FL2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }
    }

    else if ( strcmp(dispMode,"BS2") == 0 ) { /* PointLine */

      /* only atoms and bonds are taken into account */
      if ( VPf.labelsOn == 0 && VPf.framesOn == 0 ) {
	RewriteCoor(XC_BOND | XC_ATOM);
	xcDisplayFunc(xcBallStick22D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }

      /* bonds + atoms + labels are taken into account */
      if ( VPf.labelsOn == 1 && VPf.framesOn == 0 ) {
	/* Rewritecoor xat -> coor */
	RewriteCoor(XC_BOND | XC_ATOM);
	xcDisplayFunc(xcBallStick2L2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }

      /* bonds + atoms + frames are taken into account */
      if ( VPf.labelsOn == 0 && VPf.framesOn == 1 ) {
	/* Rewritecoor xat -> coor */
	RewriteCoor(XC_BOND | XC_ATOM | XC_FRAME);
	xcDisplayFunc(xcBallStick2F2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }

      /* bonds + atoms + labels + frames are taken into account */
      if ( VPf.labelsOn == 1 && VPf.framesOn == 1 ) {
	/* Rewritecoor xat -> coor */
	RewriteCoor(XC_BOND | XC_ATOM | XC_FRAME);
	xcDisplayFunc(xcBallStick2FL2D);
	hpsort_index1(tmp_nobjects, zorient, iwksp);
      }
    }

    else if ( strcmp(dispMode,"SF") == 0 ) { /* SpaceFills */
      if ( VPf.framesOn == 0 ) {
	RewriteCoor(XC_ATOM); 
      } else {
	RewriteCoor(XC_ATOM | XC_FRAME);
      }
      xcDisplayFunc(xcSpaceFill2D);
      hpsort_index1(tmp_nobjects, zorient, iwksp);
    }

  }
}


void
RewriteCoor(GLenum type)
{
  int i;
  int nobj = 0;
  /* must be bonds taken into account */
  if ( type & XC_BOND ) {
    /* rewrite BONDS --- BONDS */
    for(i=1; i<=nbonds; i++) {      
      (coor + i)->flag = BOND;
      (coor + i)->nat= *(natbond + i);
      (coor + i)->sqn= *(sqnbond + i);
      (coor + i)->bondend= *(bondend + i);
      (coor + i)->x1 = *(xbond + i);
      (coor + i)->y1 = *(ybond + i);
      (coor + i)->z1 = *(zbond + i);
      (coor + i)->x2 = *(xbond2 + i);
      (coor + i)->y2 = *(ybond2 + i);
      (coor + i)->z2 = *(zbond2 + i);
      xcUpdateCoor( vec.crdmajor, &(coor + i)->x1, &(coor + i)->y1, \
		    &(coor + i)->z1 );
      xcUpdateCoor( vec.crdmajor, &(coor + i)->x2, &(coor + i)->y2, \
		    &(coor + i)->z2 );
      /* if ATOM and BOND are on the same plane, ATOM must be drawn 
	 first so we will add very little to z-BOND */
      *(zorient + i)  = 0.5 * ((coor + i)->z1 + (coor + i)->z2) + 
	Z_OFFSET((coor + i)->flag);
    }
    nobj += nbonds;    
  }

  if ( type & XC_ATOM ) {
    /* ATOMS --- ATOMS */
    int mmm = nobj + natoms;
    for(i = nobj + 1; i <= mmm; i++) {
      (coor + i)->flag = ATOM;	
      (coor + i)->nat = *(nat + i - nobj);
      (coor + i)->sqn = i - nobj;
      (coor + i)->x1 = *(xat + i - nobj);
      (coor + i)->y1 = *(yat + i - nobj);
      (coor + i)->z1 = *(zat + i - nobj);
      xcUpdateCoor( vec.crdmajor, &(coor + i)->x1, &(coor + i)->y1, \
		    &(coor + i)->z1 );
      *(zorient + i)  = (coor + i)->z1;
    }
    nobj += natoms;
  }
  if ( type & XC_FRAME ) {
    int mmm = nobj + nframes;
    /* FRAMES --- FRAMES */
    for(i = nobj + 1; i <= mmm; i++) {      
      (coor + i)->flag = FRAME;
      (coor + i)->nat = *(frametype + i - nobj); 
      (coor + i)->sqn = i - nobj;
      (coor + i)->x1 = *(xframe + i - nobj);
      (coor + i)->y1 = *(yframe + i - nobj);
      (coor + i)->z1 = *(zframe + i - nobj);
      (coor + i)->x2 = *(xframe2 + i - nobj);
      (coor + i)->y2 = *(yframe2 + i - nobj);
      (coor + i)->z2 = *(zframe2 + i -nobj);
      xcUpdateCoor( vec.crdmajor, &(coor + i)->x1, &(coor + i)->y1, \
		    &(coor + i)->z1 );
      xcUpdateCoor( vec.crdmajor, &(coor + i)->x2, &(coor + i)->y2, \
		    &(coor + i)->z2 );
      /* if ATOM and BOND and FRAME are on the same plane, ATOM must be drawn 
	 first, than BOND, so we will add very little to z-FRAME;
	 but a little more than for BOND*/
      *(zorient + i)  = 0.5 * ((coor + i)->z1 + (coor + i)->z2) + 
	Z_OFFSET((coor + i)->flag);
    }
     nobj += nframes;
  }

  tmp_nobjects = nobj;
  printf("END RewriteCoor; tmp_nobjects=%d",tmp_nobjects);
  fflush(stdout);
}


void
xcUpdateCoor( double mat[][4], double *x, double *y, double *z)
{
  double x1, y1, z1;

  x1 = mat[0][0] * *x + mat[1][0] * *y + mat[2][0] * *z;
  y1 = mat[0][1] * *x + mat[1][1] * *y + mat[2][1] * *z;
  z1 = mat[0][2] * *x + mat[1][2] * *y + mat[2][2] * *z;
  *x = x1;
  *y = y1;
  *z = z1;
}

void makeAtomLabel2D(int ith, GLdouble Xoffset, GLdouble Yoffset) {
  char     *label;
  GLdouble width, height, sizeF;
  int      ia;
  
  ia    = (coor + ith)->sqn;
  if (ia > natoms) {
    /* 
       this should not happen, but sometimes it does with
       selection-mode (it is a bug) 
    */
    return;
  }
  sizeF = (2.0 * VPf.VPfactor);

  if (atomLabel[ia].base) {
    label  = atomLabel[ia].label;
    if (atomLabel[ia].tkfont) {
      width = (GLdouble) ( Tk_TextWidth(atomLabel[ia].tkfont, label, strlen(label)) );
    } else {
      width  = (GLdouble) (atomLabel[ia].width * strlen(label));
    }
    height = (GLdouble) atomLabel[ia].height;
  } else {
    if (atomLabel[ia].label) label = atomLabel[ia].label;
    else label = element[*(nat + ia)];

    if (globalAtomLabel.base) {
      if (globalAtomLabel.tkfont) {
	width = (GLdouble) ( Tk_TextWidth(globalAtomLabel.tkfont, label, strlen(label)) );
      } else {
	width  = (GLdouble) globalAtomLabel.width;
      }
      height = (GLdouble) globalAtomLabel.height;      
    } else {
      width  = (GLdouble) (rf.wid * strlen(label)); 
      height = (GLdouble) rf.height; 
    }
  }
  
  if ( Xoffset != CENTERED_LABEL ) {
    if (atomLabel[ia].base) glColor3fv (atomLabel[ia].bright_color);
    else glColor3fv (globalAtomLabel.bright_color);
  } else {
    if (atomLabel[ia].base) glColor3fv (atomLabel[ia].dark_color);
    else glColor3fv (globalAtomLabel.dark_color);

    Xoffset = -width / sizeF;
    Yoffset = -height / sizeF;
  }
  
  glRasterPos3d( (coor + ith)->x1 + Xoffset, 
		 (coor + ith)->y1 + Yoffset,
		 (coor + ith)->z1 );  
  
  if ( !globalAtomLabel.base && !atomLabel[ia].base ) {
    /* 
       render default raster font 
    */
    
    if (atomLabel[ia].base) {
      /* custom label */
      glListBase(atomLabel[ia].base);
      if (atomLabel[ia].do_display && !do_not_display_atomlabel[ia])
	xcFont_PrintString (label);
    } else {
      /* default label */
      if (globalAtomLabel.do_display && !do_not_display_atomlabel[ia])
	glCallList( atomlabelOffset + *(nat + ia) );
    }    
  } else {
    /*
      render one of loaded Togl fonts 
    */
    if (atomLabel[ia].base) {
      glListBase (atomLabel[ia].base);
      if (atomLabel[ia].do_display && !do_not_display_atomlabel[ia])
	xcFont_PrintString( label );
    }
    else {
      glListBase (globalAtomLabel.base);
      if (globalAtomLabel.do_display && !do_not_display_atomlabel[ia])
	xcFont_PrintString( label );    
    }
  }
}


void
xcWireFrameL2D(struct Togl *togl)
{
  int i;
    
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glLineWidth( (GLfloat) VPf.WFlinewidth );
  for(i = 1; i <= tmp_nobjects; i++) {
    if ( (coor + *(iwksp + i))->flag == BOND )
      {
	/*printf("xcWireFrameL2D> loop BOND, i=%d",*(iwksp + i));
	fflush(stdout);*/
	glBegin(GL_LINES);
	  /* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */
	  if ( VPf.unibond) glColor3fv(unibondCol);
	  else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	  glVertex3d( (coor + *(iwksp + i))->x1, 
		      (coor + *(iwksp + i))->y1, 
		      (coor + *(iwksp + i))->z1 );
	  glVertex3d( (coor + *(iwksp + i))->x2, 
		      (coor + *(iwksp + i))->y2,
		      (coor + *(iwksp + i))->z2 );
	glEnd();
      }
    if ( (coor + *(iwksp + i))->flag == ATOM )
      {
	/* printf("xcWireFrameL2D> loop ATOM, i=%d",*(iwksp + i));
	fflush(stdout); */
	makeAtomLabel2D ( *(iwksp + i), rf.wp2, rf.hp2 );
      }	
  }

  /* check if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl); 

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void
xcWireFrameF2D(struct Togl *togl)
{
  int i;
    
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  for(i = 1; i <= tmp_nobjects; i++) {
    if ( (coor + *(iwksp + i))->flag == BOND )
      {
	glLineWidth( (GLfloat) VPf.WFlinewidth );  
	glBegin(GL_LINES);
	  /* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */
	  if ( VPf.unibond) glColor3fv(unibondCol);
	  else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	  glVertex3d( (coor + *(iwksp + i))->x1, 
		      (coor + *(iwksp + i))->y1,
		      (coor + *(iwksp + i))->z1 );
	  glVertex3d( (coor + *(iwksp + i))->x2, 
		      (coor + *(iwksp + i))->y2,
		      (coor + *(iwksp + i))->z2 );
	glEnd();
      }
    if ( (coor + *(iwksp + i))->flag == FRAME )
      {
	if ( (coor + *(iwksp + i))->nat == 2 ) {
	  glLineStipple(2, FRAME_PAT);
	  glEnable(GL_LINE_STIPPLE);
	}
	glLineWidth( (GLfloat) VPf.framewidth );  
	glBegin(GL_LINES);
	  glColor3fv( framecol );
	  glVertex3d( (coor + *(iwksp + i))->x1, 
		      (coor + *(iwksp + i))->y1,
		      (coor + *(iwksp + i))->z1 );
	  glVertex3d( (coor + *(iwksp + i))->x2, 
		      (coor + *(iwksp + i))->y2,
		      (coor + *(iwksp + i))->z2 );
	glEnd(); 
	glDisable(GL_LINE_STIPPLE);
      }	
  }
  
  /* check if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void
xcWireFrameFL2D(struct Togl *togl)
{
  int i;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  for(i = 1; i <= tmp_nobjects; i++) {
    if ( (coor + *(iwksp + i))->flag == BOND )
      {
	glLineWidth( (GLfloat) VPf.WFlinewidth );  
	glBegin(GL_LINES);
	  /* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */
	  if ( VPf.unibond) glColor3fv(unibondCol);
	  else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	  glVertex3d( (coor + *(iwksp + i))->x1, 
		      (coor + *(iwksp + i))->y1,
		      (coor + *(iwksp + i))->z1 );
	  glVertex3d( (coor + *(iwksp + i))->x2,
		      (coor + *(iwksp + i))->y2,
		      (coor + *(iwksp + i))->z2 );
	glEnd();
      }
    if ( (coor + *(iwksp + i))->flag == ATOM )
      {
	makeAtomLabel2D ( *(iwksp + i), rf.wp2, rf.hp2 );
      }
    if ( (coor + *(iwksp + i))->flag == FRAME )
      {
	if ( (coor + *(iwksp + i))->nat == 2 ) {
	  glLineStipple(2, FRAME_PAT);
	  glEnable(GL_LINE_STIPPLE);
	}
	glLineWidth( (GLfloat) VPf.framewidth );  
	glBegin(GL_LINES);
	  glColor3fv( framecol );
	  glVertex3d( (coor + *(iwksp + i))->x1, 
		      (coor + *(iwksp + i))->y1,
		      (coor + *(iwksp + i))->z1 );
	  glVertex3d( (coor + *(iwksp + i))->x2,
		      (coor + *(iwksp + i))->y2,
		      (coor + *(iwksp + i))->z2 );
	glEnd(); 
	glDisable(GL_LINE_STIPPLE);
      }	
  }

  /* check if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers( togl );

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void
xcPointLineL2D(struct Togl *togl)
{
  int i;
  GLdouble loffset = VPf.PLradius / VPf.VPfactor;
    
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated(0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glLineWidth( (GLfloat) VPf.PLlinewidth );
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  for(i = 1; i <= tmp_nobjects; i++) { 
    if ( (coor + *(iwksp + i))->flag == BOND )
      {
	glBegin(GL_LINES);
	   /* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */
	   if ( VPf.unibond) glColor3fv(unibondCol);
	   else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] ); 
	   glVertex3d( (coor + *(iwksp + i))->x1, 
		       (coor + *(iwksp + i))->y1,
		       (coor + *(iwksp + i))->z1 );
	   glVertex3d( (coor + *(iwksp + i))->x2,
		       (coor + *(iwksp + i))->y2,
		       (coor + *(iwksp + i))->z2 );
	glEnd();
      }
    if ( (coor + *(iwksp + i))->flag == ATOM )
      {	
	/* first display point */
	/* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */
	glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	glPushMatrix();
	glTranslated( (coor + *(iwksp + i))->x1, 
		      (coor + *(iwksp + i))->y1,
		      (coor + *(iwksp + i))->z1 );
	glCallList( PointList );	
	glPopMatrix();
	 
	/* than atom label */

	makeAtomLabel2D ( *(iwksp + i), loffset + rf.w2, loffset + rf.h2 );
      }
  }

  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void
xcPointLineF2D(struct Togl *togl)
{
  int i;
  
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated (0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  for(i = 1; i <= tmp_nobjects; i++) { 
    if ( (coor + *(iwksp + i))->flag == BOND )
      {
	glLineWidth( (GLfloat) VPf.PLlinewidth );  
	glBegin(GL_LINES);
	   /* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */
	   if ( VPf.unibond) glColor3fv(unibondCol);
	   else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	   glVertex3d( (coor + *(iwksp + i))->x1, 
		       (coor + *(iwksp + i))->y1,
		       (coor + *(iwksp + i))->z1 );
	   glVertex3d( (coor + *(iwksp + i))->x2,
		       (coor + *(iwksp + i))->y2,
		       (coor + *(iwksp + i))->z2 );
	glEnd();
      }
    if ( (coor + *(iwksp + i))->flag == ATOM )
      {	
	/* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */
	glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	glPushMatrix();
	glTranslated( (coor + *(iwksp + i))->x1, 
		      (coor + *(iwksp + i))->y1,
		      (coor + *(iwksp + i))->z1 );
	glCallList( PointList );		
	glPopMatrix();	
      }
    if ( (coor + *(iwksp + i))->flag == FRAME )
      {
	if ( (coor + *(iwksp + i))->nat == 2 ) {
	  glLineStipple(2, FRAME_PAT);
	  glEnable(GL_LINE_STIPPLE);
	}
	glLineWidth( (GLfloat) VPf.framewidth );  
	glBegin(GL_LINES);
	  glColor3fv( framecol );
	  glVertex3d( (coor + *(iwksp + i))->x1, 
		      (coor + *(iwksp + i))->y1,
		      (coor + *(iwksp + i))->z1 );
	  glVertex3d( (coor + *(iwksp + i))->x2,
		      (coor + *(iwksp + i))->y2,
		      (coor + *(iwksp + i))->z2 );
	glEnd(); 
	glDisable(GL_LINE_STIPPLE);
      }
  }

  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void
xcPointLineFL2D(struct Togl *togl)
{
  int i;
  GLdouble loffset = VPf.PLradius / VPf.VPfactor;
    
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated (0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  for(i = 1; i <= tmp_nobjects; i++) { 
    if ( (coor + *(iwksp + i))->flag == BOND )
      {
	glLineWidth( (GLfloat) VPf.PLlinewidth );  
	glBegin(GL_LINES);
	   /* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */
	   if ( VPf.unibond) glColor3fv(unibondCol);
	   else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	   glVertex3d( (coor + *(iwksp + i))->x1, 
		       (coor + *(iwksp + i))->y1,
		       (coor + *(iwksp + i))->z1 );
	   glVertex3d( (coor + *(iwksp + i))->x2, 
		       (coor + *(iwksp + i))->y2,
		       (coor + *(iwksp + i))->z2 );
	glEnd();
      }
    if ( (coor + *(iwksp + i))->flag == ATOM )
      {	
	/* first display point */
	/* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */
	glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	glPushMatrix();
	glTranslated( (coor + *(iwksp + i))->x1, 
		      (coor + *(iwksp + i))->y1,
		      (coor + *(iwksp + i))->z1 );
	glCallList( PointList );	
	glPopMatrix();
	 
	/* than atom label */

	makeAtomLabel2D ( *(iwksp + i), loffset + rf.w2, loffset + rf.h2 );
      }
    if ( (coor + *(iwksp + i))->flag == FRAME )
      {
	if ( (coor + *(iwksp + i))->nat == 2 ) {
	  glLineStipple(2, FRAME_PAT);
	  glEnable(GL_LINE_STIPPLE);
	}
	glLineWidth( (GLfloat) VPf.framewidth );  
	glBegin(GL_LINES);
	  glColor3fv( framecol );
	  glVertex3d( (coor + *(iwksp + i))->x1, 
		      (coor + *(iwksp + i))->y1,
		      (coor + *(iwksp + i))->z1 );
	  glVertex3d( (coor + *(iwksp + i))->x2, 
		      (coor + *(iwksp + i))->y2,
		      (coor + *(iwksp + i))->z2 );
	glEnd(); 
	glDisable(GL_LINE_STIPPLE);
      }    
  }

  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void
xcBallStick1L2D(struct Togl *togl)
{
  int i;
  GLdouble x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated (0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glLineWidth( VPf.OUTlinewidth );
  glEnable( GL_CULL_FACE );
  glCullFace( GL_BACK );

  for(i = 1; i <= tmp_nobjects; i++) 
    {
      if ( (coor + *(iwksp + i))->flag == BOND )
	{
	  MakeSticks1(*(iwksp + i), 
		      &x1, &y1, &z1,  &x2, &y2, &z2,
		      &x3, &y3, &z3,  &x4, &y4, &z4);

	  glBegin(GL_POLYGON);
	    /* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */
	    if ( VPf.unibond) glColor3fv(unibondCol);
	    else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	    glVertex3d( x1, y1, z1);
	    glVertex3d( x2, y2, z2);
	    glVertex3d( x3, y3, z3);
	    glVertex3d( x4, y4, z4);
	  glEnd();

	  glBegin(GL_LINE_LOOP);
	    glColor3f( 0.0, 0.0, 0.0 );
	    glVertex3d( x1, y1, z1);
	    glVertex3d( x2, y2, z2);
	    glVertex3d( x3, y3, z3);
	    glVertex3d( x4, y4, z4);
	  glEnd();	      
	}
      else if ( (coor + *(iwksp + i))->flag == ATOM ) 
	{
	  /* first display ball */
	  glPushMatrix();
	  glTranslated( (coor + *(iwksp + i))->x1,
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	  glColor3fv (atm.col[(coor + *(iwksp + i))->sqn]);
	  xcBall2D (*(iwksp + i));
	  glPopMatrix();
	  
	  /* than atom label */

	  makeAtomLabel2D ( *(iwksp + i), CENTERED_LABEL, CENTERED_LABEL );
	}            
    }
  glDisable( GL_CULL_FACE );

  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void
xcBallStick1F2D(struct Togl *togl)
{
  int i;
  GLdouble x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated (0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glEnable( GL_CULL_FACE );
  glCullFace( GL_BACK );

  for(i = 1; i <= tmp_nobjects; i++) 
    {
      if ( (coor + *(iwksp + i))->flag == BOND )
	{
	  glLineWidth( VPf.OUTlinewidth );
	  MakeSticks1(*(iwksp + i), 
		      &x1, &y1, &z1,  &x2, &y2, &z2, 
		      &x3, &y3, &z3,  &x4, &y4, &z4);

	  glBegin(GL_POLYGON);
	    /* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */
	    if ( VPf.unibond) glColor3fv(unibondCol);
	    else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	    glVertex3d( x1, y1, z1);
	    glVertex3d( x2, y2, z2);
	    glVertex3d( x3, y3, z3);
	    glVertex3d( x4, y4, z4);
	  glEnd();

	  glBegin(GL_LINE_LOOP);
	    glColor3f( 0.0, 0.0, 0.0 );
	    glVertex3d( x1, y1, z1);
	    glVertex3d( x2, y2, z2);
	    glVertex3d( x3, y3, z3);
	    glVertex3d( x4, y4, z4);
	  glEnd();	      
	}
      if ( (coor + *(iwksp + i))->flag == ATOM ) 
	{
	  glLineWidth( VPf.OUTlinewidth );
	  glPushMatrix();
	  glTranslated( (coor + *(iwksp + i))->x1,
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	  glColor3fv (atm.col[(coor + *(iwksp + i))->sqn]);
	  xcBall2D (*(iwksp + i));
	  glPopMatrix();
  
	}            
      if ( (coor + *(iwksp + i))->flag == FRAME )
	{
	  if ( (coor + *(iwksp + i))->nat == 2 ) {
	    glLineStipple(2, FRAME_PAT);
	    glEnable(GL_LINE_STIPPLE);
	  }
	  glLineWidth( (GLfloat) VPf.framewidth );  
	  glBegin(GL_LINES);
	    glColor3fv( framecol );
	    glVertex3d( (coor + *(iwksp + i))->x1, 
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	    glVertex3d( (coor + *(iwksp + i))->x2,
			(coor + *(iwksp + i))->y2,
			(coor + *(iwksp + i))->z2 );
	  glEnd(); 
	  glDisable(GL_LINE_STIPPLE);
	}
    }
  glDisable (GL_CULL_FACE);
  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void
xcBallStick1FL2D(struct Togl *togl)
{
  int i;
  GLdouble x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated (0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glLineWidth( VPf.OUTlinewidth );
  glEnable( GL_CULL_FACE );
  glCullFace( GL_BACK );

  for(i = 1; i <= tmp_nobjects; i++) 
    {
      if ( (coor + *(iwksp + i))->flag == BOND )
	{
	  glLineWidth( VPf.OUTlinewidth );
	  MakeSticks1(*(iwksp + i), 
		      &x1, &y1, &z1,  &x2, &y2, &z2, 
		      &x3, &y3, &z3,  &x4, &y4, &z4);

	  glBegin(GL_POLYGON);
	    /* glColor3fv( atcol[(coor + *(iwksp + i))->nat] ); */
	    if ( VPf.unibond) glColor3fv(unibondCol);
	    else glColor3fv( atm.col[(coor + *(iwksp + i))->sqn] );
	    glVertex3d( x1, y1, z1);
	    glVertex3d( x2, y2, z2);
	    glVertex3d( x3, y3, z3);
	    glVertex3d( x4, y4, z4);
	  glEnd();

	  glBegin(GL_LINE_LOOP);
	    glColor3f( 0.0, 0.0, 0.0 );
	    glVertex3d( x1, y1, z1);
	    glVertex3d( x2, y2, z2);
	    glVertex3d( x3, y3, z3);
	    glVertex3d( x4, y4, z4);
	  glEnd();	      
	}
      if ( (coor + *(iwksp + i))->flag == ATOM ) 
	{
	  /* first display ball */
	  glLineWidth( VPf.OUTlinewidth );  
	  glPushMatrix();
	  glTranslated( (coor + *(iwksp + i))->x1,
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	  glColor3fv (atm.col[(coor + *(iwksp + i))->sqn]);
	  xcBall2D (*(iwksp + i));
	  glPopMatrix();
	  	  
	  /* than atom label */

	  makeAtomLabel2D ( *(iwksp + i), CENTERED_LABEL, CENTERED_LABEL );
	}            
      if ( (coor + *(iwksp + i))->flag == FRAME )
	{
	  if ( (coor + *(iwksp + i))->nat == 2 ) {
	    glLineStipple(2, FRAME_PAT);
	    glEnable(GL_LINE_STIPPLE);
	  }
	  glLineWidth( (GLfloat) VPf.framewidth );  
	  glBegin(GL_LINES);
	    glColor3fv( framecol );
	    glVertex3d( (coor + *(iwksp + i))->x1, 
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	    glVertex3d( (coor + *(iwksp + i))->x2,
			(coor + *(iwksp + i))->y2,
			(coor + *(iwksp + i))->z2 );
	  glEnd();
	  glDisable(GL_LINE_STIPPLE);
	}
    }
  glDisable (GL_CULL_FACE);

  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void
xcBallStick2L2D(struct Togl *togl)
{
  int i;
  
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated (0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glLineWidth( VPf.OUTlinewidth );
  glEnable( GL_CULL_FACE );
  glCullFace( GL_BACK );

  for(i = 1; i <= tmp_nobjects; i++) 
    {
      if ( (coor + *(iwksp + i))->flag == BOND )	
	{	
	  /* MakeSticks2(*(iwksp + i), (coor + *(iwksp + i))->nat); */
	  MakeSticks2(*(iwksp + i), (coor + *(iwksp + i))->sqn, BALL);
	}
      if ( (coor + *(iwksp + i))->flag == ATOM ) 
	{
	  /* first display ball */
	  glPushMatrix();
	  glTranslated( (coor + *(iwksp + i))->x1,
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	  glColor3fv (atm.col[(coor + *(iwksp + i))->sqn]);
	  xcBall2D (*(iwksp + i));
	  glPopMatrix();

	  /* than atom label */

	  makeAtomLabel2D ( *(iwksp + i), CENTERED_LABEL, CENTERED_LABEL );
	}      
    }
  glDisable( GL_CULL_FACE );

  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void
xcBallStick2F2D(struct Togl *togl)
{
  int i;
  
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated (0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glLineWidth( VPf.OUTlinewidth );
  glEnable( GL_CULL_FACE );
  glCullFace( GL_BACK );

  for(i = 1; i <= tmp_nobjects; i++) 
    {
      if ( (coor + *(iwksp + i))->flag == BOND )	
	{	  
	  /* MakeSticks2(*(iwksp + i), (coor + *(iwksp + i))->nat); */
	  glLineWidth( VPf.OUTlinewidth );
	  MakeSticks2(*(iwksp + i), (coor + *(iwksp + i))->sqn, BALL);
	}
      if ( (coor + *(iwksp + i))->flag == ATOM ) 
	{
	  glLineWidth( VPf.OUTlinewidth );
	  glPushMatrix();
	  glTranslated( (coor + *(iwksp + i))->x1,
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	  glColor3fv (atm.col[(coor + *(iwksp + i))->sqn]);
	  xcBall2D (*(iwksp + i));
	  glPopMatrix();	  
	}      
      if ( (coor + *(iwksp + i))->flag == FRAME )
	{
	  if ( (coor + *(iwksp + i))->nat == 2 ) {
	    glLineStipple(2, FRAME_PAT);
	    glEnable(GL_LINE_STIPPLE);
	  }
	  glLineWidth( (GLfloat) VPf.framewidth );  
	  glBegin(GL_LINES);
	    glColor3fv( framecol );
	    glVertex3d( (coor + *(iwksp + i))->x1, 
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	    glVertex3d( (coor + *(iwksp + i))->x2,
			(coor + *(iwksp + i))->y2,
			(coor + *(iwksp + i))->z2 );
	  glEnd(); 
	  glDisable(GL_LINE_STIPPLE);
	}
    }
  glDisable (GL_CULL_FACE);

  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


void
xcBallStick2FL2D(struct Togl *togl)
{
  int i;
  
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  if (VPf.perspective) 
    glTranslated (0.0, 0.0, persp.shiftZ);
  else
    glTranslated (0.0, 0.0, -ort.size);

  /* fog & antialias */
  xcFog       (togl, VPf.fog, VPf.perspective);
  xcAntiAlias (VPf.antialias);

  glLineWidth( VPf.OUTlinewidth );
  glEnable( GL_CULL_FACE );
  glCullFace( GL_BACK );

  for(i = 1; i <= tmp_nobjects; i++) 
    {
      if ( (coor + *(iwksp + i))->flag == BOND )	
	{	  
	  /* MakeSticks2(*(iwksp + i), (coor + *(iwksp + i))->nat); */
	  glLineWidth( VPf.OUTlinewidth );
	  MakeSticks2(*(iwksp + i), (coor + *(iwksp + i))->sqn, BALL);
	}
      if ( (coor + *(iwksp + i))->flag == ATOM ) 
	{
	  /* first display ball */
	  glLineWidth( VPf.OUTlinewidth );
	  glPushMatrix();
	  glTranslated( (coor + *(iwksp + i))->x1,
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	  glColor3fv (atm.col[(coor + *(iwksp + i))->sqn]);
	  xcBall2D (*(iwksp + i));
	  glPopMatrix();

	  /* than atom label */

	  makeAtomLabel2D ( *(iwksp + i), CENTERED_LABEL, CENTERED_LABEL );
	}      
      if ( (coor + *(iwksp + i))->flag == FRAME )
	{
	  if ( (coor + *(iwksp + i))->nat == 2 ) {
	    glLineStipple(2, FRAME_PAT);
	    glEnable(GL_LINE_STIPPLE);
	  }
	  glLineWidth( (GLfloat) VPf.framewidth );  
	  glBegin(GL_LINES);
	    glColor3fv( framecol );
	    glVertex3d( (coor + *(iwksp + i))->x1, 
			(coor + *(iwksp + i))->y1,
			(coor + *(iwksp + i))->z1 );
	    glVertex3d( (coor + *(iwksp + i))->x2,
			(coor + *(iwksp + i))->y2,
			(coor + *(iwksp + i))->z2 );
	  glEnd(); 
	  glDisable(GL_LINE_STIPPLE);
	}
    }
  glDisable (GL_CULL_FACE);

  /* chack if coordinate sistem should be displayed */
  if ( VPf.xyzOn ) xcDisplayXYZ();

  glFlush();
  Togl_SwapBuffers(togl);

  /* every snapshot movie making */
  if ( makeMovie.doit && makeMovie.mode == MOVIE_MODE_EVERY_SNAPSHOT && !makeMovie.printing ) {
    fprintf(stderr,"making snapshot\n");
    createMoviePPMFrame(togl);
  }  
}


/* this procedure turn labels On */
void
xcTurnLabelsOn(void)
{
  /* if dimType is XC_3D all that's needed to do is to set VPf.labelsOn to 1,
   * but if dimType is XC_2D, than it's a little harder
   */

  VPf.labelsOn = 1;  

  /* XC_2D */
  if ( dimType == XC_2D )  WhichDispFunc();
}


/* this procedure turn labels Off */
void
xcTurnLabelsOff(void)
{
  /* if dimType is XC_3D all that's needed to do is to set VPf.labelsOn to 0,
   * but if dimType is XC_2D, than it's a little harder
   */

  VPf.labelsOn = 0;

  /* XC_2D */
  if ( dimType == XC_2D )  WhichDispFunc();
}


/* this procedure turn frames On */
void
xcTurnFramesOn(void)
{
  /* if dimType is XC_3D all that's needed to do is to set VPf.framesOn to 1,
   * but if dimType is XC_2D, than it's a little harder
   */

  VPf.framesOn = 1;

  /* XC_2D */
  if ( dimType == XC_2D )  WhichDispFunc();
}


/* this procedure turn labels Off */
void
xcTurnFramesOff(void)
{
  /* if dimType is XC_3D all that's needed to do is to set VPf.framesOn to 0,
   * but if dimType is XC_2D, than it's a little harder
   */

  VPf.framesOn = 0;

  /* XC_2D */
  if ( dimType == XC_2D )  WhichDispFunc();
}
 

/* ========================================================================= */
/* this proc is used by xcTurnLabelsOn & xcTurnLabelsOff to determine which  */
/* displayFunc to use======================================================= */
void
WhichDispFunc(void)
{
  printf("In WhichDispFunc, dimType=%d",dimType);
  fflush(stdout);
  
  /* are we in WF mode */
  if ( xcDisplay == xcWireFrame2D ||
       xcDisplay == xcWireFrameL2D ||
       xcDisplay == xcWireFrameF2D ||
       xcDisplay == xcWireFrameFL2D ) {
    xcAssignDisplayFunc("WF");
  }
  
  else if ( xcDisplay == xcPointLine2D ||
	    xcDisplay == xcPointLineL2D ||
	    xcDisplay == xcPointLineF2D ||
	    xcDisplay == xcPointLineFL2D )
    xcAssignDisplayFunc("PL");

  else if ( xcDisplay == xcPipeBall2D )
    xcAssignDisplayFunc("PB");
  
  else if ( xcDisplay == xcBallStick12D ||
	    xcDisplay == xcBallStick1L2D ||
	    xcDisplay == xcBallStick1F2D ||
	    xcDisplay == xcBallStick1FL2D )
    xcAssignDisplayFunc("BS1");
  
  else if ( xcDisplay == xcBallStick22D ||
	    xcDisplay == xcBallStick2L2D ||
	    xcDisplay == xcBallStick2F2D ||
	    xcDisplay == xcBallStick2FL2D )
    xcAssignDisplayFunc("BS2");

  else if ( xcDisplay == xcSpaceFill2D )
    xcAssignDisplayFunc("SF");
}


/* ========================================================================= */
/* this is used by "XC_UpdateStrCmd", to update DisplayFunc                  */
/* ========================================================================= */
void
UpdateDispFunc(void)
{
  printf("In UpdateDispFunc, dimType=%d",dimType);
  fflush(stdout);
  if ( dimType == XC_2D ) {
    WhichDispFunc();
    /* xcMakeBallLists(); */
  }
  else if ( dimType == XC_3D ) {
    /*
      if (is.spacefillmode) xcMakeSpaceFill3DLists();
      if (is.ballmode)      xcMakeBall3DLists();
      if (is.stickmode)     xcMakeStick3DLists();
    */
    /*    if (VPf.framesOn)     xcMakeFrame3DLists();*/
    /*      if (VPf.labelsOn) { */
    /*        xcMakeBallLabel3DList(); */
    /*        xcMakeSpaceLabel3DList(); */
/*      } */
  }
}    
    

/* ========================================================================= */
/* this is used by "XC_UpdateStrCmd", to update Projection  (glOrtho)        */
/* ========================================================================= */
void
UpdateProjection(void)
{
  if ( dimType == XC_2D ) {
    if ( xcDisplay == xcWireFrame2D ||
	 xcDisplay == xcWireFrameL2D ||
	 xcDisplay == xcWireFrameF2D ||
	 xcDisplay == xcWireFrameFL2D ||
	 xcDisplay == xcPointLine2D ||
	 xcDisplay == xcPointLineL2D ||
	 xcDisplay == xcPointLineF2D ||
	 xcDisplay == xcPointLineFL2D )
      xcMakeProjection2D("WF");  /* WF and PL have the same "projections" */
  
    else if ( xcDisplay == xcPipeBall2D ||
	      xcDisplay == xcBallStick12D ||
	      xcDisplay == xcBallStick1L2D ||
	      xcDisplay == xcBallStick1F2D ||
	      xcDisplay == xcBallStick1FL2D ||
	      xcDisplay == xcBallStick22D ||
	      xcDisplay == xcBallStick2L2D ||
	      xcDisplay == xcBallStick2F2D ||
	      xcDisplay == xcBallStick2FL2D )
      xcMakeProjection2D("BS1");  /* BS1 and BS2 have the same "projections" */
    else if ( xcDisplay == xcSpaceFill2D )
      xcMakeProjection2D("SF");
  }
  else if ( dimType == XC_3D ) {
    if (is.spacefillmode) xcMakeProjection3D("space");
    if (is.ballmode || is.stickmode) xcMakeProjection3D("balls");
  }
}   
    


