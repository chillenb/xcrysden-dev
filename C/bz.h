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
 * Source: $XCRYSDEN_TOPDIR/C/bz.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

#define BZ_MAXPOINTS    100
#define BZ_MAXALLPOINTS 500

#define BZ_PRIMCELL 0
#define BZ_CONVCELL 1

#define BZ_OFFSET   0.5 /* offset for BZ in canvas */

#define BZ_CENTERPOINT 1
#define BZ_EDGEPOINT   2
#define BZ_LINEPOINT   4
#define BZ_POLYPOINT   8
#define BZ_LINECENTER  16

#define BZ_THIN_POLYGON_WIDTH 1
#define BZ_BOLD_POLYGON_WIDTH 3
#define BZ_THIN_POINT_SIZE    3.0
#define BZ_BOLD_POINT_SIZE    4.0
#define BZ_CENTERPOINT_SIZE   5.0

#define BZ_ISS 100 /* for xc_bz iss <canvasName> */

typedef struct {
  float recvec[3][3]; /* reciprocal vectors */  
  float rot_recvec[3][3];
  float vp_recvec[9][2];
  float reccoor[BZ_MAXALLPOINTS][3];
  int   recvec_is_rendered; /* are recvectors rendered */
  int   npoly; /* number of polygons */
  int   nvert[WIGNERSEITZ_MAXPOLY]; /*number of vertices for each polygon*/
  float max;   /* maximum value of coordinates (either X,Y,Z) */
  float poly[WIGNERSEITZ_MAXPOLY][WIGNERSEITZ_MAXVERTEX][3]; /*coord of poly */
  float norm[WIGNERSEITZ_MAXPOLY][WIGNERSEITZ_MAXVERTEX][3];
  int   polyindex[BZ_MAXALLPOINTS]; /* index of poly for Z-orientation */
  int   n_edgepoint; /* number of edgepoint */
  int   n_linepoint;
  int   n_polypoint;
  int   npoint;  /* number of special points */
  int   nZpoint; /* number of special points + selected-lines; 
                    used for Z orientation */
  float point[BZ_MAXALLPOINTS][3]; /* original coords of points */
  int   linepoint_index[BZ_MAXALLPOINTS][2]; /* which polygons share 
                                                the linepoint */
  float rotmat[3][3]; /* rotation matrix */
  float rot_poly[WIGNERSEITZ_MAXPOLY][WIGNERSEITZ_MAXVERTEX][3]; /* rotation coordinates of polygons */
  float rot_norm[WIGNERSEITZ_MAXPOLY][WIGNERSEITZ_MAXVERTEX][3];
  float rot_points[BZ_MAXALLPOINTS][3]; /* rotation coords of points */
  float vp_poly[WIGNERSEITZ_MAXPOLY][WIGNERSEITZ_MAXVERTEX][2]; /* viewport coordinates of polygons */
  float vp_points[BZ_MAXALLPOINTS][2]; /* viewport coordinates of points */
  int   pointtype[BZ_MAXALLPOINTS]; /* type of point */
  int   isrendered; /* is BZ already rendered */
  int   nselected; /* number of selected points */                   
  int   selected_is_rendered[BZ_MAXALLPOINTS]; /* is selected point rendered */
  int   selectedID[BZ_MAXALLPOINTS]; /* which is selected points -> ID od ordinary point */
  int   n_sline; /* number of selected lines (sline) */
  float sline[BZ_MAXALLPOINTS][6]; /* coorinates of sline */
  int   slinetype[BZ_MAXALLPOINTS]; /* type of sline */
  int   sline_index[BZ_MAXALLPOINTS][2]; /* which polygons share the sline; used just for BZ_EDGEPOINT<->BZ_LINEPOINT lines */
  int   sline_is_rendered[BZ_MAXALLPOINTS]; /* is sline rendered */
  float rot_sline[BZ_MAXALLPOINTS][6]; /* rotation coordinates of sline */
  float vp_sline[BZ_MAXALLPOINTS][4]; /* viewport coordinates of sline */
  float slcenter[BZ_MAXALLPOINTS][3]; /* coordinates of sline center */
  float rot_slcenter[BZ_MAXALLPOINTS][3]; /* rot coords of sline center */
} BRILLOUINZONE;
BRILLOUINZONE bz[2];

