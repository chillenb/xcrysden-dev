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
 * Source: $XCRYSDEN_TOPDIR/C/xcfunc.h
 * ------                                                                    *
 * Copyright (c) 1996-2003 by Anton Kokalj                                   *
 *****************************************************************************

*/

/*****************************************************************************/
/* XC_CPP_* are pre-processor flags */
#ifndef XC_CPP_STRUCT
#   include "struct.h"
#endif

#ifndef XC_CPP_GLPARAM
#   include "xcGLparam.h"
#endif

/* ===================================================================== *
 *                   FUNCTION PROTOTYPES FOR XCRYSDEN                    *
 * ===================================================================== */

/* 3D.c */
extern void LoadNull(GLdouble matrix[][4]);
extern void LoadNull44One(GLdouble matrix[][4]);
extern void LoadIdentity(GLdouble matrix[][4]);
extern void GetMajorMat(void);
extern void MajorMatToVec(void);
extern void GetOldMat(void);
extern void GetRotXYMat(GLdouble fiX, GLdouble fiY);
extern void GetRotX(GLdouble fiX);
extern void GetRotY(GLdouble fiY);
extern void GetRotZ(GLdouble fiZ);
extern void MakeCylinderCoor(void);


/* --- auxilproc.c: --- */
/* auxilproc.c */
extern double dist6(double x1, double x2, double y1, double y2, double z1, double z2);
extern double dist3(double x, double y, double z);
extern double distdv(double *vec);
extern int normalizepv(double *x, double *y, double *z);
extern int iround(double x);
extern void Rotate(double *x, double *y, double *z, double cosfi, double sinfi);
extern void xcMat44Copyd(double mat1[][4], double mat2[][4], int n, int m);
extern int normalizepvf(float *x, float *y, float *z);
extern int normalizepvfv(float *vec);
extern void RevertVectorfv(float *vec);
extern float dist3f(float x, float y, float z);
extern float distfv(float *v);
extern int iroundf(float x);
extern void xcDeleteBitFlags(unsigned int *var, unsigned int flags);
extern void MatVecMult33f(float mat[3][3], float vec[3], float new[3]);
extern void MatMult33f(float a[3][3], float b[3][3], float new[3][3]);
extern void MatCopy33f(float a[3][3], float b[3][3]);
extern void MatToZero33f(float a[3][3]);
extern int IsSamePointvf(float *a, float *b, float limit);
extern int IsEqualf(float a, float b, float limit);
extern int WriteArgv(char *argv[], int *argc, char *fmt, ...);
extern void VecSum2f(float vec1[], float vec2[], int i, int j, float res[]);
extern void VecSum3f(float vec1[], float vec2[], float vec3[], int i, int j, int k, float res[]);
extern void VecProduct3f(float v0[], float v1[], float res[]);
extern float CramerRule2x2Detf(float d[2][3], int i);
extern void xcRotate3fv(float fi, float *v, float *m, float *p);
extern double det3x3d(const double *v0, const double *v1, const double *v2);
extern float det3x3f(const float *v0, const float *v1, const float *v2);

/* --- colorplane.c: --- */
extern void ColorPlane(int obj, int cb, int fn, int nx, int ny);
extern void GetIsoLine2D_Attributes(int obj, int cb, int fn);
extern float xcLinf(float value);
extern float xcLogf(float value);
extern float xcLog10f(float value);
extern float xcSqrtf(float value);
extern float xcRoot3f(float value);
extern float xcGaussf(float value);
extern float xcSlaterf(float value);
extern float xcExpf(float value);
extern float xcExp2f(float value);

/* --- detnsplit.c: --- */
extern int DetNBonds(int *nsplit);

/* hpsort.c */
extern int hpsort_index(size_t count, const double *data, int *p);
extern int hpsort_index1(size_t count, const double *data, int *p);
extern int hpsort_index_f(size_t count, const float *data, int *p);
extern int hpsort_index1_f(size_t count, const float *data, int *p);

/* --- isorender.c: --- */
extern void xcRenderColorplane(int obj);
extern void xcRenderIsosurf(int obj);
extern void xcRenderIsoLine2D(int obj);

/* --- lighting.c: --- */
extern void CalcLightPosition (int il, GLdouble size);
extern void LoadBlendfunc_And_Frontface(void);
extern void LoadIsoMaterial(int type);
extern void LoadLights(void);
extern void LoadLightModel(int type);
extern void LoadStructMaterial(void);
extern void LoadVoronoiMaterial(int voronoi_type);
extern void LoadCageOrVecMaterial(int type);


/* --- loadValues.c: --- */
extern double GetDefault(int var);
extern double GetValue(int var);
extern void LoadDefaultValues(void);
#ifdef TOGL_H
   extern void LoadNewValue(struct Togl *togl, int var, double value1, 
			    double value2, double value3, double value4);
   extern void ResetVar(struct Togl *togl, int var);
#endif

/* --- polygonise.c: --- */

extern void WrapperSurfSmoothing(int ilevel);
extern void MarchingCubes(float isolevel, int ilevel, int smooth_nstep, float smooth_weight, int algorithm, int shade_model, int normals_model);

#ifdef ISOSURF_H
extern int Polygonise(GRIDCELL grid, float isolevel, GRID ijk, TRIANGLE *triangles, TRIG_INFO *trig_info);
extern XYZ VertexInterp(float isolevel, XYZ p1, XYZ p2, float valp1, float valp2);
extern int PolygoniseTetrahedral(GRIDCELL grid, float iso, GRID ijk, TRIANGLE *triangl, TRIG_INFO *trig_info);
#endif

/* polygonise_auxil.c */
#ifdef ISOSURF_H
extern void VertexRecognise(GRID grd, int ntriangl, TRIANGLE *triangles, int *triangl_status, TRIG_INFO *tri2verIN, XYZ *vertex, int *vertex_status, VERT2TRIG_INFO (*ver2triIN)[20 ], int *nver2triIN, GRID_INFO *grid_info, int *nvertex);
extern void SurfSmoothing(int nstep, float weight, XYZ *vertex, int *vertex_status, VERT2TRIG_INFO (*ver2triIN)[20 ], TRIG_INFO *tri2verIN, int *nver2triIN, int nvertex);
extern void SurfNmlAver(int ntriangl, int nvertex, XYZ *nml, XYZ *vertex, TRIG_INFO *tri2verIN, float *sign, int is_tetrahedral);
#endif
#ifdef XC_CPP_XYZ
extern XYZ TriangleWeightNormalS(XYZ p0, XYZ p1, XYZ p2, float sign);
extern XYZ TriangleWeightNormal(XYZ p0, XYZ p1, XYZ p2);
extern XYZ VertexNormal(XYZ v1, XYZ v2);
extern float distXYZ(XYZ v);
extern int normalizeXYZ(XYZ *v);
#endif

/* gridNormals.c */
#ifdef ISOSURF_H
extern void gridNormals(GRIDVERTEX ***g, GRID grd);
extern void gradient_SurfNml(GRIDVERTEX ***g, GRID grd, XYZ spanning_vec[3], int ntriangl, TRIG_INFO *tri2verIN, int nvertex, XYZ *vertex, XYZ *normal);
extern void gradient_SurfNmlNew(float origin[], float incr_vec[3][3], float ***g, GRID  gr, int   ntriangl, TRIG_INFO *tri2verIN, int  nvertex, XYZ  vertex[], XYZ  normal[]);
#endif

/* readisodata.c */
extern void WriteBinVertexFP(void);
#ifdef XC_CPP_XYZ
extern int ReadIsoData(int framenum[ISODATA_MAXFILES][ISODATA_MAXSTACK][2], int ntimes);
extern int ReadBlock0(int type, int s0, int s1, int s2, int s3);
#endif

/* --- readstrf.c: --- */
/* readstrf.c */
extern int ReadStructFile(FILE *fp, const char *file, int format, int mode);
extern void MakeBonds(void);
extern void MallocCoor(void);
extern int ReadXSF(FILE *fp);
extern int ReadVec(FILE *fp, double vec[][4]);
extern void FindMaxRad(void);
extern void ParseAtomName(char *atomName);


/* --- remakestr.c: --- */
#ifdef TOGL_H
   extern void ReDisplay(struct Togl *togl, GLuint remake3D);
#endif
extern void ReMakeStr(void);

/* --- vectors.c: --- */
extern void  VecRotateX(double cosfi, double sinfi);
extern void  VecRotateXY(double cosfiX, double sinfiX, double cosfiY, double sinfiY);
extern void  VecRotateXY(double cosfiX, double sinfiX, double cosfiY, double sinfiY);
extern void  VecRotateY(double cosfi, double sinfi);
extern void  VecRotateZ(double cosfi, double sinfi);
extern void  vecMatToVec(double major[4][4], double vec[16]);
extern void vecVecToMat(double vec[16], double major[4][4]);
extern void  VecRotTo_XY(void);
extern void  VecRotTo_XZ(void);
extern void  VecRotTo_YZ(void);
extern void  xcRotatefv(float fi, float u[]);
extern void  VecRotTo_AB(void);
extern void  VecRotTo_AC(void);
extern void  VecRotTo_BC(void);

/* --- xcAppInit.c: --- */
#ifdef _TCL
extern int Tcl_AppInit(Tcl_Interp *interp);
extern int  XC_DisplayMode2DCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int  XC_DisplayMode3DCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int  XC_DrawStyle3DCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int  XC_GetDefaultCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int  XC_GetValueCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int  XC_MesaContextCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int  XC_NewValueCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int  XC_OpenStrCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int  XC_PointSizeCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int  XC_ResetVarCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int  XC_ShadeModel3DCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
#endif /* _TCL */
#ifdef TOGL_H
extern int XC_RotateCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_TranslateCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_Dump2EpsCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
#endif
#ifdef _TCL
extern int  XC_UpdateStrCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_DisplayCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
#endif 
extern int xcRotateXY(double fiX, double fiY);



/* --- xcAtomAdd.c: --- */
#ifdef _TCL
extern int  XC_AtomAddCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, const char *argv[]);
#endif


/* --- xcDisplayFunc.c: --- */
extern GLuint xcGenLists(GLsizei i);
extern GLuint findList1(int lindex, GLdouble *paramArray, int size);
extern GLuint makeModelPtr1(int lindex, GLdouble sizeArray[], int count);
#ifdef TOGL_H
extern void xcDisplayFunc(void (*Func)(struct Togl *togl));
extern int xcToglDisplayFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv);
extern void xcDummyDisplay(struct Togl *togl);
#endif
extern void xcGenDispList(void);
extern void xcMakePointList(void);
#ifdef TOGL_H
extern void xcWireFrame2D(struct Togl *togl);
extern void xcPointLine2D(struct Togl *togl);
#endif
extern void xcBall2D(int iatom);
extern void xcSmallBall2D(int iatom);
extern void xcBigBall2D(int iatom);
#ifdef TOGL_H
extern void xcBallStick12D(struct Togl *togl);
extern void xcPipeBall2D(struct Togl *togl);
extern void xcBallStick22D(struct Togl *togl);
extern void xcSpaceFill2D(struct Togl *togl);
#endif
extern void EnableOr2D_Or3D(void);
extern void xcWireSphere(GLdouble radius);
extern void xcSolidSphere(GLdouble radius);
extern void xcWireCylinder(GLdouble radius, GLdouble height);
extern void xcSolidCylinder(GLdouble radius, GLdouble height);
extern void xcSolidBond(GLdouble radius, GLdouble height, int bondFlag);
extern void xcSolidCone(GLdouble baseradius, GLdouble topradius, GLdouble height);
extern void xcRenderSolidBalls3D(void);
extern void xcRenderWireBalls3D(void);
extern void xcRenderBonds(GLenum type);
extern void xcRenderSolidSpaceFills3D(void);
extern void xcRenderWireSpaceFills3D(void);
extern void xcRenderSolidFrame3D(void);
extern void xcRenderWireFrame3D(void);
extern void xcRenderLineFrame3D(void);
extern void xcMakeBallLabel3D(void);
extern void xcMakeSpaceLabel3D(void);
extern void xcMaybeDestroyLists(void);
extern void xcMakeProjection2D(const char *mode);
extern void xcMakeProjection3D(const char *mode);
extern void xcDisplayXYZ(void);
#ifdef TOGL_H
extern void xcDisplay3D(struct Togl *togl);
extern void xcClearScreen(struct Togl *togl);
extern void xcChangeBackground(struct Togl *togl, GLclampf bckg[4]);
extern void xcTurnXYZOff(struct Togl *togl);
#endif

/* --- xcDisplayFunc2.c: --- */
extern void UpdateDispFunc(void);
extern void UpdateProjection(void);
extern void RewriteCoor(GLenum type);
extern void xcTurnLabelsOff(void);
extern void xcTurnLabelsOn(void);
extern void xcTurnFramesOff(void);
extern void xcTurnFramesOn(void);
extern void makeAtomLabel2D(int ith, GLdouble Xoffset, GLdouble Yoffset);
#ifdef TOGL_H
   extern void xcBallStick1F2D(struct Togl *togl);
   extern void xcBallStick1FL2D(struct Togl *togl);
   extern void xcBallStick1L2D(struct Togl *togl);
   extern void xcBallStick2F2D(struct Togl *togl);
   extern void xcBallStick2FL2D(struct Togl *togl);
   extern void xcBallStick2L2D(struct Togl *togl);
   extern void xcPointLineF2D(struct Togl *togl);
   extern void xcPointLineFL2D(struct Togl *togl);
   extern void xcPointLineL2D(struct Togl *togl);
   extern void xcWireFrameF2D(struct Togl *togl);
   extern void xcWireFrameFL2D(struct Togl *togl);
   extern void xcWireFrameL2D(struct Togl *togl);
#endif

#ifdef _TCL
/* --- xcExit.c: --- */
extern int XC_ExitCmd(ClientData clientData,Tcl_Interp *interp,
		      int argc, const char *argv[]);


/* --- xcGLparam.c: --- */
extern int  XC_GetGLparamCmd(ClientData clientData, Tcl_Interp *interp,
			     int argc, const char *argv[]);
extern int  XC_SetGLparamCmd(ClientData clientData, Tcl_Interp *interp,
			     int argc, const char *argv[]);
extern int xcSplitList( int code, Tcl_Interp *interp, const char *argv[], GetGlParam *var);
extern void xcTclListError( int code, Tcl_Interp *interp, const char **argv );


/* --- xcIsoSurf.c: --- */
extern int XC_IsoCmd(ClientData clientData, Tcl_Interp *interp,
		     int argc, const char *argv[]);
extern int XC_IsodataCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, const char *argv[]);
extern int XC_IsofilesCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, const char *argv[]);
extern int XC_IsoplaneCmd(ClientData clientData, Tcl_Interp *interp,
			   int argc, const char *argv[]);
extern int XC_IsopointsCmd(ClientData clientData, Tcl_Interp *interp,
			   int argc, const char *argv[]);
extern int XC_IsosignCmd(ClientData clientData, Tcl_Interp *interp,
			 int argc, const char *argv[]);
extern int XC_IsostackCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, const char *argv[]);
extern int  XC_IsosurfCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, const char *argv[]);
extern int  XC_IsoplaneCmd(ClientData clientData, Tcl_Interp *interp,
			   int argc, const char *argv[]);
#endif /* _TCL */
extern void xcIsoError(void);
extern int IsIsoPlane123( int type );


/* --- xcLabels.c: --- */
extern void makeAtomLabels(void);
extern void makeCrdList(void);
extern void makeRasterFont(void);
extern void makeTemp3D2DList(void);
extern void makeXYZLabels(void);


/* --- xcMesaContext.c: --- */


/* --- xcMesaEvent.c: --- */
/* extern void XC_MesaEvent(ClientData clientData, XEvent *eventPtr);
   extern void xcSwapBuffers(void);*/


/* --- xcMesaWin.c: --- */
/* extern int XC_MesaWinCmd(ClientData clientData, Tcl_Interp *interp, 
   int argc, const char *argv[]);*/


/* --- xcSelect.c: --- */
/* xcSelect.c */
#ifdef _TCL
extern int XC_SelectCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_DeselectCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
#endif
extern int xcSelectSqn(int x, int y);
extern void GetSelCylinderPar(double x21, double y21, double z21, double *xrvb, double *yrvb, double *zrvb, double *fibond, double *bondl);
#ifdef TOGL_H
extern void xcPointLineSel(struct Togl *togl);
extern void xcBallStick1Sel(struct Togl *togl);
extern void xcBallStick2Sel(struct Togl *togl);
#endif
extern void xcRenderSelAtoms3D(void);
extern void xcRenderSelBonds3D(void);

/* --- xcballstick.c: --- */
extern void  MakeArcPoints(void);
extern int MakeSticks1(int i,
		       GLdouble *x1, GLdouble *y1, GLdouble *z1, 
		       GLdouble *x2, GLdouble *y2, GLdouble *z2, 
		       GLdouble *x3, GLdouble *y3, GLdouble *z3, 
		       GLdouble *x4, GLdouble *y4, GLdouble *z4);
extern int MakeSticks2(int i, int col, int flag); 


/* xcdebug.c */
extern void xcdebug(const char *text);
extern void xcErrDebug(const char *text);
extern void breakpoint(const char *text);


/* --- xcviewport.c: --- */
extern void xcViewPort(void);
extern void MaybeClipAndMakeProjection(void);
extern void Screen2Model_Coords(int xs, int ys, float *xm, float *ym);

/* --- voronoi.c: --- */
extern void xcRenderVoronoi(void);


#ifdef _TCL
/* --- xcBz.c --- */
extern int XC_BzCmd(ClientData clientData,Tcl_Interp *interp,
		    int argc, const char *argv[]);

/* --- xcSuperCell --- */
extern int  XC_SuperCellCmd(ClientData clientData,Tcl_Interp *interp, int argc, const char *argv[]);
#endif /* _TCL */
extern void (*xcSuperCell)(void);
extern void xcSuperCellFunc( void (*Func)(void) );
#ifdef XC_CPP_VECTOR
   extern void SetUnitCellCage( double vec[4][4], CellCage *cage );
#endif

/* --- xcPrimitives.c --- */
#ifdef XC_CPP_VECTOR
   extern void xcSolidCage( CellCage cage );
   extern void xcWireCage( CellCage cage );
   extern void xcSolidVector( RenderVectors vec );
#endif
extern void xcParallelogram( int type, float p[4][3], float nml[3] );


/* --- setOpenGLState.c --- */
#ifdef XC_CPP_GLPARAM
   extern void SetCageOGLState( int dim, BLENDFUNC blend, GLint *shade_model, 
				GLboolean *two_side, GLboolean *cull_face );
   extern void DisableCageOGLState( int dim, GLint shade_model, 
				    GLboolean two_side, GLboolean cull_face );
#endif


/* memory.c */
extern void *xcMalloc(size_t size);
extern void *xcCalloc(size_t nmemb, size_t size);
extern void *xcRealloc(void *ptr, size_t size);
extern void xcFree(FREE_ARG ptr);
extern void xcError(char error_text[]);
extern void xcFree_Vectorf(float *v);
extern void xcFree_Matrixf(float **m);
extern void xcFree_Tensor3f(float ***t);
extern void xcFree_Tensor3i(int ***t);
extern void xcFree_ReallocatedTensor3f_00to11(float ***rt);
extern float *xcMallocVectorf(long n);
extern float **xcMallocMatrixf(long nr, long nc);
extern float ***xcMallocTensor3f(long nr, long nc, long nd);
extern int ***xcMallocTensor3i(long nr, long nc, long nd);
extern short int ***xcMallocTensor3si(long nr, long nc, long nd);
extern float ***xcReallocTensor3f_00to11(float ***t, long nr, long nc, long nd);
extern void *cryMem_Malloc_W(size_t size, char *where);

/* --- splineInt.c --- */
extern void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a);
extern void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
extern void xcRegularSplineInt2(float gridX[], float gridY[], float **Fn,
				int nX, int nY, float ndegree, float **newFn);
extern void xcRegSplineInt(float xa[], float ya[], float y2a[], 
			   int i, float x, float *y);
extern void xcRegularSplineInt3(float gridX[], float gridY[], float gridZ[], 
				float ***Fn, int nX, int nY, int nZ,
				float ndegree, float ***newFn);
extern void spline3(float x[], int n, float ***y, int j, int k, 
		    float yp1, float ypn, float y2[]);
extern void xcRegSpline3Int(float xa[], float ***ya, int i, int  j, int k, \
			    float y2a[], float x, float *y);



/* --- isoMalloc.c --- */
#ifdef XC_CPP_ISOSURF
   extern void xcFree_PLANEVERTEX(PLANEVERTEX **m);
   extern void xcFree_GRIDVERTEX(GRIDVERTEX ***t);
   extern void xcFree_LINE(LINE *t);
   extern PLANEVERTEX **xcMallocPLANEVERTEX(long nr, long nc);
   extern GRIDVERTEX ***xcMallocGRIDVERTEX(long nr, long nc, long nd);
   extern LINE *xcMallocLINE(long nr);
   extern VERTEX *xcRealloc2xBiggerVertex(VERTEX *vec, int *size);
   extern TRIANGLE *xcRealloc2xBiggerTriangl(TRIANGLE *vec, int *size);
   extern LINE *xcRealloc2xBiggerLINE(LINE *vec, int *size);
#endif


/* --- readPlvertex.c --- */
extern void ReadPlvertex123( int type, int islide );


/* --- isoline.c ---*/
int IsoLine2D(int obj, int cb, int fn, int nx, int ny);

/* --- paraSize.c --- */
extern double DetermineParapipedSized(double *vec1, double *vec2, double *vec3, double *orig);
extern float DetermineParapipedSize(float *vec1, float *vec2, float *vec3, float *orig);
extern float DetermineParalleSize(float *vec1, float *vec2, float *orig);
extern double papipedMaxDiagSize(double v0[3], double v1[3], double v2[3]);


/* --- xcIsoSpaceSel.c --- */
extern void IsoSpaceSel_Parallelogram(void);
extern void IsoSpaceSel_3D(void);
#ifdef _TCL 
extern int XC_IsoSpaceSelCmd(ClientData clientData, Tcl_Interp *interp,
			     int argc, const char *argv[]);

/* --- xcWigner.c --- */
extern int XC_WignerCmd(ClientData clientData,Tcl_Interp *interp,
			int argc, const char *argv[]);
#endif /* _TCL */

/* --- cells.c --- */
extern void CellTypes(void);


/* --- xcTogl.c --- */
#ifdef TOGL_H
/* xcTogl.c */
extern int xcToglCreateFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv);
extern int xcToglReshapeFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv);
extern int xcToglDestroyFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv);
extern int xcToglTimerFunc(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv);
extern int XC_B1MotionCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_ShiftB1MotionCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_B2MotionCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_ButtonReleaseCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
#endif


/* datagrid.c */
extern void CloseGridList(void);
#ifdef XC_CPP_STRUCT
extern struct DATAGRID *FindDataGrid(int index);
#endif
extern int GetNumberOfGridBlocks(void);
extern void NewGridList(int gridtype, FILE *gridFP);
extern FILE *MaybeOpenDataGridFile(char *mode);
extern void CloseDataGridFile(void);
extern void SetDataGridCommentLine(char *line);
extern int ReadDataGrid(FILE *fp, FILE *gridFP, int gridtype, char *ident);
extern int ReadBandGrid(FILE *fp, FILE *gridFP, int gridtype, char *ident);
extern int WriteDataGrid(FILE *fp, int gridtype, const char *ident, int obj);


/* --- xcColorScheme.c --- */
extern void LoadAtmCol(const int flag);
#ifdef _TCL
extern int XC_ColorschemeCmd(ClientData clientData,Tcl_Interp *interp,
			     int argc, const char *argv[]);

/* --- sInfo.c --- */
extern void Set_sInfoArray( Tcl_Interp *interp );


/* --- xcIsoDataGrid.c --- */
extern int XC_IsoDataGridCmd(ClientData clientData, Tcl_Interp *interp,
			     int argc, const char *argv[]);

/* --- xcReadXSF.c --- */
extern int XC_ReadXSFCmd(ClientData clientData, Tcl_Interp *interp,
			 int argc, const char *argv[]);

/* --- xcF3toI4.c --- */
extern int XC_F3toI4Cmd(ClientData clientData, Tcl_Interp *interp,
			int argc, const char *argv[]);


/* --- xcFractCoor.c --- */
extern int XC_FractCoorCmd(ClientData clientData, Tcl_Interp *interp,
			   int argc, const char *argv[]);
#endif /* _TCL */
void GetFractionalCoor(double ivec[][4], float coor[4], float frcoor[]);


/* --- isosurface.c ---*/
#ifdef XC_CPP_ISOSURF 
   extern ISOSURFACE *FindIsoSurf(int index);
   extern void AddToIsoSurfList(ISOSURFACE *g);
   extern int NewIsoSurf(void);
#endif

/* --- clip.c --- */
#ifdef XC_CPP_XYZ
   extern int ClipFacet(XYZ *p, XYZ n, XYZ p0);
   extern int ClipSimple(XYZ *p, XYZ n, XYZ p0);
#endif


/* --- forces.c --- */
extern void MallocForceVectors(void);
#ifdef XC_CPP_VECTOR
extern void SetForceVectorsCoor( ForceVector *fvPtr, double (*fvec)[3], RenderVectors *rvec );
#endif
extern void xcRenderVectorForces(void);


/* --- xcForces.c --- */
#ifdef _TCL
extern int XC_ForcesCmd(ClientData clientData,Tcl_Interp *interp, int argc, const char *argv[]);
#endif
extern void BuildForceVectors(ForceVector *fvPtr);

/* gl2psPrintTogl.c */
#ifdef TOGL_H
/* gl2psPrintTogl.c */
extern int CRY_gl2psPrintToglCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);

/* ppmPrintTogl.c */
extern int CRY_Dump2PpmCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int Togl_DumpToPpmFile(Togl *togl, const char *filename);
extern int Togl_DumpToEpsFile(Togl *togl, const char *filename, int inColor, void (*user_redraw)(const Togl *));

/* togl_ppm.c */
extern int Togl_DumpToPpmFile(struct Togl *togl, const char *filename);
#endif


/* lcasi.c */
extern void cryError(int status, char *function, char *message);
extern unsigned char **cryRegPeriodInterpolator_uc_LCASI2D(int n[2], int degree[2], unsigned char **src);
extern void cryRegPeriodInterpolate_uc_LCASI1D(int N, int DEGREE, unsigned char *src, unsigned char *dst);
extern void cryRegPeriodInterpolate_uc_LCASI2D(int N1, int N2, int DEGREE1, int DEGREE2, unsigned char **src, unsigned char **dst);


/* lcasif.c */
extern float **cryRegPeriodInterpolator_f_LCASI2D(int n[2], int degree[2], float **src);
extern void cryRegPeriodInterpolate_f_LCASI1D(int N, int DEGREE, float *src, float *dst);
extern void cryRegPeriodInterpolate_f_LCASI2D(int N1, int N2, int DEGREE1, int DEGREE2, float **src, float **dst);

/* xcFont.c */
#ifdef TOGL_H
extern int XC_SetFontCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_SetAtomLabelCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_ClearAtomLabelCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_QueryFontCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
#endif
extern void xcFont_PrintString(const char *s);
extern void xcTkFontFreeAll(void);


/* cryNewContext.c */
#ifdef TOGL_H
extern NEW_WIN_CONTEXT *FindWinContextByTogl(struct Togl *togl);
extern NEW_WIN_CONTEXT *FindWinContext(int index);
extern void DestroyWinContext(struct Togl *togl);
extern void AddToWinContext(NEW_WIN_CONTEXT *g);
extern int NewWinContext(void);
extern void cryNewToglInit(struct Togl *togl);
#endif

/* isoInterpolate.c */
extern void isoInterpolate(int degree);

/* crySetProj.c */
#ifdef TOGL_H
extern void crySetProjection(NEW_WIN_CONTEXT *wc, struct Togl *togl);
#endif

/* cryTransform.c */
#ifdef TOGL_H
extern int cryRotateXY(NEW_WIN_CONTEXT *wc, double fiX, double fiY);
#endif


/* xcMolSurf.c */
#ifdef MOLSURF_H
extern MOL_SURF *FindMolSurf(int index);
extern int MolSurfGrid(MOL_SURF *mols, float ***molgrid);
#endif
extern int NewMolSurf(void);
#ifdef ISOSURF_H
extern void MolSurfGridAtomicColor(ISOSURFACE *iso, float mols_d);
#endif

/* crySurf.c */
#if defined (MOLSURF_H) && defined (TOGL_H)
extern int FindSurfaceOfWindow(NEW_WIN_CONTEXT *wc, MOL_SURF *mols);
#endif


/* fs.c */
#if defined (XC_CPP_STRUCT) && defined (MOLSURF_H)
extern int fsReadBand(struct DATAGRID *grid, MOL_SURF *mols);
#endif
#if defined (TOGL_H) && defined (ISOSURF_H)
extern void CropBz(NEW_WIN_CONTEXT *wc, ISOSURFACE *iso);
extern void CropBz_New(NEW_WIN_CONTEXT *wc, ISOSURFACE *iso);
#endif


/* crySurfArgs.c */
#ifdef MOLSURF_H
extern int VerifyReadSurfOpt(const char *option, size_t n, struct READSURF_OPT *opts, char *argv[]);
#endif
#if defined(_TCL) && defined(MOLSURF_H) && defined(TOGL_H)
extern int cryReadSurfARGS(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[], NEW_WIN_CONTEXT *wc, MOL_SURF *mols, struct READSURF_OPT *opt);
extern int ReadFS_ARGS(Tcl_Interp *interp, int argc, const char *argv[], NEW_WIN_CONTEXT *wc, MOL_SURF *mols);
#endif

/* fog.c */
#ifdef __gl_h_
/* fog.c */
#ifdef TOGL_H
extern void xcFog(struct Togl *togl, GLboolean make_fog, GLboolean perspective);
#endif
extern void xcAntiAlias(GLboolean make_antialias);
#endif

/* writeXSF.c */
#ifdef _TCL 
extern int XC_WriteXSFCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern int XC_WriteBandXSFCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
#endif

/* hbonds.c */
extern void make_H_Bonds(void);
extern void xcRenderHbonds3D(void);
/*extern VEC3d *xcReallocVEC3d(VEC3d *vec, int old_nmemb, int new_nmemb);*/

/* xcHBonds.c */
#ifdef _TCL 
int XC_HBondsCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
#endif

/* linear.c */
extern float ***cryGeneralGridRegPeriodInterpolator_f_Linear3D(int n[3], int degree[3], float ***src);

/* realTimeMovie.c */
#ifdef TOGL_H
extern int CRY_RealTimeMovieCb(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]);
extern void createMoviePPMFrame(struct Togl *togl);
#endif

/* fft3d.c */
extern float ***general_grid_fft_interpolator_tensor3f(int n[3], int degree[3], float ***src);
extern double *fft3d_interpolator(int n1, int n2, int n3, int n1new, int n2new, int n3new, double *func);
extern void fft3d_interpolate(int n1, int n2, int n3, int n1new, int n2new, int n3new, double *func, double *result);

/* shank3d.c */
extern float ***general_grid_shankland_interpolator_tensor3f(int n[3], int degree[3], float ***src);
extern double *shankland3d_interpolator(int n1, int n2, int n3, int n1new, int n2new, int n3new, double *func, double (*roughnessFunc)(double m2, double sigma), double log10_damp);
extern void shankland3d_interpolate(int n1, int n2, int n3, int n1new, int n2new, int n3new, double *func, double (*roughnessFunc)(double m2, double sigma), double log10_damp, double *result);
