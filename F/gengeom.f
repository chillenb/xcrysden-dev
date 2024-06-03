      PROGRAM GENGEOM
C     *******************************************
C     ***  GENERAL GEOMETRIC PURPOSE PROGRAM  ***
C     *******************************************

c     =======================================================================
c     USAGE: 
c     gengeom  MODE1  MODE2  MODE3  IGRP  NXDIR  NYDIR  NZDIR  OUTPUT  INPUT
c        0       1      2      3      4     5      6      7      8       9
c     =======================================================================

c     program read data from ftn34, but if INPUT ARGUMENT is present it 
c     rather read data from gengeom_formated INPUT file instead

C     NXDIR,NYDIR,NZDIR --> n. of repeatitions in each directions

C     *****
C     IGR - SEQUENTIAL GROUP NUMBER (USED ONLY FOR CRYSTALS)
C     posebej je treba paziti na R/H !!!!!
C     *****

c     how the vectors are packed to matrix::
c     --------------------------------------
c     DIRECT VECTORS:: 
c     ----------------
c     VECTOR1 GOES TO FIRST COLOUMN, VETCOR2 TO SECOND & VECTOR3 TO 
c     THIRD COLOUMN;;; DV(COL,ROW);  vec(# ,x/y/z)
c                                    =============
c           |11   21   31|          |x1  x2  x3|
c     VEC = |12   22   32|  => DV = |y1  y2  y3|
c           |13   23   33|          |z1  z2  z3|
c
c     RECIPROCAL VECTORS:: (in transpose manner)
c     --------------------
c     VECTOR1 GOES TO FIRST ROW, VETCOR2 TO SECOND & VECTOR3 TO 
c     THIRD ROW;;;                   vec(x/y/z, #)
c                                    =============
c            |11   21   31|          |x1* y1* z1*|
c     IVEC = |12   22   32|  =>IDV = |x2* y2* z2*|
c            |13   23   33|          |x3* y3* z3*|
c
c     so:                 |1  0  0|
c          VEC x IVEC ==  |0  1  0|
c                         |0  0  1|

      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      CHARACTER*256 MODE1, MODE2, MODE3,
     *     ANX, ANY, ANZ, FILE2, INPUT
      REAL*8 PC(3,1),AC(3,2),BC(3,2),CC(3,2),FC(3,4),IC(3,2),RC(3,3),
     *     HC(3,3),CSTMC(3,4),P2A(3,3),P2B(3,3),P2C(3,3),P2F(3,3),
     *     P2I(3,3),R2H(3,3),P2AI(3,3),P2BI(3,3),P2CI(3,3),P2FI(3,3),
     *     P2II(3,3),R2HI(3,3),P2H(3,3),P2HI(3,3),
     *     DV(3,3),SOP(3,3,48),TRX(48),TRY(48),TRZ(48),
     *     DVC(3,3),CON13,CON23,HTTC(3,13),RTTC(3,13),MINTOL,
     *     DVEC(3,3),TMPVEC(3,3),X(NAC),Y(NAC),Z(NAC),
     *     IDVEC(3,3),IDV(3,3),IDVC(3,3),
     *     RCR(3,3),RTTCR(3,13),S_AC(3,2),S_BC(3,2),S_CC(3,2),
     *     S_FC(3,4),S_IC(3,2),S_RC(3,3),S_HC(3,3)
      INTEGER C2I, NAT(NAC)
      LOGICAL L, lreverse, IsReverse, notfound, Lprimvec

      include 'mode.inc'

      COMMON/FTN34/ SOP,TRX,TRY,TRZ,DV,IDIM,IGROUP,NSYMOP
      COMMON/APOS/ AC,BC,CC,FC,IC,HC,RC,DVC,CSTMC
      COMMON/DIR/ NXDIR,NYDIR,NZDIR
      COMMON/IVEC/ IDVEC
      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR
      COMMON/BOHRR/ BOHR,IMODE3
      COMMON/KEYWORD/ notfound, Lprimvec

      PARAMETER (CON13=1.0d0/3.0d0, CON23=2.0d0/3.0d0, MINTOL=1.d-6,
     $     CON13N=-1.0d0/3.0d0, CON23N=-2.0d0/3.0d0,
     $     TOLMIN=1.d-6)

C     VARIABLE NAMES:
C     DV{1|2|3}.............3 direct vectors in PRIMITIV cell
C     DVC{1|2|3}............3 direct vectors in CONVENTIONAL cell
C     TR{X|Y|Z}.............X,Y,Z COMPONENT OF TRANSL. DUE TO POINT SYM. OPER.
C     IDIM..................DIMENSION OF THE SYSTEM
C     IGROUP................GROUP NUMBER (LOOK BELOW)
C     SOP...................SYM. OPERATORS (MATRICES)

C     ***********************************************
C     TRANSFORMATION of PRIMITIV CELL to CONVENTIONAL
C     ***********************************************
C     ===================
C     IGROUP DEFINITIONS:
C     ===================
C     GROUP NUMBER USED IN ftn34
C     N.  symbol   N.of atoms   Variable
C     1.....P   -->    1           PC
C     2.....A   -->    2           AC
C     3.....B   -->    2           BC
C     4.....C   -->    2           CC
C     5.....F   -->    4           FC
C     6.....I   -->    2           IC
C     7.....R   -->    3           RC --> RTTC
C     REMARK: ALL HEXAGONAL CELLS ARE OF TYPE 1 IN FTN34 !!!!! 
C     SO I'LL USE 8 FOR HEXAGONAL CELL
C     8.....H   -->    3           HC --> HTTC
C     for TRIGONAL CELL that are not RHOMOHEDRAL, I will use additional flag
C     9.....TRIGONAL_NOT_RHOMOHEDRAL
C     10....custom; used just when calling CONVATOM to tell to read CSTMC
C     ----------------------------------------
C     ----------------------------------------
C     TRANSFORMATION MATRICES:
C     ------------------------
C     EXAMPLE:
C     P2C means: PRIMITIV to C-centered
C     P2C and the rest refer to coor. of point, so
C     P2CI and the rest refer to transf. of vectors
C     P2CI means "inverse of P2C"

      DATA PC/
     1     0.0d0, 0.0d0, 0.0d0/     
C     CONVENTIONAL CELL; COORD. OF LATTICE POINTS WITHIN CELL
      DATA S_AC/
     1     0.0d0, 0.0d0, 0.0d0,
     2     0.0d0, 0.5d0, 0.5d0/
      DATA S_BC/
     1     0.0d0, 0.0d0, 0.0d0,
     2     0.5d0, 0.0d0, 0.5d0/
      DATA S_CC/
     1     0.0d0, 0.0d0, 0.0d0,
     2     0.5d0, 0.5d0, 0.0d0/
      DATA S_FC/
     1     0.0d0, 0.0d0, 0.0d0,
     2     0.5d0, 0.5d0, 0.0d0,
     3     0.0d0, 0.5d0, 0.5d0,
     4     0.5d0, 0.0d0, 0.5d0/
      DATA S_IC/
     1     0.0d0, 0.0d0, 0.0d0,
     2     0.5d0, 0.5d0, 0.5d0/      
      DATA S_HC/
     1     0.0d0, 0.0d0, 0.0d0,
     2     CON23, CON13, 0.0d0,
     3     CON13, CON23, 0.0d0/
C     RHOMBOHEDRALY CENTERED (OBVERSE SETTINGS) - hexagonal axes
C     CRYSTALXX uses this     ^^^^^^^
      DATA S_RC/
     1     0.0d0, 0.0d0, 0.0d0,
     2     CON23, CON13, CON13,
     3     CON13, CON23, CON23/      

C     RHOMBOHEDRALY CENTERED (REVERSE SETTINGS) - hexagonal axes
C     WIEN uses this          ^^^^^^^
      DATA RCR/
     1     0.0d0, 0.0d0, 0.0d0,
     2     CON13, CON23, CON13,
     3     CON23, CON13, CON23/      

C     HEXAGONALY TRIPLE CENTERED CELL - HEXAGONAL SHAPE 
      DATA HTTC/
     1     1.0d0, 0.0d0, 0.0d0,   1.0d0, 1.0d0, 0.0d0,
     2     0.0d0, 1.0d0, 0.0d0,  -1.0d0, 0.0d0, 0.0d0,
     3    -1.0d0,-1.0d0, 0.0d0,   0.0d0,-1.0d0, 0.0d0,
     4     0.0d0, 0.0d0, 0.0d0,
     5     CON23, CON13, 0.0d0,   CON13, CON23, 0.0d0,
     6    CON13N, CON13, 0.0d0,  CON23N,CON13N, 0.0d0,
     7    CON13N,CON23N, 0.0d0,   CON13,CON13N, 0.0d0/
C     RHOMBOHEDRAL - described as TRIPLE HEXAGONAL CELL - HEXAGONAL SHAPE
c     OBVERSE SETTINGS
c     ^^^^^^^
      DATA RTTC/
     1     1.0d0, 0.0d0, 0.0d0,   1.0d0, 1.0d0, 0.0d0,
     2     0.0d0, 1.0d0, 0.0d0,  -1.0d0, 0.0d0, 0.0d0,
     3    -1.0d0,-1.0d0, 0.0d0,   0.0d0,-1.0d0, 0.0d0,
     4     0.0d0, 0.0d0, 0.0d0,
     5     CON23, CON13, CON13,   CON13, CON23, CON23,
     6    CON13N, CON13, CON13,  CON23N,CON13N, CON23,
     7    CON13N,CON23N, CON13,   CON13,CON13N, CON23/
c     REVERSE SETTINGS
c     ^^^^^^^
      DATA RTTCR/
     1     1.0d0, 0.0d0, 0.0d0,   1.0d0, 1.0d0, 0.0d0,
     2     0.0d0, 1.0d0, 0.0d0,  -1.0d0, 0.0d0, 0.0d0,
     3    -1.0d0,-1.0d0, 0.0d0,   0.0d0,-1.0d0, 0.0d0,
     4     0.0d0, 0.0d0, 0.0d0,
     5     CON13, CON23, CON13,   CON23, CON13, CON23,
     6    CON23N, CON23, CON13,  CON13N,CON23N, CON23,
     7    CON23N,CON13N, CON13,   CON23,CON23N, CON23/

C     TRANSFORMATION MATRICES IN CRYSTAL95
      DATA P2A/
     1     1.0d0, 0.0d0, 0.0d0,
     2     0.0d0, 0.5d0,-0.5d0,
     3     0.0d0, 0.5d0, 0.5d0/
      DATA P2AI/
     1     1.0d0, 0.0d0, 0.0d0,
     2     0.0d0, 1.0d0, 1.0d0,
     3     0.0d0,-1.0d0, 1.0d0/

      DATA P2B/
     1     0.5d0, 0.0d0, 0.5d0,
     2     0.0d0, 1.0d0, 0.0d0,
     3    -0.5d0, 0.0d0, 0.5d0/
      DATA P2BI/
     1     1.0d0, 0.0d0,-1.0d0,
     2     0.0d0, 1.0d0, 0.0d0,
     3     1.0d0, 0.0d0, 1.0d0/

      DATA P2C/
     1     0.5d0,-0.5d0, 0.0d0,
     2     0.5d0, 0.5d0, 0.0d0,
     3     0.0d0, 0.0d0, 1.0d0/
      DATA P2CI/
     1     1.0d0, 1.0d0, 0.0d0,
     2    -1.0d0, 1.0d0, 0.0d0,
     3     0.0d0, 0.0d0, 1.0d0/

      DATA P2F/
     1     0.0d0, 0.5d0, 0.5d0,
     2     0.5d0, 0.0d0, 0.5d0,
     3     0.5d0, 0.5d0, 0.0d0/
      DATA P2FI/
     1    -1.0d0, 1.0d0, 1.0d0,
     2     1.0d0,-1.0d0, 1.0d0,
     3     1.0d0, 1.0d0,-1.0d0/

      DATA P2I/
     1    -0.5d0, 0.5d0, 0.5d0,
     2     0.5d0,-0.5d0, 0.5d0,
     3     0.5d0, 0.5d0,-0.5d0/
      DATA P2II/
     1     0.0d0, 1.0d0, 1.0d0,
     2     1.0d0, 0.0d0, 1.0d0,
     3     1.0d0, 1.0d0, 0.0d0/

      DATA R2H/
     1     CON23,CON13N,CON13N,
     2     CON13, CON13,CON23N,
     3     CON13, CON13, CON13/
      DATA R2HI/
     1     1.0d0, 0.0d0, 1.0d0,
     2    -1.0d0, 1.0d0, 1.0d0,
     3     0.0d0,-1.0d0, 1.0d0/
C     HEXAGONAL CELL P --> TRIPLE HEXAGONAL CELL H1;
C     INTERNATIONAL TABLES OF CRYSTALOGRAPHY, 1993; TABLE 5-1
      DATA P2H/
     1     CON23,CON13N, 0.0d0,
     2     CON13, CON13, 0.0d0,
     3     0.0d0, 0.0d0, 1.0d0/     
      DATA P2HI/
     1     1.0d0, 1.0d0, 0.0d0,
     2    -1.0d0, 2.0d0, 0.0d0,
     3     0.0d0, 0.0d0, 1.0d0/      
C     === END_OF_DATA_BLOCK ===

      call MATCOPY(S_AC,AC,3,2)
      call MATCOPY(S_BC,BC,3,2)
      call MATCOPY(S_CC,CC,3,2)
      call MATCOPY(S_FC,FC,3,4)
      call MATCOPY(S_IC,IC,3,2)
      call MATCOPY(S_RC,RC,3,3)
      call MATCOPY(S_HC,HC,3,3)

C     CONVERSION FROM BOHRS TO ANGS.
      BOHR=0.529177d0

      NCSTMC=0

c     number of command line arguments must not be lower than 8
      NARG=IARGC()
      IF(NARG.LT.8)THEN
         WRITE(*,*) 
     1        'usage:  gengeom <mode1> <mode2> <mode3> <igrp> <nxdir> <n
     2ydir> <nzdir> <output> [input]'
         STOP
      ENDIF

C     ***********************************
C     *** READ MODE & N?DIR ARGUMENTS ***
C     ***********************************
C     look in 'mode.inc' for definitions !!!!
      CALL GETARG(IARG_MODE1,MODE1)
      CALL GETARG(IARG_MODE2,MODE2)
      CALL GETARG(IARG_MODE3,MODE3)
      CALL GETARG(IARG_NXDIR,ANX)
      CALL GETARG(IARG_NYDIR,ANY)
      CALL GETARG(IARG_NZDIR,ANZ)
      IMODE1=C2I(MODE1)
      IMODE2=C2I(MODE2)
      IMODE3=C2I(MODE3)
      NXDIR =C2I(ANX)
      NYDIR =C2I(ANY)
      NZDIR =C2I(ANZ)
c     determine aproximate number of atoms (COUNT JUST ONE ATOM PER CELL)
      NATOMS=(NXDIR+2) * (NYDIR+2) *(NZDIR+1)  !overestimated, thats good

      CALL GETARG(IARG_OUTPUT,FILE2)
      IFILE2=INDEX(FILE2,' ')
      OPEN(UNIT=12,FILE=FILE2(1:IFILE2-1),STATUS='UNKNOWN')
C      file2='tmp.tmp'
C      OPEN(UNIT=12,FILE='tmp.tmp',STATUS='UNKNOWN')

c     --- DEBUG_BEGIN ---
c      print *,'IMODE1=',IMODE1,'; IMODE2=',IMODE2,'; IMODE3=',IMODE3
c      print *,'NXDIR=',NXDIR,'; NYDIR=',NYDIR,'; NZDIR=',NZDIR
c      print *,'IFILE2=',IFILE2,'; OUTPUT::',FILE2(1:IFILE2-1)
c     --- DEBUG_END ---

c     *********************************
c     take care of imode3; imode3 comprise several modes, that are based on
c     following criteria: 1258, 1 is mode3, 2 is mode4, 5 is mode5, 8 is mode6
c     so far only MODE3 & MODE4 are used!!!!
      IMODE4=MOD(IMODE3,10)
      IMODE3=(IMODE3-IMODE4)/10

C     *********************************
C     ***     write INFO record     ***
      WRITE(12,'(1x,a)') 'INFO'
      WRITE(12,'(a,3x,i4,1x,i4,1x,i4)') 'nunit',NXDIR,NYDIR,NZDIR
      IF(IMODE2.EQ.M2_CELL) WRITE(12,'(a,3x,a)') 'unit','cell'
      IF(IMODE2.EQ.M2_TR_ASYM_UNIT)
     1     WRITE(12,'(a,3x,a)') 'unit','tr_asym'
      IF(IMODE1.EQ.M1_PRIM .OR. IMODE1.EQ.M1_PRIM3) THEN
         WRITE(12,'(a,3x,a)') 'celltype', 'primcell'
      ELSE
         WRITE(12,'(a,3x,a)') 'celltype', 'convcell'
      ENDIF
      IF(IMODE1.EQ.M1_HEXA_SHAPE .OR. IMODE1.EQ.M1_PRIM3) THEN
         WRITE(12,'(a,3x,a)') 'shape','hexagonal'
      ELSE
         WRITE(12,'(a,3x,a)') 'shape','parapipedal'
      ENDIF
      WRITE(12,'(1x,a)') 'END_INFO'
C     *********************************

      IGROUP_input=0
      IF(NARG.EQ.9)THEN
C     *********************************
C     *** is INPUT ARGUMENT PRESENT ***
C     ***      read from INPUT      ***
C     *********************************
         IDUM1=0
         IDUM2=0
         CALL GETARG(IARG_INPUT,INPUT)
         I_INPUT=INDEX(INPUT,' ')
         OPEN(UNIT=11,FILE=INPUT(1:I_INPUT),STATUS='OLD')
c     READF1 has two common's that are:
c     /MULTAT/ X,Y,Z,NAT,DVEC,NATR
c     /IVEC/ IDVEC

c     backward compatibility
         CALL READF1('DIM-GROUP',IGROUP,IDIM)
c     --
         CALL READF1('MOLECULE',IGROUP,IDIM)
         CALL READF1('POLYMER',IGROUP,IDIM)
         CALL READF1('SLAB',IGROUP,IDIM)
         CALL READF1('CRYSTAL',IGROUP,IDIM)

         CALL READF1('PRIMVEC',IGROUP,IDIM)
         if (notfound) then
            Lprimvec=.false.         
         else
            Lprimvec=.true.
            CALL MATCOPY(DVEC,DV,3,3)
         endif

         CALL READF1('CONVVEC',IGROUP,IDIM)
         if ((.not.notfound) .and. (.not.Lprimvec) ) then
            call Matcopy(dvec,dv,3,3)
            write(12,'(1x,a)') 'PRIMVEC'
            write(12,'(3(1x,f15.10)/3(1x,f15.10)/3(1x,f15.10))') 
     1           ((dvec(j,i),i=1,3),j=1,3)
         else if (notfound .and. (.not.Lprimvec)) then
            WRITE(6,*) 'PRIMVEC and CONVEC keywords not found !!!'
            STOP            
         endif
         CALL MATCOPY(DVEC,DVC,3,3)

         if (Lprimvec) then
            CALL READF1('PRIMCOORD',IDUM1,IDUM2) !here NAT is read
         else 
c     assuming PRIMCOORD == CONVCOORD
            CALL READF1('CONVCOORD',IDUM1, IDUM2)
         endif

         CALL GetCenteredPoints(DV, DVC, IDIM, CSTMC, NCSTMC)
         IGROUP_input=IGROUP    !to save original igroup read from input
         IGROUP=10              !10 means custom (input file was read)

c     DEBUG_BEGIN
c         if (IGROUP_input.EQ.7)
c     $        print *,'DEBUG: IsReverse = ', IsReverse(CSTMC,RCR)
c         print *,'DEBUG: NCSTMC::', NCSTMC
c         do i=1,4
c            print *,'DEBUG: CSTMC::',(CSTMC(j,i),j=1,3)
c         enddo
c     DEBUG_END

         if (ncstmc.gt.4) then
c     Largest number of special centered points has FCC == 4; if NCSTMC exceed
c     4, something is wrong
            STOP 'more then four centered points; conventional vectors 
     $           are probably badly chosen'
         ENDIF
         CLOSE(11)
      ELSE
C     *************************************************
C     *****  READ   FTN34;  all modes need this *******
C     *************************************************
         NATR=0
c     READFTN34 has one common:
C     /FTN34/ IDIM,IGROUP,DV,NSYMOP,SOP,TRX,TRY,TRZ
         CALL READFTN34(NATR, IMODE4) !NATR is read

C     ****************************
C     get the CONVENTIONAL VECTORS
C     ****************************
         call ZEROMAT(DVC, 3, 3)      
C     P->P
         IF(IGROUP.EQ.1) THEN
            DO J=1,3
               DO I=1,3
                  DVC(I,J)=DV(I,J)
               ENDDO
            ENDDO
         ENDIF
C     P->A
         IF(IGROUP.EQ.2) CALL MATMULT(DV,3,3,P2AI,3,DVC)
C     P->B
         IF(IGROUP.EQ.3) CALL MATMULT(DV,3,3,P2BI,3,DVC)
C     P->C                       
         IF(IGROUP.EQ.4) CALL MATMULT(DV,3,3,P2CI,3,DVC)
C     P->F            
         IF(IGROUP.EQ.5) CALL MATMULT(DV,3,3,P2FI,3,DVC)
C     P->I                       
         IF(IGROUP.EQ.6) CALL MATMULT(DV,3,3,P2II,3,DVC)
C     PR->H; primitiv-rhombohedrall to triple-hexagonal R cell
         IF(IGROUP.EQ.7) CALL MATMULT(DV,3,3,R2HI,3,DVC)
C     PH->H; primitiv hexagonal --> triple hexagonal H cell
         IF(IGROUP.EQ.8 .OR. IGROUP.EQ.9) 
     1        CALL MATMULT(DV,3,3,P2HI,3,DVC)

C     WRITE CONVENTIONAL DIRECT VECTORS 
         CALL WVEC('CONVVEC',DVC)
      ENDIF

c     ///////////////////////////////////////////////////////
c     // below is common to ftn34 and XCrySDen format file //
c     ///////////////////////////////////////////////////////

C     NOW CALCULATE POSITIONS OF ATOMS NEEDED WITHIN CONV.CELL
      CALL CONVATOM(IGROUP,NARG,NCSTMC,IMODE4)  !primcoord & convcoord are writen to UNIT12

C     make reciprocal vectors (TRANSPOSE form)
      CALL ZeroMat(IDV, 3, 3)
      CALL ZeroMat(IDVC, 3, 3)
c     initialise tmpvec matrix; tmpvec is dummy anyway
      call InvertSqMat123(DV, IDV, IDIM)
      call InvertSqMat123(DVC, IDVC, IDIM)

c     write PRIMITIV & CONVENTIONAL RECIPROCAL VECTORS in normal not
c     transpose form
      CALL WIVEC('RECIP-PRIMVEC',IDV) 
      CALL WIVEC('RECIP-CONVVEC',IDVC)

c     write WIGNER-SEITZ CELL
      IF(IDIM.EQ.3 .AND. IMODE1.NE.M1_INFO) THEN
         CALL WignerSeitz(DV,   'WIGNER-SEITZ-PRIMCELL')
         CALL WignerSeitz(DVC,  'WIGNER-SEITZ-CONVCELL')
         Call MatTranspose(IDV, TMPVEC, 3, 3)
         CALL WignerSeitz(TMPVEC,  'BRILLOUIN-ZONE-PRIMCELL')
         Call MatTranspose(IDVC, TMPVEC, 3, 3)            
         CALL WignerSeitz(TMPVEC, 'BRILLOUIN-ZONE-CONVCELL')
      ELSEIF(IDIM.EQ.2 .AND. IMODE1.NE.M1_INFO) THEN
         CALL WignerSeitz2D(DV,   'WIGNER-SEITZ-PRIMCELL')
         Call MatTranspose(IDV, TMPVEC, 3, 3)
         CALL WignerSeitz2D(TMPVEC,  'BRILLOUIN-ZONE-PRIMCELL')
      ELSEIF(IDIM.EQ.1 .AND. IMODE1.NE.M1_INFO) THEN
         CALL WignerSeitz1D(DV,   'WIGNER-SEITZ-PRIMCELL')
         CALL WignerSeitz1D(IDV,  'BRILLOUIN-ZONE-PRIMCELL')
      ENDIF            
      
C     ****************************
C     deside upon IMODE1
      IF(IMODE1.EQ.M1_INFO) GOTO 1111 !info mode
      GOTO(150,250,250,250) IMODE1
C     ============================
 150  CONTINUE
C     *******************************************************
C     *    CALCULATE COORDINATES FOR PRIMITIV CELL(S);      *
C     *                              ^^^^^^^^               *
C     *                   IMODE1 = M1_PRIM                  *
C     *******************************************************

c     ASSING primitiv-vectors to DVEC & IDVEC
      CALL MATCOPY(DV,DVEC,3,3)  !dvec is used in MULTATOM; common/MULTAT/
      CALL MATCOPY(IDV,IDVEC,3,3)

C     IF IDIM=0 JUST PRINT THE COORDINATES
      IF(IDIM.EQ.0)THEN
         CALL MOLECULE
      ELSE
c         print *,'NATOMS',NATOMS
         CALL MULTATOM(PC,1,1,NATOMS,.false.,IDIM,IMODE2)
      ENDIF

      GOTO 1111
 
 250  CONTINUE
C     *********************************************************
C     **  calculate coordinates for CONVENTIONAL cell(s) &&  **
C     **  also M1_PRIM3; look below ^^^^^^^^^^^^             **
C     *********************************************************
c
C     **********************************************************
C     for TRIGONAL/HEXAGONAL/RHOMBOHEDRAL SYSTEMS: 
C     **  IMODE1=M1_CONV ........ PARAPIPEDAL SHAPED TRIPLE hexagonal CELL
C     **  IMODE1=M1_HEXA_SHAPE .. HEXAGONAL SHAPED TRIPLE hexagonal CELL
c
C     for TRIGONAL_NOT_RHOMBO/HEXAGONAL systems
C     **  IMODE1=M1_PRIM3 --> HEXAGONAL SHAPED PRIMITIV CELL
c                                              ^^^^^^^^
C     **********************************************************

C     *** so far, only crystals can be converted from prim. to conv. cell ***
      IF(IDIM.NE.3)
     1     STOP 'conversion to conventional cell for non-crystal'

c     ASSING conventional-direct-vectors to DVEC
      CALL MATCOPY(DVC,DVEC,3,3) !dvec is used in MULTATOM & MULTHEXA
      CALL MATCOPY(IDVC,IDVEC,3,3)

C     HEXAGONAL SHAPED PRIMITIV RHOMOHEDRAL CELL
c      IF(IMODE1.EQ.M1_PRIM3 .AND. IGROUP.EQ.7) !maybe turn this off
c     1     CALL MULTATOM(RC,3,8,NATOMS*3,.true.,3,IMODE2)       

C     HEXAGONAL SHAPED PRIMITIV TRIGONAL_NOT_RHOMBO/HEXAGONAL CELL
C     t.k.: maybe TRIGONAL is not allowed for this
      IF(IMODE1.EQ.M1_PRIM3 .AND.
     $     (IGROUP.EQ.8 .OR. IGROUP.EQ.9 .OR.
     $     IGROUP_input.EQ.8 .OR. IGROUP_input.EQ.9)) then
         CALL MULTATOM(HC,3,8,NATOMS*3,.true.,3,IMODE2)
         goto 1111
      ENDIF

C     WHEATER (NOT HEXA/RHOMBO/TRIGONAL_NOT_RHOMBO)
      L=.FALSE.
      IF(IGROUP.EQ.1) CALL MULTATOM(PC,1,IGROUP,NATOMS,L,3,IMODE2)
      IF(IGROUP.EQ.2) CALL MULTATOM(AC,2,IGROUP,NATOMS*2,L,3,IMODE2)
      IF(IGROUP.EQ.3) CALL MULTATOM(BC,2,IGROUP,NATOMS*2,L,3,IMODE2)
      IF(IGROUP.EQ.4) CALL MULTATOM(CC,2,IGROUP,NATOMS*2,L,3,IMODE2)
      IF(IGROUP.EQ.5) CALL MULTATOM(FC,4,IGROUP,NATOMS*4,L,3,IMODE2)
      IF(IGROUP.EQ.6) CALL MULTATOM(IC,2,IGROUP,NATOMS*2,L,3,IMODE2)
      IF(IGROUP.EQ.10) then
         if(IGROUP_input.EQ.7 .and. IMODE1.EQ.M1_HEXA_SHAPE) then
c     rhombohedral
            lreverse=IsReverse(CSTMC, RCR)
            if(lreverse) then
c     reverse settings
               CALL MULTHEXA(RTTCR,IGROUP,NATOMS*13) 
               goto 1111
            endif
         elseif(IGROUP_input.GE.8 .and. IMODE1.EQ.M1_HEXA_SHAPE) then
c     hexagonal
            IGROUP=IGROUP_input
         else
            CALL MULTATOM(CSTMC,NCSTMC,IGROUP,NATOMS*NCSTMC,L,3,IMODE2)
            goto 1111
         endif
      ENDIF

C     TRIPLE HEXA/RHOMBO
C     M1_HEXA_SHAPE & M2_TR_ASYM_UNIT together are not allowed
      IF(IMODE1.EQ.M1_HEXA_SHAPE .AND. IMODE2.EQ.M2_TR_ASYM_UNIT)
     *  STOP 'M1_HEXA_SHAPE & M2_TR_ASYM_UNIT can not be used together'
      IF(IMODE1.EQ.M1_CONV .AND. IGROUP.EQ.7) 
     *     CALL MULTATOM(RC,3,7,NATOMS*3,L,3,IMODE2)
      IF(IMODE1.EQ.M1_CONV .AND. (IGROUP.EQ.8 .OR. IGROUP.EQ.9))
     *     CALL MULTATOM(HC,3,8,NATOMS*3,L,3,IMODE2) 
      IF(IMODE1.EQ.M1_HEXA_SHAPE .AND. IGROUP.EQ.7)
     *     CALL MULTHEXA(RTTC,IGROUP,NATOMS*13) 
      IF(IMODE1.EQ.M1_HEXA_SHAPE .AND. (IGROUP.EQ.8 .OR. IGROUP.EQ.9))
     *     CALL MULTHEXA(HTTC,IGROUP,NATOMS*13) 
 
 1111 CONTINUE
      CLOSE(12)
      END


c     INTEGER FUNCTION DETGROUP(IGR)
C     FUNCTION DETERMINES IF GROUP IS OF HEXAGONAL TYPE
C     8 FOR HEXAGONAL TYPE
c     DETGROUP=1
c     IF(IGR.GE.143.AND.IGR.LE.194)DETGROUP=8
c     RETURN
c     END


      SUBROUTINE READFTN34(NATR, IMODE4)
C     THIS SUBROUTINE READ THE ftn34
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*80 CDUM,AGRP
      REAL*8 DV(3,3),SOP(3,3,48),TRX(48),TRY(48),TRZ(48)
      INTEGER C2I
      COMMON/FTN34/ SOP,TRX,TRY,TRZ,DV,IDIM,IGROUP,NSYMOP
      COMMON/BOHRR/ BOHR,IMODE3

      include 'mode.inc'

c     --- DEBUG_BEGIN ---
c      print *,'in READFTN34; IMODE4=',IMODE4
c     --- DEBUG_END ---

      if (IMODE4 .eq. M4_CR95) then
         READ(34,*) CDUM
         READ(34,*) IDIM,IGROUP
      elseif (IMODE4 .eq. M4_CR98) then
         READ(34,*) IDIM,IGROUP,ITYPE
      else
         STOP 'unknown IMODE4, must be M4_CR95 or M4_CR98'
      endif

C     TAKE CARE IF CELL OF TYPE 1 IS HEXAGONAL TYPE
      IF(IGROUP.EQ.1) THEN
         CALL GETARG(IARG_IGRP,AGRP)
c     agrp is '8' for HEXAGONAL systems and '9' for TRIGONAL_NOT_RHOMBOHEDRAL
         IGROUP=C2I(AGRP)  
      ENDIF
      WRITE(12,*) 'DIM-GROUP'
      WRITE(12,*) IDIM,IGROUP
C     VECTOR1 GOES TO FIRST COLOUMN, VETCOR2 TO SECOND & VECTOR3 TO 
C     THIRD COLOUMN;;; DV(COL,ROW)
C           |11   21   31|          |x1  x2  x3|
C     VEC = |12   22   32|  => DV = |y1  y2  y3|
C           |13   23   33|          |z1  z2  z3|
      READ(34,*) ((DV(I,J),J=1,3),I=1,3)
      IF(IMODE3.EQ.M3_BOHR)THEN !convert from bohrs to angstroms
         DO I=1,3
            DO J=1,3
               DV(J,I) = BOHR * DV(J,I)
            ENDDO
         ENDDO
      ENDIF
      if(idim.lt.3)dv(3,3) = 1.0d0
      if(idim.lt.2)dv(2,2) = 1.0d0
      if(idim.lt.1)dv(1,1) = 1.0d0
      WRITE(12,*) 'PRIMVEC'
      
      DO I=1,3
         WRITE(12,'(3(1x,f15.10)/3(1x,f15.10)/3(1x,f15.10))')
     $        (DV(I,J),J=1,3)
      ENDDO
      READ(34,*) NSYMOP
      DO I=1,NSYMOP
C     READ SYMM. OPERATOR
         DO J=1,3
            READ(34,*) (SOP(JJ,J,I),JJ=1,3)
         ENDDO
C     READ CARTESIAN COMP. OF THE TRANSALTIUON DUE TO SYM. OPER.
         READ(34,*) TRX(I),TRY(I),TRZ(I)
         IF(IMODE3.EQ.M3_BOHR)THEN !convert from bohrs to angstroms
            trx(i)=trx(i)*BOHR
            try(i)=try(i)*BOHR
            trz(i)=trz(i)*BOHR
         ENDIF
      ENDDO

C     READ NATR;
      READ(34,*) NATR
C     NAT,X,Y,Z WILL ARE READ in CONVATOM subroutine !!!!!!

c      print *,'READFTN34_END'

      RETURN
      END

      SUBROUTINE CONVATOM(IGROUP,NARG,NCSTMC,IMODE4)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      CHARACTER*80 FILE2
      INTEGER NAT(NAC)
      REAL*8 AC(3,2),BC(3,2),CC(3,2),FC(3,4),IC(3,2),RC(3,3),HC(3,3),
     $     CSTMC(3,4),
     *     DVC(3,3),SOP(3,3,48),TRX(48),TRY(48),TRZ(48),
     *     XC(NAC,4),YC(NAC,4),ZC(NAC,4),
     *     X(NAC),Y(NAC),Z(NAC),DVEC(3,3)
      COMMON/APOS/ AC,BC,CC,FC,IC,HC,RC,DVC,CSTMC
      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR
      COMMON/BOHRR/ BOHR,IMODE3
 
      include 'mode.inc'

c     DEBUG_BEGIN
c      print *,'DEBUG: CONVATOM'
c      print *,'DEBUG: NCSTMC::', NCSTMC
c      do i=1,4
c         print *,'DEBUG: CSTMC::',(CSTMC(j,i),j=1,3)
c      enddo
c     DEBUG_END

C     READ NAT,X,Y,Z FROM UNIT 11 (FILE1) & PREPARE ATOMS FOR DIFFERENT 
C     CENTERED CELL

c      print *,'IN CONVATOM'

      DO I=1,NATR
         if (narg.eq.8) then
            READ(34,*) NAT(I),X(I),Y(I),Z(I)
            IF(IMODE3.EQ.M3_BOHR) THEN
               X(I) = X(I) * BOHR
               Y(I) = Y(I) * BOHR
               Z(I) = Z(I) * BOHR
            ENDIF
         endif
      enddo

      if (IMODE4 .eq. M4_CR98) then
c     generate all equivalent atoms from non-equivalent
         call GenEquPos    !GenEquPos is still bugy
      endif

c     write PRIMCOOR just if we have read data from unit34
      if (narg.eq.8) then
         WRITE(12,*) 'PRIMCOORD'
         WRITE(12,*) NATR,' 1'
      endif
      
      DO I=1,NATR
         if(narg.eq.8) WRITE(12,'(i3,3x,3(f15.10,2x))')
     $        NAT(I),X(I),Y(I),Z(I)
c     print *,'i;nat,x,y,z>>',i,NAT(I),X(I),Y(I),Z(I)
C     ALL DIFFERENT GROUPS HAVE (0,0,0)
         XC(I,1)=X(I)
         YC(I,1)=Y(I)
         ZC(I,1)=Z(I)         
         GOTO(101,201,301,401,501,601,701,801,801,1001) IGROUP

 201     CONTINUE
C     ***  P->A
         CALL GETCCOOR(DVC,AC,2,X,Y,Z,XC,YC,ZC,I,NATR)
         GOTO 1001

 301     CONTINUE
C     ***  P->B
         CALL GETCCOOR(DVC,BC,2,X,Y,Z,XC,YC,ZC,I,NATR)
         GOTO 1001

 401     CONTINUE
C     ***  P->C
         CALL GETCCOOR(DVC,CC,2,X,Y,Z,XC,YC,ZC,I,NATR)
         GOTO 1001

 501     CONTINUE
C     ***  P->F
c         print *,'before DO'
         DO IN=2,4
            CALL GETCCOOR(DVC,FC,IN,X,Y,Z,XC,YC,ZC,I,NATR)
         ENDDO
c         print *,'after DO'
         GOTO 1001

 601     CONTINUE
C     ***  P->I
         CALL GETCCOOR(DVC,IC,2,X,Y,Z,XC,YC,ZC,I,NATR)
         GOTO 1001

 701     CONTINUE
C     ***  PR->R; primitiv rhombohedral --> rhombohedraly centered
         DO IN=2,3
c            print *,'RC>>>>',(rc(j,in),j=1,3)
            CALL GETCCOOR(DVC,RC,IN,X,Y,Z,XC,YC,ZC,I,NATR)
         ENDDO
         GOTO 1001

C     ***  P->H
 801     CONTINUE
         DO IN=2,3
            CALL GETCCOOR(DVC,HC,IN,X,Y,Z,XC,YC,ZC,I,NATR)
         ENDDO
C     *** P->CSTMC
 1001    CONTINUE
         DO IN=1,NCSTMC
            CALL GETCCOOR(DVC,CSTMC,IN,X,Y,Z,XC,YC,ZC,I,NATR)
         ENDDO
 101     CONTINUE
      ENDDO

C     WRITE COORD. OF ATOMS WITHIN CONV. CELL IN XYZ COORD.
      WRITE(12,*) 'CONVCOORD'
C     NCELL.......N. of ATOMS PER CONV. CELL
      if (igroup.lt.10) then
         NCELL=NCDET(IGROUP)
      else
         NCELL=NCSTMC  ! ncstmc was determined by GetCenteredPoints
      endif
c      print *,'DEBUG: NCELL=',ncell
      WRITE(12,*) NATR,NCELL
      DO I=1,NATR
         DO II=1,NCELL
            WRITE(12,'(i3,3x,3(f15.10,2x))')
     $           NAT(I),XC(I,II),YC(I,II),ZC(I,II)
         ENDDO
      ENDDO

c      print *,'CONVATOM_END'

      RETURN
      END


      SUBROUTINE GenEquPos
      IMPLICIT REAL*8 (a-h,o-z)
      include 'param.inc'
      INTEGER NAT(NAC)
      REAL*8 DV(3,3),SOP(3,3,48),TRX(48),TRY(48),TRZ(48),
     *     X(NAC),Y(NAC),Z(NAC),DVEC(3,3)
      REAL*8 v(3), vr(3), dvinv(3,3), a(3), d(3)
      REAL*8 rx(NAC), ry(NAC), rz(NAC)
      LOGICAL not_equal, equal
      COMMON/FTN34/ SOP,TRX,TRY,TRZ,DV,IDIM,IGROUP,NSYMOP
      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR
      
      PARAMETER (TOLMIN=1.0d-6)

      call InvertSqMat123(dv, dvinv, 3)
c      print *,'DV::',((dv(i,j),i=1,3),j=1,3)
c      print *,'IDV::',((dvinv(i,j),i=1,3),j=1,3)

c     to fractional
      do i=1,natr
         vr(1)=x(i)
         vr(2)=y(i)
         vr(3)=z(i)
         call ZeroMat(a,3,1)
         call MatMult(dvinv, 3, 3, vr, 1, a)
         rx(i)=a(1)
         ry(i)=a(2)
         rz(i)=a(3)
      enddo

      n_non_a=natr
      do i=1,n_non_a
         na = nat(i)
         v(1) = x(i)
         v(2) = y(i)
         v(3) = z(i)
         do j=1,nsymop
            vr(1) = trx(j)
            vr(2) = try(j)
            vr(3) = trz(j)
            call MatMult(sop(1,1,j),3,3,v,1,vr)
c            print *,i,j,(vr(k),k=1,3)
c     to fractional
            call ZeroMat(a,3,1)
            call MatMult(dvinv, 3, 3, vr, 1, a)
c     to [-0.5,0.5]
            do i3=1,idim
               a(i3)=a(i3) - dnint(a(i3))
            enddo
c     to Cartesian
            call ZeroMat(vr,3,1)
            call MatMult(dv, 3, 3, a, 1, vr)

c     check if we have new atom
            not_equal=.true.
            k=1
            do while (k .le. natr .and. not_equal)
               if ( nat(k) .ne. na ) goto 1
               d(1)=rx(k)-a(1)
               d(2)=ry(k)-a(2)
               d(3)=rz(k)-a(3)
               do i3=1,idim
                  d(i3)=d(i3) - dnint(d(i3))
               enddo
               if ( abs(d(1)) + abs(d(2)) + abs(d(3)) .lt. TOLMIN ) 
     $              not_equal=.false.
 1             continue
               k=k+1
            enddo
c            print *,'vr_point:',i,na,vr(1),vr(2),vr(3),.not.not_equal
            if (not_equal) then
c               print *,'DIFF:', abs(d(1)) + abs(d(2)) + abs(d(3))
               natr=natr+1
               nat(natr) = na
               x(natr) = vr(1)
               y(natr) = vr(2)
               z(natr) = vr(3)
               rx(natr)= a(1)
               ry(natr)= a(2)
               rz(natr)= a(3)
c     debug--debug
c               print *,'new_position: ', a(1), a(2), a(3)
            endif
         enddo
      enddo

      return
      END


      INTEGER FUNCTION NCDET(IGROUP)
C     DETERMINE HOW MENY ATOMS WILL BE IN CONV. CELL
      IF(IGROUP.EQ.1 .OR. IGROUP.EQ.9) NCDET=1
      IF((IGROUP.GT.1.AND.IGROUP.LT.5).OR.(IGROUP.EQ.6)) NCDET=2
      IF(IGROUP.EQ.5) NCDET=4
      IF(IGROUP.EQ.7.OR.IGROUP.EQ.8) NCDET=3
      RETURN
      END


      SUBROUTINE GETCCOOR(A33,B33,BROW,X,Y,Z,XC,YC,ZC,NA,NATR)
      include 'param.inc'
      REAL*8 B33(3,*),A33(3,3),RA(3,3),COOR(3),X(NATR),Y(NATR),Z(NATR),
     *     XC(NAC,4),YC(NAC,4),ZC(NAC,4),RX(3)
      INTEGER BROW

C     ********
C     ATTENTION: HERE WE HAVE MULT. OFF "SAME-LAYING" ELEMENTS,
C     BECAUSE B33(*,BROW) IS ROW VECTOR INSTEAD OF COLOUMN !!!!
C     ********   
c      print *,'DVC::',((a33(i,j),j=1,3),i=1,3)
c     print *,'FC::',(b33(i,brow),i=1,3)
      DO J=1,3
         COOR(J)=0.0d0
         DO I=1,3
            COOR(J)=COOR(J)+A33(I,J)*B33(I,BROW)
         ENDDO
      ENDDO
      XC(NA,BROW)=X(NA)+COOR(1)
      YC(NA,BROW)=Y(NA)+COOR(2)
      ZC(NA,BROW)=Z(NA)+COOR(3)
      
c     *****************************************************
c     COORDINATES must be within [-1,1], check that !!!
c     *****************************************************
      call InvertSqMat123( A33, RA, 3 )
      do i=1,3
         rx(i) = xc(na,brow)*ra(1,i) + yc(na,brow)*ra(2,i) +
     $        zc(na,brow)*ra(3,i)
         if ( rx(i) .gt. 1.0d0 ) then
            rx(i) = rx(i) - 1.0d0
         else if ( rx(i) .lt. -1.0d0 ) then
            rx(i) = rx(i) + 1.0d0
         endif         
      enddo
c      print *, 'DEBUG: fractC #',na,'::',rx(1),rx(2),rx(3)
      xc(na,brow) = rx(1)*a33(1,1) + rx(2)*a33(2,1) + rx(3)*a33(3,1)
      yc(na,brow) = rx(1)*a33(1,2) + rx(2)*a33(2,2) + rx(3)*a33(3,2)
      zc(na,brow) = rx(1)*a33(1,3) + rx(2)*a33(2,3) + rx(3)*a33(3,3)

      RETURN
      END


C      SUBROUTINE READARG(NATOM,F1,F2,NX,NY,NZ)
C      COMMON/DIR/ NXDIR,NYDIR,NZDIR
C      CHARACTER*80 ANX,ANY,ANZ,FILE1,FILE2
C      INTEGER C2I,F1,F2
C     this command read: NXDIR, NYDIR, NZDIR arguments
C     THIS SUBROUTINE READ THE REST OF COMMAND LINE ARGUMENTS
C     & ESTIMATE THE NUMBER OF ATOMS THAT WILL BE PRODUCED
C     WITH MULTIPLICATION (COUNT JUST ONE ATOM PER CELL)
C
C      CALL GETARG(F1,FILE1)
C      OPEN(UNIT=11,FILE=FILE1,STATUS='OLD')
C      CALL GETARG(F2,FILE2)
C      OPEN(UNIT=12,FILE=FILE2,STATUS='UNKNOWN')
C      CALL GETARG(NX,ANX)
C      CALL GETARG(NY,ANY)
C      CALL GETARG(NZ,ANZ)
C      NXDIR=C2I(ANX)
C      NYDIR=C2I(ANY)
C      NZDIR=C2I(ANZ)
c      print *,'anx,any,anz',anx,any,anz
C     NUMBER OF ATOMS/CELLS (ONE PER CELL)
C      NATOM=(NXDIR+2)*(NYDIR+2)*(NZDIR+1)
c      print *,'nx,ny,nz>',nxdir,nydir,nzdir
C      RETURN
C      END



C     THIS SUBROUTINE IS USED WHEN WE ARE DEALING WITH MOLECULE:
C     IT WRITE A COORDINATES OF MOLECULE/CLUSTER
      SUBROUTINE MOLECULE
      include 'param.inc'
      INTEGER NAT(NAC)
      REAL*8 X(NAC),Y(NAC),Z(NAC),DVEC(3,3)
      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR
      WRITE(12,*) 'ATOMS'
      DO I=1,NATR
         WRITE(12,21) NAT(I),X(I),Y(I),Z(I)
      ENDDO
 21   FORMAT(I5,3F16.10)
      RETURN
      END
         
      LOGICAL Function IsReverse(CSTMC,RCR)
      REAL*8 CSTMC(3,4), RC(3,3), RCR(3,3), TOLMIN
      LOGICAL isequal
      TOLMIN=1.0d-6
c     R(x/y/z,#point)
      do ir=1,3
         isequal=.FALSE.
c     can we find RCR(*,ir) point among CSTMC(*,#) points
         do ic=1,3
            do j=1,3
               if ( abs(CSTMC(j,ic)-RCR(j,ir)) .lt. TOLMIN )
     $              isequal=.TRUE.
            enddo
         enddo
         if (.not.isequal) then
c     RCR(*,ir) point wasn't find -> OBVERSE SETTINGS
            IsReverse=.FALSE.
            return
         endif
      enddo
      IsReverse=.TRUE.
      return
      End
