      PROGRAM ATOMLAB
C     PROGRAM FOR FINDING LABEL OF ATOM TO BE SUBSITUTED/REMOVED
C     ATOMSUBS, ATOMREMO OPTIONS

C     USAGE: atomlab <type> <file1> <file2> <cr98>
C     
C     type == 1  return atomlabel
C     type == 2  take coordinates of added atom and translate it to cell 0
C
C     FILE1: gengeom file (UNIT 11)
C     FILE2: file where NAT X Y Z are writen (STDIN)
C
C     cr98: fourth argument is present if CRYSTAL version is 98 or newer; 
C           read coorinates rather from unit33
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      INTEGER NAT(NAC)
      REAL*8 X(NAC),Y(NAC),Z(NAC),DVEC(3,3),IDVEC(3,3)
      REAL*8 A(NAC),B(NAC),C(NAC),MINDIS2,MINTOL
      CHARACTER*80 FILE1,FILE2
      CHARACTER ARG1
      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR
      COMMON/IVEC/ IDVEC
      COMMON/BOHRR/ BOHR,IMODE3 !this is not used, but is needed because 
                                !this program uses readf1 routine

      PARAMETER (MINTOL=1.D-4)

c     this is just because of common/bohrr/
      IMODE3=0
      BOHR=0

C     TYPE OF RUN
C     1 ... return atom label
C     2 ... return atom coordinates after translation to basic cell
      CALL GETARG(1,ARG1)

c     2nd argument

      CALL GETARG(2,FILE1)
      OPEN(UNIT=11,FILE=FILE1,STATUS='OLD')
C     read GENGEOM FILE
      CALL READF1('DIM-GROUP',IGROUP,IDIM)
      CALL READF1('PRIMVEC',0,IDIM)
      CALL READF1('RECIP-PRIMVEC',0,IDIM)
      CALL READF1('PRIMCOORD',0,IDIM)
      CLOSE(11)

c     3rd argument

      CALL GETARG(3,FILE2)
      OPEN(UNIT=13,FILE=FILE2,STATUS='OLD')
C     read FILE2
      READ(13,*) ATN,XS,YS,ZS
      CLOSE(13)

c     4th argument

C     does FOURTH ARGUMENT exists
      if (Iargc() .eq. 4) then
C     read XYZ coordinates form unit 33 (XYZ format file)
         Call ReadXYZ(NATR, NAT, X, Y, Z)
      endif

C     from XYZ coordinates to ABC coordinates of primitiv cell
      DO I=1,NATR
         IF(IDIM .eq. 1) A(I)=X(I)*IDVEC(1,1)
         IF(IDIM .eq. 2) THEN
            A(I)=X(I)*IDVEC(1,1)+Y(I)*IDVEC(2,1)
            B(I)=X(I)*IDVEC(1,2)+Y(I)*IDVEC(2,2)
         ENDIF
         IF(IDIM .eq. 3) THEN
            A(I)=X(I)*IDVEC(1,1)+Y(I)*IDVEC(2,1)+Z(I)*IDVEC(3,1)
            B(I)=X(I)*IDVEC(1,2)+Y(I)*IDVEC(2,2)+Z(I)*IDVEC(3,2)
            C(I)=X(I)*IDVEC(1,3)+Y(I)*IDVEC(2,3)+Z(I)*IDVEC(3,3)
         ENDIF
      ENDDO

      IF(IDIM .eq. 1) THEN 
         AS=XS*IDVEC(1,1)
      ELSEIF(IDIM .eq. 2) THEN
         AS=XS*IDVEC(1,1)+YS*IDVEC(2,1)
         BS=XS*IDVEC(1,2)+YS*IDVEC(2,2)      
      ELSEIF(IDIM .eq. 3) THEN
         AS=XS*IDVEC(1,1)+YS*IDVEC(2,1)+ZS*IDVEC(3,1)
         BS=XS*IDVEC(1,2)+YS*IDVEC(2,2)+ZS*IDVEC(3,2)
         CS=XS*IDVEC(1,3)+YS*IDVEC(2,3)+ZS*IDVEC(3,3)
      ENDIF
C     now we must relate (AS,BS,CS) to one of (A(i),B(i),C(i))
C     translate (AS,BS,CS) in the first cell
      IF(IDIM .gt. 0) AS=AS-DBLE(NINT(AS))      !AS is between [-0.5,0.5]
      IF(IDIM .gt. 1) BS=BS-DBLE(NINT(BS))
      IF(IDIM .gt. 2) CS=CS-DBLE(NINT(CS))
      
      IF(IDIM .eq. 1) XS=AS*DVEC(1,1)
      IF(IDIM .eq. 2) THEN
         XS=AS*DVEC(1,1)+BS*DVEC(2,1)
         YS=AS*DVEC(1,2)+BS*DVEC(2,2)
      ENDIF
      IF(IDIM .eq. 3) THEN
         XS=AS*DVEC(1,1)+BS*DVEC(2,1)+CS*DVEC(3,1)
         YS=AS*DVEC(1,2)+BS*DVEC(2,2)+CS*DVEC(3,2)
         ZS=AS*DVEC(1,3)+BS*DVEC(2,3)+CS*DVEC(3,3)
      ENDIF
      
C     IF TYPE OF RUN IS 2 --> RETURN (AS,BS,CS)
      IF(ARG1.EQ.'2')THEN         
         WRITE(6,'(3(F18.10,2x))') XS,YS,ZS
         GOTO 1000
      ENDIF

c      print *,'as,bs,cs',as,bs,cs
      MINDIS2=3.0d0 !this is maximum DSQRT(distance) within one cell
      DO I=1,NATR
C     in CRYSTAL95 all coordinates of non-equivalent atoms in primitiv
C     cell are within [-0.5,0.5], but because of roundoff error something
C     that is actually a little lower than 0.5 became a litle greter than 0.5
C     and instead of ~0.5 we get ~-0.5 -> CHECK THAT
C     THIS WILL BE CHECKED BY NEAREST DISTANCE (THIS IS SIMILAR TO MC 
C     BOUNDARY CONDITIONS)
         
         IF (IDIM .eq. 0) THEN
            ADIS2=(X(i)-XS) * (X(i)-XS)
            BDIS2=(Y(i)-YS) * (Y(i)-YS)
            CDIS2=(Z(i)-ZS) * (Z(i)-ZS)
         ENDIF
         IF (IDIM .eq. 1) THEN
            ADIS2=DISTANCE2(A(i),AS)
            BDIS2=(Y(i)-YS) * (Y(i)-YS)
            CDIS2=(Z(i)-ZS) * (Z(i)-ZS)
         ENDIF
         IF (IDIM .eq. 2) THEN
            ADIS2=DISTANCE2(A(i),AS)
            BDIS2=DISTANCE2(B(i),BS)
            CDIS2=(Z(i)-ZS) * (Z(i)-ZS)
         ENDIF
         IF (IDIM .eq. 3) THEN
            ADIS2=DISTANCE2(A(i),AS)
            BDIS2=DISTANCE2(B(i),BS)
            CDIS2=DISTANCE2(C(i),CS)
         ENDIF
         DIS2=ADIS2+BDIS2+CDIS2
c         print *,'dis2=',dis2
         IF(DIS2.LT.MINDIS2)THEN
            MINDIS2=DIS2
            LABEL=I
         ENDIF
      ENDDO
c      print *,mindis2,mintol
C     DSQRT(MINDIS2) MUST BE LOWER THAN MINTOL, OTHERWISE WE DIDN'T FIND
C     APPROPRIATE ATOM
      NAT100 = NAT(LABEL) - INT(NAT(LABEL)/100)*100 ! ATN is within [0..99], so should be NAT100
      IF(DSQRT(MINDIS2).LT.MINTOL.AND.NAT100.EQ.ATN)THEN
         WRITE(6,*) LABEL
      ELSE
         STOP 'WARNING: no appropriate atom was found !!!!'
      ENDIF
 1000 CONTINUE
      END


c     get the square of distance
      REAL*8 FUNCTION DISTANCE2(A1,A2)
      REAL*8 A1,A2
      IF(DABS(A1-A2).GT.0.5d0)THEN
         DISTANCE2=(1.0d0-DABS(A1-A2))*(1.0d0-DABS(A1-A2))
      ELSE
         DISTANCE2=(A1-A2)*(A1-A2)
      ENDIF
      RETURN
      END

