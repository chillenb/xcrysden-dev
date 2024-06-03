c     naredil bom takole: nml=AxB v (x,y,z) reprezentaciji, nato pa
c     gre nml v (h*,k,*l*) reprezentacijo
c
c     |x|     |a|     |a|       |x|
c     |y| = A |b| ==> |b| = A^-1|y|  
c     |z|     |c|     |c|       |z|
c     
c     A = matrix of reciprocal vectors ==> A^-1 is a matrix of direct vectors

c     ================
      PROGRAM CALPLANE 
      CHARACTER*80 file1,file2
      REAL dv1(3),dv2(3),dv3(3),p1(3),p2(3),p3(3)
      REAL vec1(3),vec2(3),nml(3)
      CHARACTER*7 word
      PARAMETER (MINTOL=0.001)
      INTEGER ia(3)
      COMMON/INDX/ IA
C     This program calculates selected plane 
c     file1...direct vectors coord.
c     file2...coordinates of 3 selected points in XYZ

      CALL GETARG(1,FILE1)
      CALL GETARG(2,FILE2)
      OPEN(11,FILE=FILE1,STATUS='OLD')
      OPEN(12,FILE=FILE2,STATUS='OLD')

 20   CONTINUE
      READ(11,'(1x,a)') WORD
      IF(WORD.NE.'CONVVEC')GOTO20 
c     read direct vectors
      READ(11,*) (DV1(i),i=1,3)
      READ(11,*) (DV2(i),i=1,3)
      READ(11,*) (DV3(i),i=1,3)

c     translate to reciprocal
C     read 3 sel. points
      READ(12,*) (P1(i),i=1,3)
      READ(12,*) (P2(i),i=1,3)
      READ(12,*) (P3(i),i=1,3)
      CLOSE(11)
      CLOSE(12)

c     vec1 is (P1)-(P2); vec2 is (P3)-(P2)
      do i=1,3
         vec1(i)=P1(i)-P2(i)      
         vec2(i)=P3(i)-P2(i)
      enddo
      
      CALL FVecProduct(vec1,vec2,nml)
c     translate to [h*k*l*] basis
      AH=nml(1)*DV1(1) + nml(2)*DV2(1) + nml(3)*DV3(1)
      AK=nml(1)*DV1(2) + nml(2)*DV2(2) + nml(3)*DV3(2)
      AL=nml(1)*DV1(3) + nml(2)*DV2(3) + nml(3)*DV3(3)
      
C     ROUND the components of normal vector --> Just good enough
C     for unit cells without BASIS 

      CALL GETINDX(AH,AK,AL)
C     choose Miller indexes so that there will be a minority of negative
C     numbers
      NNEG=0
      NNUL=0
      DO I=1,3
         IF(IA(I).LT.0)NNEG=NNEG+1
         IF(IA(I).EQ.0)NNUL=NNUL+1
      ENDDO
      IF(((NNUL.EQ.2).AND.(NNEG.EQ.1)).OR.
     *     (NNEG.GE.2))THEN
         DO I=1,3
            IA(I)=-IA(I)
         ENDDO
      ENDIF
C     first unrounded M. indexes, and then rounded indexes
      WRITE(6,*) AH,AK,AL,(IA(I),I=1,3)

      END


c     =======================================
      SUBROUTINE FVecProduct(vec1,vec2,resvec)
      IMPLICIT REAL (A-H,O-Z)
      REAL vec1(3),vec2(3),resvec(3)
      
      resvec(1) = vec1(2)*vec2(3) - vec2(2)*vec1(3)
      resvec(2) = vec1(3)*vec2(1) - vec2(3)*vec1(1)
      resvec(3) = vec1(1)*vec2(2) - vec2(1)*vec1(2)
      RETURN
      END


c     ============================
      SUBROUTINE GETINDX(AH,AK,AL)
      INTEGER IFC(2),IA(3)
      REAL A(3),B(2),MIN
      COMMON/INDX/ IA

      A(1)=REAL(AH)
      A(2)=REAL(AK)
      A(3)=REAL(AL)

c     tolerance
      TOL1=0.001
      TOL2=0.001

c     determine the lowest number
      MIN=ABS(A(1))+ABS(A(2))+ABS(A(3))
      DO I=1,3
         IF(ABS(A(I)).LT.MIN.AND.ABS(A(I)).GT.0.0)THEN
            MIN=A(I)
            II=I
         ENDIF
      ENDDO
c     devide all numbers by the smallest one
c     two larger numbers are copied to b()
      N=0
      DO I=1,3
         A(I)=A(I)/MIN
         IF(I.NE.II)THEN
            N=N+1
            B(N)=A(I)
         ENDIF
      ENDDO
      
      ifc(1)=1
      ifc(2)=1
c     find an integer multiplier for b()
      DO 92 J=1,2
         I=0
 91      CONTINUE
         I=I+1
         B1=ABS(REAL(I)*B(J))
         IF(ABS(B1).LT.1.E-5)GOTO92 !B(J) is zero
         B2=ABS(B1-ABS(REAL(NINT(REAL(I)*B(J)))))
         IF(B2.LT.TOL2)THEN
            IFC(J)=I
         ELSE
c     some combination of triplet-numbers can produce realy huge numbers,
c     so this loop can never and, lets say that the greater allowed I is 10000
            if(i.gt.10000)then
               tol2=2.*tol2
               i=0
            endif
            GOTO91
         ENDIF
 92   CONTINUE
      
c     what is the smallest common multiplier for b()
      IDEL=1
      IMIN=ABS(IFC(1))
      IF(ABS(IFC(2)).LT.IMIN)IMIN=ABS(IFC(2))
      DO I=2,IMIN
         AD=REAL(IFC(1))/REAL(I)-REAL(INT(IFC(1)/I))
         BD=REAL(IFC(2))/REAL(I)-REAL(INT(IFC(2)/I))
         IF(AD.LT.TOL1.AND.BD.LT.TOL1)IDEL=IDEL*I
      ENDDO
c     update result
      IFAC=IFC(1)*IFC(2)/IDEL      
      DO J=1,3
         IA(J)=NINT(A(J)*IFAC) !Miller indexes
      ENDDO

      RETURN
      END
