      SUBROUTINE TTCBOXES(TTC)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 TTCVEC(3,6,8),TTC(3,13),TTC12(3,2,8),TTCPOINT(3,7,8)
      INTEGER TTCBOX(5,8),NTTCVEC(8),TTCPAIR(2,6,8)
      REAL*8 MINTOL,NULL

      COMMON/TTCBOXS/ TTCVEC,NTTCVEC,TTCPAIR,TTCPOINT

      PARAMETER (CON13=1.0d0/3.0d0, CON23=2.0d0/3.0d0, MINTOL=1.d-6, 
     1     NULL=0.0d0, TWO=2.0d0, TWON=-2.0d0)

C     ************************************************************************
C     TTCVEC --- vectors for all 8 boxes
C     TTC    --- coordinates for RTTC/HTTC points
C     TTC12  --- to TTC points we add one point in front and one on the back
C     TTCBOX --- indexes of TTC points, which are in box
C     NTTCVEC--- number of vectors in box
C     TTCPAIR--- which pair of points are basis for vector
C     ************************************************************************
C     ***** to each "ttcbox box" we add first and last point, 
C     here are coordinates *****
      DATA TTC12/
     1     NULL, TWO, NULL,
     1     TWO, NULL, NULL,
     2     TWON,TWON, NULL,
     2     TWO, NULL, NULL,
     3     TWON,TWON, NULL,
     3     TWO, TWO, NULL,
     4     NULL,TWON, NULL,
     4     TWO, TWO, NULL,
     5     NULL,TWON, NULL,
     5     TWON, NULL, NULL,
     6     TWO, TWO, NULL,
     6     TWON, NULL, NULL,
     7     TWO, TWO, NULL,
     7     TWON,TWON, NULL,
     8     NULL, TWO, NULL,
     8     TWON,TWON, NULL/

C     *** indexes of "BOXs" POINTS ***
      DATA TTCBOX/
     1     3,4,5,6,1,
     2     5,6,1,0,0,
     3     5,6,1,2,0,
     4     6,1,2,0,0,
     5     6,1,2,3,4,
     6     2,3,4,0,0,
     7     2,3,4,5,0,
     8     3,4,5,0,0/
C     *** how many points are in each box ***
      DATA NTTCVEC/6,4,5,4,6,4,5,4/

C     *** from TTC & TTC12 make TTCPOINTs ***
C     *** insert first point ***
      DO I=1,8
         DO K=1,3
            TTCPOINT(K,1,I)=TTC12(K,1,I)
         ENDDO
      ENDDO
C     *** TTC -> TTCPOINT ***
      DO I=1,8   !four boxes 
         DO J=2,NTTCVEC(I)       !vectors per boxes
            DO K=1,3  !xyz coor of point        
               TTCPOINT(K,J,I)=TTC(K,TTCBOX(J-1,I))
            ENDDO
         ENDDO
      ENDDO
C     *** insert last point ***
      DO I=1,8
         DO K=1,3
            TTCPOINT(K,NTTCVEC(I)+1,I)=TTC12(K,2,I)
         ENDDO
      ENDDO

C     FROM POINTS -> VECTORS
c      print *,'TTCVEC>'
      DO I=1,8   !four boxes 
         DO J=1,NTTCVEC(I)       !vectors per boxes
            DO K=1,3  !xyz coor of point        
               TTCVEC(K,J,I)=TTCPOINT(K,J+1,I)-TTCPOINT(K,J,I)
            ENDDO
            write(6,*) (TTCVEC(k,j,i),k=1,3)
         ENDDO
c         print *,'------------------------------'
      ENDDO

c      print *,'TTCPOINT'
c      do i=1,8
c         do j=1,7
c            print *,(ttcpoint(k,j,i),k=1,3)
c         enddo
c         print *,'-------------------------------'         
c      enddo

C     **** VECPAIR --- prvi vector je vector, ki skupaj z C tvori ravnino,
C     ****             drugi pa pove katera stran je "prava"
      DATA TTCPAIR/
C     *** "BOX1" PAIRS ***
     1     1,-2,
     2     2,3,
     3     3,4,
     4     4,5,
     5     5,-6,
     6     6,5,     
C     *** "BOX2" PAIRS ***     
     1     1,-2,
     2     2,3,
     3     3,-4,
     4     4,3,
     5     0,0,
     6     0,0,     
C     *** "BOX3" PAIRS ***
     1     1,-2,
     2     2,3,
     3     3,4,
     4     4,-5,
     5     5,4,
     6     0,0,     
C     *** "BOX4" PAIRS ***
     1     1,-2,
     2     2,3,
     3     3,-4,
     4     4,3,
     5     0,0,
     6     0,0,     
C     *** "BOX5" PAIRS
     1     1,-2,
     2     2,3,
     3     3,4,
     4     4,5,
     5     5,-6,
     6     6,5,     
C     *** "BOX6" PAIRS
     1     1,-2,
     2     2,3,
     3     3,-4,
     4     4,3,
     5     0,0,
     6     0,0,
C     *** "BOX7" PAIRS
     1     1,-2,
     2     2,3,
     3     3,4,
     4     4,-5,
     5     5,4,
     6     0,0,
C     *** "BOX8" PAIRS
     1     1,-2,
     2     2,3,
     3     3,-4,
     4     4,3,
     5     0,0,
     6     0,0/     

      RETURN
      END

      
C     tocka1(p1,p2,p3) pove kje je skatla, tocka2(r1,r2,r3) pa katero tocko
C     ocenjujemo
      LOGICAL FUNCTION TTCBOX(NBOX,P1,P2,P3,R1,R2,R3,DVEC)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 TTCPOINT(3,7,8),TTCVEC(3,6,8),DVEC(3,3),IDVEC(3,3)
      INTEGER NTTCVEC(8),TTCPAIR(2,6,8),TMPVECP1,TMPVECP2
      LOGICAL ISINTTC,LTTC(2,8),LDONE

      COMMON/TTCBOXS/ TTCVEC,NTTCVEC,TTCPAIR,TTCPOINT
      COMMON/IVEC/ IDVEC  !inverse matrix of DVEC(3,3)

      PARAMETER(CON13=1.0d0/3.0d0, CON23=2.0d0/3.0d0)

      TTCBOX=.TRUE.
      LDONE=.FALSE.

C     x,y,z are "relative coord." of point in box
      x=r1-p1
      y=r2-p2
      z=r3-p3
C     transform (x,y,z) to DVEC basis
      xx=x*IDVEC(1,1)+y*IDVEC(2,1)+z*IDVEC(3,1)
      yy=x*IDVEC(1,2)+y*IDVEC(2,2)+z*IDVEC(3,2)
      zz=x*IDVEC(1,3)+y*IDVEC(2,3)+z*IDVEC(3,3)
c      print *,'xx,yy,zz>',xx,yy,zz

      LTTC(1,1)=YY.GT.1.0d0
      LTTC(2,1)=XX.GT.1.0d0
      LTTC(1,2)=XX.LT.-1.0d0.AND.YY.LT.-1.0d0
      LTTC(2,2)=XX.GT.1.0d0
      LTTC(1,3)=XX.LT.-1.0d0.AND.YY.LT.-1.0d0
      LTTC(2,3)=XX.GT.1.0d0.AND.YY.GT.1.0d0
      LTTC(1,4)=YY.LT.-1.0d0
      LTTC(2,4)=XX.GT.1.0d0.AND.YY.GT.1.0d0
      LTTC(1,5)=YY.LT.-1.0d0
      LTTC(2,5)=XX.LT.-1.0d0
      LTTC(1,6)=XX.GT.1.0d0.AND.YY.GT.1.0d0
      LTTC(2,6)=XX.LT.-1.0d0
      LTTC(1,7)=XX.GT.1.0d0.AND.YY.GT.1.0d0
      LTTC(2,7)=XX.LT.-1.0d0.AND.YY.LT.-1.0d0
      LTTC(1,8)=YY.GT.1.0d0
      LTTC(2,8)=XX.LT.-1.0d0.AND.YY.LT.-1.0d0

C     *** first check 3th part of the box **
      I=NTTCVEC(NBOX)
      PX=P1+TTCPOINT(1,I,NBOX)*DVEC(1,1)+
     1     TTCPOINT(2,I,NBOX)*DVEC(2,1)+TTCPOINT(3,I,NBOX)*DVEC(3,1)
      PY=P2+TTCPOINT(1,I,NBOX)*DVEC(1,2)+
     1     TTCPOINT(2,I,NBOX)*DVEC(2,2)+TTCPOINT(3,I,NBOX)*DVEC(3,2)
      PZ=P3+TTCPOINT(1,I,NBOX)*DVEC(1,3)+
     1     TTCPOINT(2,I,NBOX)*DVEC(2,3)+TTCPOINT(3,I,NBOX)*DVEC(3,3)
C     is (R1,R2,R3) in the 3TH part of the box
      if(lttc(2,nbox))then
c         print *,'PART> 3th',xx,yy
         ldone=.true.
         if(.not.isinttc(nbox,ttcpair(1,i,nbox),ttcpair(2,i,nbox),
     1        px,py,pz,r1,r2,r3,dvec)) ttcbox=.false.
      endif   

      DO I=1,NTTCVEC(NBOX)-1
c         print *,'NTTCVEC> i=',i,'; xx,yy,zz>',xx,yy,zz
         TMPVECP1=TTCPAIR(1,I,NBOX)
         TMPVECP2=TTCPAIR(2,I,NBOX)  
c         print *,'vec1=',tmpvecp1
c         print *,'vec2=',tmpvecp2         
C     IS POINT(R1,R2,R3) ON THE RIGHT SIDE OF THE PLANE?????
         PX=P1+TTCPOINT(1,I,NBOX)*DVEC(1,1)+
     1        TTCPOINT(2,I,NBOX)*DVEC(2,1)+TTCPOINT(3,I,NBOX)*DVEC(3,1)
         PY=P2+TTCPOINT(1,I,NBOX)*DVEC(1,2)+
     1        TTCPOINT(2,I,NBOX)*DVEC(2,2)+TTCPOINT(3,I,NBOX)*DVEC(3,2)
         PZ=P3+TTCPOINT(1,I,NBOX)*DVEC(1,3)+
     1        TTCPOINT(2,I,NBOX)*DVEC(2,3)+TTCPOINT(3,I,NBOX)*DVEC(3,3)
c         print *,'PXYZ> ',px,py,pz
C     is (R1,R2,R3) in the 1ST part of the box; 
C     for 1ST part i must be 1 -> this means: first point and first vector
         if(i.eq.1.and.lttc(1,nbox))then
c            print *,'PART> 1st',xx,yy
            ldone=.true.
            if(.not.isinttc(nbox,tmpvecp1,tmpvecp2,
     1           px,py,pz,r1,r2,r3,dvec)) ttcbox=.false.
         endif
C     is (R1,R2,R3) in the 2ND part of the box:
         if(i.gt.1.and.(.not.ldone))then
c            print *,'PART> 2nd',xx,yy
            if(.not.isinttc(nbox,tmpvecp1,tmpvecp2,
     1           px,py,pz,r1,r2,r3,dvec)) ttcbox=.false.
         endif
      ENDDO                     !DO I=1,NVEC(NBOX)

      RETURN
      END

 
      LOGICAL FUNCTION ISINTTC(NBOX,NVC1,NVC2,P1,P2,P3,R1,R2,R3,dvec)
      REAL*8 P1,P2,P3,R1,R2,R3
      REAL*8 DVEC(3,3),MINTOL
      REAL*8 A,B,C    !components of normal vector
      REAL*8 D,KFC   ! Ax + By + Cz + D = 0
      REAL*8 AVEC,BVEC,CVEC,AVEC2,BVEC2,CVEC2
C     if nvec2<0 negative=.true.; third vector is turned upside down
      LOGICAL NEGATIVE  
      REAL*8 TTCPOINT(3,7,8),TTCVEC(3,6,8)
      INTEGER NTTCVEC(8),TTCPAIR(2,6,8)
      
      COMMON/TTCBOXS/ TTCVEC,NTTCVEC,TTCPAIR,TTCPOINT
      
      PARAMETER (MINTOL=1.0d-6)
     
      NVEC1=NVC1
      NVEC2=NVC2
c      print *,'nvec1=',nvec1,', nvec2=',nvec2
      NEGATIVE=.FALSE.
      ISINTTC=.FALSE.
      IF(NVEC2.LT.0.)THEN
         NVEC2=-NVEC2
         NEGATIVE=.TRUE.
      ENDIF

C     plvec is in coloumn major mode, but dvec is in row-major mode
C     *** transform plvec1 in XYZ basis ***
      AVEC=DVEC(1,1)*TTCVEC(1,NVEC1,NBOX)+
     1     DVEC(2,1)*TTCVEC(2,NVEC1,NBOX)+DVEC(3,1)*TTCVEC(3,NVEC1,NBOX)
      BVEC=DVEC(1,2)*TTCVEC(1,NVEC1,NBOX)+
     1     DVEC(2,2)*TTCVEC(2,NVEC1,NBOX)+DVEC(3,2)*TTCVEC(3,NVEC1,NBOX)
      CVEC=DVEC(1,3)*TTCVEC(1,NVEC1,NBOX)+
     1     DVEC(2,3)*TTCVEC(2,NVEC1,NBOX)+DVEC(3,3)*TTCVEC(3,NVEC1,NBOX)
      A=BVEC*DVEC(3,3)-DVEC(3,2)*CVEC
      B=CVEC*DVEC(3,1)-DVEC(3,3)*AVEC
      C=AVEC*DVEC(3,2)-DVEC(3,1)*BVEC

C     *** transform ttcvec2 in XYZ basis ***
      AVEC2=DVEC(1,1)*TTCVEC(1,NVEC2,NBOX)+
     1     DVEC(2,1)*TTCVEC(2,NVEC2,NBOX)+DVEC(3,1)*TTCVEC(3,NVEC2,NBOX)
      BVEC2=DVEC(1,2)*TTCVEC(1,NVEC2,NBOX)+
     1     DVEC(2,2)*TTCVEC(2,NVEC2,NBOX)+DVEC(3,2)*TTCVEC(3,NVEC2,NBOX)
      CVEC2=DVEC(1,3)*TTCVEC(1,NVEC2,NBOX)+
     1     DVEC(2,3)*TTCVEC(2,NVEC2,NBOX)+DVEC(3,3)*TTCVEC(3,NVEC2,NBOX)
      IF(A*AVEC2+B*BVEC2+C*CVEC2.LT.0.0d0)THEN
         A=-A
         B=-B
         C=-C
      ENDIF

c      print *,'TTCVEC:',TTCVEC(1,NVEC1,NBOX),TTCVEC(2,NVEC1,NBOX),
c     1     TTCVEC(3,NVEC1,NBOX)*DVEC(3,1)
c      print *,'DVEC',DVEC(3,1),DVEC(3,2),DVEC(3,3)
c      print *,'p1,p2,p3=',p1,p2,p3
c      print *,'r1,r2,r3=',r1,r2,r3
      D=-A*P1-B*P2-C*P3
c      print *,'a,b,c,d::',a,b,c,d
      KFC=A*R1+B*R2+C*R3+D
c      print *,'kfc1=',kfc,', kfc2=',A*R1+B*R2+C*R3-(-A*P1-B*P2-C*P3)
      IF(KFC.GE.(0.0d0-MINTOL).AND.(.NOT.NEGATIVE)) ISINTTC=.TRUE.
      IF(KFC.LE.(0.0d0+MINTOL).AND.NEGATIVE) ISINTTC=.TRUE.
c      print *,'ISINTTC=',isinttc

      RETURN
      END


      
