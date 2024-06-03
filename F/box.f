C     ************************************************************************
C     ** rutina pove, ce je tocka(R1,R2,R3) na isti strani ravnine kot DVEC[k]
C     ** *** ravnino dolocata vector DVEC[i],DVEC[j] in tocka(P1,P2,P3)  *****
      LOGICAL FUNCTION ISINSIDE(DVEC,I,J,K,P1,P2,P3,R1,R2,R3)
      REAL*8 DVEC(3,3),P1,P2,P3,R1,R2,R3,KFC,MINTOL
      REAL*8 A,B,C    !components of normal vector
      REAL*8 D   ! Ax + By + Cz + D = 0
C     if k<0 negative=.true.; third vector is turned upside down
      LOGICAL NEGATIVE  
      PARAMETER (MINTOL=1.0d-3)
 
c      print *,p1,p2,p3,r1,r2,r3,i,j,k
c      print *,'k=',k
      NEGATIVE=.FALSE.
      ISINSIDE=.FALSE.
      kk=k
      IF(KK.LT.0)THEN
         KK=-KK
         NEGATIVE=.true.
      ENDIF

      A=DVEC(I,2)*DVEC(J,3)-DVEC(J,2)*DVEC(I,3)
      B=DVEC(I,3)*DVEC(J,1)-DVEC(J,3)*DVEC(I,1)
      C=DVEC(I,1)*DVEC(J,2)-DVEC(J,1)*DVEC(I,2)
      
      IF(A*DVEC(KK,1)+B*DVEC(KK,2)+C*DVEC(KK,3).LT.0.0d0)THEN
         A=-A
         B=-B
         C=-C
      ENDIF

      D=-A*P1-B*P2-C*P3
      KFC=A*R1+B*R2+C*R3+D

      IF(KFC.GE.(0.0d0-MINTOL).AND.(.NOT.NEGATIVE)) ISINSIDE=.TRUE.
      IF(KFC.LE.(0.0d0+MINTOL).AND.NEGATIVE) ISINSIDE=.TRUE.
c      print *,'INSIDE=',isinside

      RETURN
      END


      SUBROUTINE BOXES
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PLPOINT(3,4,6),PLVEC(3,3,6)
      INTEGER NPOINT(6),NVEC(6),VECPAIR(2,3,6)
      REAL*8 MINTOL, NULL
      COMMON/BOXS/ PLVEC,NVEC,VECPAIR,PLPOINT

      PARAMETER (CON13=1.0d0/3.0d0, CON23=2.0d0/3.0d0, MINTOL=1.d-6, 
     1     NULL=0.0d0)
      
C     NPOINT(4) TELLS NUMBER OF POINTS FOR EACH "BOX"
      DATA NPOINT/4,4,3,4,4,3/
C     NVEC(4) tells number of vectors for each "box"
      DATA NVEC/3,3,2,3,3,2/

      DATA PLPOINT/
C     *** "BOX1" POINTS ***
     1     1.0d0, 0.5d0, NULL,
     2     CON23, CON13, NULL,
     3     CON13, CON23, NULL,
     4     0.5d0, 1.0d0, NULL,
C     *** "BOX2" POINTS ***
     1     1.0d0, 0.5d0, NULL,
     2     CON23, CON13, NULL,
     3     CON13, CON23, NULL,
     4     NULL, 0.5d0, NULL,
C     *** "BOX3" POINTS ***
     1     0.5d0, 1.0d0, NULL,
     2     CON13, CON23, NULL,
     3     NULL, 0.5d0, NULL,
     4     NULL, NULL, NULL,   !"BOX3" has just three points
C     *** "BOX4" POINTS ***
     1     0.5d0, NULL, NULL,
     2     CON23, CON13, NULL,
     3     CON13, CON23, NULL,
     4     0.5d0, 1.0d0, NULL,
C     *** "BOX5" POINTS ***
     1     0.5d0, NULL, NULL,
     2     CON23, CON13, NULL,
     3     CON13, CON23, NULL,
     4     NULL, 0.5d0, NULL,
C     *** "BOX6" POINTS ***
     1     1.0d0, 0.5d0, NULL,
     2     CON23, CON13, NULL,
     3     0.5d0, NULL, NULL,
     4     NULL, NULL, NULL/ !"BOX6" has just three points

c      do i=1,6
c         do j=1,4
c            print *,(plpoint(k,j,i),k=1,3)
c         enddo
c         print *,'---------------------'
c      enddo

C     FROM POINTS -> VECTORS
      DO I=1,6   !four boxes 
         DO J=1,NVEC(I)       !vectors per boxes
            DO K=1,3  !xyz coor of point        
               PLVEC(K,J,I)=PLPOINT(K,J+1,I)-PLPOINT(K,J,I)
            ENDDO
            write(6,*) (plvec(k,j,i),k=1,3)
         ENDDO
c         print *,'------------------------------'
      ENDDO
           
C     **** VECPAIR --- prvi vector je vector, ki skupaj z C tvori ravnino,
C     ****             drugi pa pove katera stran je "prava"
      DATA VECPAIR/
C     *** "BOX1" PAIRS ***
     1     1,2,
     2     2,3,
     3     3,-1,
C     *** "BOX2" PAIRS ***     
     1     1,2,
     2     2,-3,
     3     3,2,
C     *** "BOX3" PAIRS ***
     1     1,2,
     2     2,-1,
     3     0,0,   !BOX3 HAS JUST TWO PAIRS     
C     *** "BOX4" PAIRS ***
     1     1,2,
     2     2,-3,
     3     3,2,
C     *** "BOX5" PAIRS
     1     1,3,
     2     2,3,
     3     3,-2,
C     *** "BOX6" PAIRS
     1     1,2,
     2     2,-1,
     3     0,0/   !BOX6 HAS JUST TWO PAIRS   

      RETURN
      END

      
C     tocka1(p1,p2,p3) pove kje je skatla, tocka2(r1,r2,r3) pa katero tocko
C     ocenjujemo
      LOGICAL FUNCTION BOX(NBX,P1,P2,P3,R1,R2,R3,DVEC)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PLPOINT(3,4,6),PLVEC(3,3,6),DVEC(3,3),IDVEC(3,3)
      INTEGER NVEC(6),VECPAIR(2,3,6),TMPVECP1,TMPVECP2
      LOGICAL ISINBOX

      COMMON/BOXS/ PLVEC,NVEC,VECPAIR,PLPOINT
      COMMON/IVEC/ IDVEC  !inverse matrix of DVEC(3,3)

      PARAMETER(CON13=1.0d0/3.0d0, CON23=2.0d0/3.0d0)

c      print *,'IN BOX !!!'
      BOX=.TRUE.
      ISIDE=1
      NBOX=NBX
C     if nbox < 0 we must turn around the box
      IF(NBOX.LT.0)THEN
         NBOX=-NBOX
         ISIDE=-1
      ENDIF
         
C     if we deal with NBOX=2.or.NBOX=4, we will need this
      x=r1-p1
      y=r2-p2
      z=r3-p3
C     transform (x,y,z) to DVEC basis
      xx=x*IDVEC(1,1)+y*IDVEC(2,1)+z*IDVEC(3,1)
      yy=x*IDVEC(1,2)+y*IDVEC(2,2)+z*IDVEC(3,2)
      zz=x*IDVEC(1,3)+y*IDVEC(2,3)+z*IDVEC(3,3)
c      print *,'xx,yy,zz>',xx,yy,zz

      DO I=1,NVEC(NBOX)
c         print *,'INVEC=',i
         TMPVECP1=VECPAIR(1,I,NBOX)
         TMPVECP2=ISIDE*VECPAIR(2,I,NBOX)  

c         print *,'vec1=',tmpvecp1
c         print *,'vec2=',tmpvecp2         

C     IS POINT(R1,R2,R3) ON THE RIGHT SIDE OF THE PLANE
C     POINT1(PX,PY,PZ) MUST BE CALCULATED
c     take into account that plpoint(3,*,*) is always ZERO !!!
         PX=P1+PLPOINT(1,I,NBOX)*DVEC(1,1)+PLPOINT(2,I,NBOX)*DVEC(2,1)
         PY=P2+PLPOINT(1,I,NBOX)*DVEC(1,2)+PLPOINT(2,I,NBOX)*DVEC(2,2)
         PZ=P3+PLPOINT(1,I,NBOX)*DVEC(1,3)+PLPOINT(2,I,NBOX)*DVEC(2,3)
c         print *,'PXYZ> ',px,py,pz
C     ****************************
C     *** handle BOXES 1,3,5,6 ***
C     ****************************
         if(nbox.ne.2.and.nbox.ne.4)then               
            if(.not.isinbox(nbox,tmpvecp1,tmpvecp2,
     1           px,py,pz,r1,r2,r3,dvec)) box=.false.
         endif
C     **********************************
C     *** *** handle BOXES 2 & 4 *** ***
C     **********************************
         if(nbox.eq.2.or.nbox.eq.4)then
C     *** "BOX2" ****
            if(nbox.eq.2)then
C     is r1,r2,r3 in the 3th part of the box; 
C     for 3th part i must be 1 -> this means: first point and first vector
               if(i.eq.1.and.xx.ge.CON23)then
c                  print *,'PART> 3th',xx
                  if(.not.isinbox(nbox,tmpvecp1,tmpvecp2,
     1             px,py,pz,r1,r2,r3,dvec)) box=.false.
               endif
C     is r1,r2,r3 in the 2nd part of the box
               if(i.eq.2.and.xx.ge.CON13.and.
     1              xx.lt.CON23)then
c                  print *,'PART> 2nd',xx
                  if(.not.isinbox(nbox,tmpvecp1,tmpvecp2,
     1                 px,py,pz,r1,r2,r3,dvec)) box=.false.
               endif
C     is r1,r2,r3 in the first part of the box
               if(i.eq.3.and.xx.lt.CON13)then
c                  print *,'PART> 1st'
                  if(.not.isinbox(nbox,tmpvecp1,tmpvecp2,
     1                 px,py,pz,r1,r2,r3,dvec)) box=.false.
               endif   
C     *** "BOX4" ***
            elseif(nbox.eq.4)then
C     is r1,r2,r3 in the 3th part of the box;
C     for 3th part i must be 3 -> this means third point & 3th vector
               if(i.eq.3.and.yy.ge.CON23)then
c                  print *,'PART> 3th',i,yy
                  if(.not.isinbox(nbox,tmpvecp1,tmpvecp2,
     1                 px,py,pz,r1,r2,r3,dvec)) box=.false.
               endif
C     is r1,r2,r3 in the 2nd part of the box
               if(i.eq.2.and.yy.ge.CON13.and.
     1              yy.lt.CON23)then
c                  print *,'PART> 2nd',i,yy
                  if(.not.isinbox(nbox,tmpvecp1,tmpvecp2,
     1                 px,py,pz,r1,r2,r3,dvec)) box=.false.
               endif
C     is r1,r2,r3 in the first part of the box
               if(i.eq.1.and.yy.lt.CON13)then
c                  print *,'PART> 1st',i,yy
                  if(.not.isinbox(nbox,tmpvecp1,tmpvecp2,
     1                 px,py,pz,r1,r2,r3,dvec)) box=.false.
               endif   
            endif               !BOX4
         endif                  !BOX4.OR.BOX2
      enddo                     !DO I=1,NVEC(NBOX)
     
      RETURN
      END

 
      LOGICAL FUNCTION ISINBOX(NBOX,NVC1,NVC2,P1,P2,P3,R1,R2,R3,dvec)
      REAL*8 P1,P2,P3,R1,R2,R3
      REAL*8 DVEC(3,3),MINTOL
      REAL*8 A,B,C    !components of normal vector
      REAL*8 D,KFC   ! Ax + By + Cz + D = 0
      REAL*8 AVEC,BVEC,CVEC,AVEC2,BVEC2,CVEC2
C     if k<0 negative=.true.; third vector is turned upside down
      LOGICAL NEGATIVE  
      REAL*8 PLPOINT(3,4,6),PLVEC(3,3,6)
      INTEGER NVEC(6),VECPAIR(2,3,6)
      
      COMMON/BOXS/ PLVEC,NVEC,VECPAIR,PLPOINT

      PARAMETER (MINTOL=1.0d-6)
     
      NVEC1=NVC1
      NVEC2=NVC2
c      print *,'nvec1=',nvec1,'nvec2=',nvec2
      NEGATIVE=.FALSE.
      ISINBOX=.FALSE.
      IF(NVEC2.LT.0.)THEN
         NVEC2=-NVEC2
         NEGATIVE=.TRUE.
      ENDIF

C     plvec is in coloumn major mode, but dvec is in row-major mode
C     *** transform plvec1 in XYZ basis ***
      AVEC=DVEC(1,1)*PLVEC(1,NVEC1,NBOX)+DVEC(2,1)*PLVEC(2,NVEC1,NBOX)+
     1     DVEC(3,1)*PLVEC(3,NVEC1,NBOX)
      BVEC=DVEC(1,2)*PLVEC(1,NVEC1,NBOX)+DVEC(2,2)*PLVEC(2,NVEC1,NBOX)+
     1     DVEC(3,2)*PLVEC(3,NVEC1,NBOX)
      CVEC=DVEC(1,3)*PLVEC(1,NVEC1,NBOX)+DVEC(2,3)*PLVEC(2,NVEC1,NBOX)+
     1     DVEC(3,3)*PLVEC(3,NVEC1,NBOX)
      A=BVEC*DVEC(3,3)-DVEC(3,2)*CVEC
      B=CVEC*DVEC(3,1)-DVEC(3,3)*AVEC
      C=AVEC*DVEC(3,2)-DVEC(3,1)*BVEC

C     *** transform plvec2 in XYZ basis ***
      AVEC2=DVEC(1,1)*PLVEC(1,NVEC2,NBOX)+DVEC(2,1)*PLVEC(2,NVEC2,NBOX)+
     1     DVEC(3,1)*PLVEC(3,NVEC2,NBOX)
      BVEC2=DVEC(1,2)*PLVEC(1,NVEC2,NBOX)+DVEC(2,2)*PLVEC(2,NVEC2,NBOX)+
     1     DVEC(3,2)*PLVEC(3,NVEC2,NBOX)
      CVEC2=DVEC(1,3)*PLVEC(1,NVEC2,NBOX)+DVEC(2,3)*PLVEC(2,NVEC2,NBOX)+
     1     DVEC(3,3)*PLVEC(3,NVEC2,NBOX)
      IF(A*AVEC2+B*BVEC2+C*CVEC2.LT.0.0d0)THEN
         A=-A
         B=-B
         C=-C
      ENDIF

c      print *,'PLVEC:',PLVEC(1,NVEC1,NBOX),PLVEC(2,NVEC1,NBOX),
c     1     PLVEC(3,NVEC1,NBOX)*DVEC(3,1)
c      print *,'DVEC',DVEC(3,1),DVEC(3,2),DVEC(3,3)
c      print *,'p1,p2,p3=',p1,p2,p3
c      print *,'r1,r2,r3=',r1,r2,r3
      D=-A*P1-B*P2-C*P3
c      print *,'a,b,c,d::',a,b,c,d
      KFC=A*R1+B*R2+C*R3+D
c      print *,'kfc1=',kfc,', kfc2=',A*R1+B*R2+C*R3-(-A*P1-B*P2-C*P3)
      IF(KFC.GE.(0.0d0-MINTOL).AND.(.NOT.NEGATIVE)) ISINBOX=.TRUE.
      IF(KFC.LE.(0.0d0+MINTOL).AND.NEGATIVE) ISINBOX=.TRUE.
c      print *,"ISINBOX=",isinbox

      RETURN
      END


      
