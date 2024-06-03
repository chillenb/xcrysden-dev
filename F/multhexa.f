      SUBROUTINE MULTHEXA(TTC,IGROUP,NATOM)
C     TTC...DATA FOR HTTC OR RTTC
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
c     **********************************************************************
c     in order to allow compiling on non-g77 compilers; replace NATOM arrays
c     with INATOM parameter
c     **********************************************************************
      REAL*8 TTC(3,13),TOL,XC(INATOM),YC(INATOM),ZC(INATOM),
     *     AAB(2),BAB(2),X(NAC),Y(NAC),Z(NAC),DVEC(3,3),A(3)
      INTEGER ITTC(13,6),NTTC(6),CONTYP(13),FLAG(INATOM),NAT(NAC)
      LOGICAL EQUAL,ISINSIDE,TTCBOX,LPRINT
      REAL*8 AC(3,2),BC(3,2),CC(3,2),FC(3,4),IC(3,2),RC(3,3),HC(3,3),
     1     DVC(3,3),IDVC(3,3)
      
      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR
      COMMON/DIR/ NXDIR,NYDIR,NZDIR
      COMMON/APOS/ AC,BC,CC,FC,IC,HC,RC,DVC

      PARAMETER (TOL=1.0d-8)

c      print *,'MULTHEXA'
C     THERE IS FIVE ARBITRARY DIFFERENT PIECES OF HTTC/RTTC 
C     FROM THIS FIVE PIECES IT'S POSSIBLE TO PUT TOGETHER THE CRYSTAL
C     PIECE #1: 4,5,6th --> end; end is last point (point 13)
C     PIECE #2: 3,4,5,6th --> end
C     PIECE #3: 1--> end
C     PIECE #4: 1,5,6 --> end
C     PIECE #5: 5,6 --> end
     
C     there is also 6th piece, whose purpose is filling the crystal so,
C     that when cutting on out borders, it will be cut-out nicely 
C     PIECE #6: 7 --> end

C     NTTC...TELLS HOW MANY ATOMS IS IN EVERY PIECE
      DATA NTTC/10,11,13,10,9,7/
      DATA ITTC/
     1     4, 5, 6, 7, 8, 9,10,11,12,13, 0, 0, 0,
     2     3, 4, 5, 6, 7, 8, 9,10,11,12,13, 0, 0,
     3     1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,
     4     1, 5, 6, 7, 8, 9,10,11,12,13, 0, 0, 0,
     5     5, 6, 7, 8, 9,10,11,12,13, 0, 0, 0, 0,
     6     7, 8, 9,10,11,12,13, 0, 0, 0, 0, 0, 0/
C     THIS IS FOR CELL-FRAMES BASED CONNECTIVITY
C     CONNECTIVITY TYPES:
C     0...NO CONNECTIVITY
C     1...OUTER CONNECTIVITY
C     2...INNER CONNECTIVITY
C     3...O/I CONNECTIVITY
C     TO EACH POSITION WE ASSIGN CORRESPONDING FLAG

c     DATA CONTYP/3,1,3,1,3,1,2,0,0,0,0,0,0/
c     Tue Apr  8 14:32:28 METDST 1997;; changed to
      DATA CONTYP/1,3,1,3,1,3,2,0,0,0,0,0,0/

C     THIS IS FOR PROPAGATION DIRECTION OF CRYSTAL GROWTH
C     IF PROPAGATION IS IN A,B PLANE THEN PROPAGATION DIRECTIONS ARE:
C     A+(A+B)->AAB & B+(A+B)->BAB
      DATA AAB/2.0d0,1.0d0/
      DATA BAB/1.0d0,2.0d0/

      CALL TTCBOXES(TTC)

c     idv will be used for calculating fractional coorinates
      call InvertSqMat123(dvec, idvc, 3)

C     ************************************************
C     ******           ATOMS SECTION            ******
C     ************************************************
C     WRITE "ATOMS" TO FILE2, TO MARK THE BEGINNING OF "ATOMS" SECTION
      WRITE(12,*) 'ATOMS'

C     IF NXDIR=0 OR NYDIR=0 -> CALL MULTATOM
      IF((NXDIR.EQ.0.OR.NYDIR.EQ.0).AND.IGROUP.EQ.7)
     *     CALL MULTATOM(RC,3,IGROUP,NATOM,.false.,3,1)
      IF((NXDIR.EQ.0.OR.NYDIR.EQ.0).AND.IGROUP.EQ.8)
     *     CALL MULTATOM(HC,3,IGROUP,NATOM,.false.,3,1)

      DNX=DBLE(NXDIR)
      DNY=DBLE(NYDIR)
      DNZ=DBLE(NZDIR)

c      print *,'MULTHEXA'
      N=0
      DO I=0,NXDIR+1
         DO J=0,NYDIR+1
            DO K=0,NZDIR              
C     DETERMINE WHICH PIECE TO ADD
               IF(I.EQ.0.OR.I.GT.NXDIR.OR.
     1              J.EQ.0.OR.J.GT.NYDIR)THEN
                  IPIECE=6               
               ELSEIF(I.EQ.NXDIR.AND.J.EQ.NYDIR)THEN
                  IPIECE=3
               ELSEIF(I.EQ.1.AND.J.LT.NYDIR)THEN
                  IPIECE=1
               ELSEIF(I.EQ.NXDIR.AND.J.LT.NYDIR)THEN
                  IPIECE=4
               ELSEIF(I.LT.NXDIR.AND.J.EQ.NYDIR)THEN
                  IPIECE=2
               ELSE
                  IPIECE=5
               ENDIF
c               print *,'ipiece=',ipiece               
C     CALCULATE THE ORIGIN OF THE PIECE TO BE ADDED
               DI=DBLE(I-1)
               DJ=DBLE(J-1)
               DK=DBLE(K)
c               print *,'di,dj,dk',di,dj,dk
               AORIG=AAB(1)*DI+BAB(1)*DJ
               BORIG=AAB(2)*DI+BAB(2)*DJ
               CORIG=DK
               
               DO II=1,NTTC(IPIECE)
                  ACMP=AORIG+TTC(1,ITTC(II,IPIECE))
                  BCMP=BORIG+TTC(2,ITTC(II,IPIECE))
                  CCMP=CORIG+TTC(3,ITTC(II,IPIECE))

C     CHECK IF "CCMP" IS "OUT OF RANGE"; THAT IS, WE WANT THAT
C     OUR CRYSTAL IS CUT OUT NICELY ON THE OUT-BORDERS
C     IF "OUT OF RANGE" --> DELETE POINT
                  XPOS=ACMP*DVEC(1,1)+BCMP*DVEC(2,1)+CCMP*DVEC(3,1)
                  YPOS=ACMP*DVEC(1,2)+BCMP*DVEC(2,2)+CCMP*DVEC(3,2)
                  ZPOS=ACMP*DVEC(1,3)+BCMP*DVEC(2,3)+CCMP*DVEC(3,3)
C     *** *** FRAMES *** ***
                  if((i.gt.0.and.i.le.nxdir).and.
     1                 (j.gt.0.and.j.le.nydir))then
                     N=N+1
                     XC(N)=XPOS
                     YC(N)=YPOS
                     ZC(N)=ZPOS
C     FOURTH COORDINATE IS FOR CONNECTIVITY
                     FLAG(N)=CONTYP(ITTC(II,IPIECE))                  
                  endif

                  DO JJ=1,NATR
                     lprint=.true.
c     *****************************************************
c     COORDINATES must be within [-1,1], check that !!!
c     *****************************************************
                     rx=x(jj)*idvc(1,1)+y(jj)*idvc(2,1)+z(jj)*idvc(3,1)
                     ry=x(jj)*idvc(1,2)+y(jj)*idvc(2,2)+z(jj)*idvc(3,2)
                     rz=x(jj)*idvc(1,3)+y(jj)*idvc(2,3)+z(jj)*idvc(3,3)
                     A(1)=ACMP - AORIG + RX
                     A(2)=BCMP - BORIG + RY
                     A(3)=CCMP - CORIG + RZ
                     do ia=1,3
                        if ( a(ia) .gt. 1.0d0 ) then
                           a(ia) = a(ia) - 1.0d0
                        else if ( a(ia) .lt. -1.0d0 ) then
                           a(ia) = a(ia) + 1.0d0
                        endif
                     enddo   
                     a(1)=a(1) + AORIG
                     a(2)=a(2) + BORIG
                     a(3)=a(3) + CORIG
                     XX=a(1)*DVEC(1,1)+a(2)*DVEC(2,1)+a(3)*DVEC(3,1)
                     YY=a(1)*DVEC(1,2)+a(2)*DVEC(2,2)+a(3)*DVEC(3,2)
                     ZZ=a(1)*DVEC(1,3)+a(2)*DVEC(2,3)+a(3)*DVEC(3,3)

c     print *,'NEW ATOM'
C     ***********************************************************************
C     WE WANT THAT CRYSTAL IS CUT OUT NICELY ON THE OUT-BORDER
C     IF .not. "OUT OF RANGE" -> GO
C     *** *** AB planes *** ***
                     if(k.eq.0.or.k.eq.nzdir)then
C     *** lower AB plane ***
                        if(.not.isinside(dvec,1,2,3,
     1                       0.0d0,0.0d0,0.0d0,xx,yy,zz)) 
     *                       lprint=.false.
C     *** upper AB plane ***
                        if(.not.isinside(dvec,1,2,-3,
     1                       dnz*dvec(3,1),
     2                       dnz*dvec(3,2),
     3                       dnz*dvec(3,3),xx,yy,zz)) lprint=.false.
                     endif
                     
                     if(i.eq.0) dii=di+1.0d0
                     if(i.gt.0.and.i.le.nxdir) dii=di
                     if(j.eq.0) djj=dj+1.0d0
                     if(j.gt.0.and.j.le.nydir) djj=dj
                     if(i.gt.nxdir) DI=DI-1.0d0
                     if(j.gt.nydir) DJ=DJ-1.0d0

                     AORIG2=AAB(1)*DII+BAB(1)*DJJ
                     BORIG2=AAB(2)*DII+BAB(2)*DJJ
                     p1=aorig2*dvec(1,1)+borig2*dvec(2,1)
                     p2=aorig2*dvec(1,2)+borig2*dvec(2,2)
                     p3=aorig2*dvec(1,3)+borig2*dvec(2,3)
C     *** AC & BC PLANES *** 
                     if(i.le.1.and.j.le.1)then
c     print *,'BOX1'
                        if(.not.ttcbox(1,p1,p2,p3,xx,yy,zz,dvec))
     1                       lprint=.false. 
                     endif
                     if((i.gt.1.and.i.lt.nxdir).and.j.le.1)then
c                           print *,'BOX2'
                        if(.not.ttcbox(2,p1,p2,p3,xx,yy,zz,dvec))
     1                       lprint=.false.
                     endif
                     if(i.ge.nxdir.and.j.le.1)then
c     print *,'BOX3'
                        if(.not.ttcbox(3,p1,p2,p3,xx,yy,zz,dvec))
     1                       lprint=.false.
                     endif
                     if(i.ge.nxdir.and.(j.gt.1.and.j.lt.nydir))then
c     print *,'BOX4'
                        if(.not.ttcbox(4,p1,p2,p3,xx,yy,zz,dvec))
     1                       lprint=.false.
                     endif
                     if(i.ge.nxdir.and.j.ge.nydir)then
c     print *,'BOX5>',p1,p2,p3,'i,j=',i,j
                        if(.not.ttcbox(5,p1,p2,p3,xx,yy,zz,dvec))
     1                       lprint=.false.
                     endif
                     if((i.gt.1.and.i.lt.nxdir).and.j.ge.nydir)then
c     print *,'BOX6'
                        if(.not.ttcbox(6,p1,p2,p3,xx,yy,zz,dvec))
     1                       lprint=.false.
                     endif 
                     if(i.le.1.and.j.ge.nydir)then
c     print *,'BOX7'
                        if(.not.ttcbox(7,p1,p2,p3,xx,yy,zz,dvec))
     1                       lprint=.false.
                     endif
                     if(i.le.1.and.(j.gt.1.and.j.lt.nydir))then
c     print *,'BOX8'
                        if(.not.ttcbox(8,p1,p2,p3,xx,yy,zz,dvec))
     1                       lprint=.false.
                     endif
                     
C     *** if (lprint) then print coor 
                     IF(LPRINT)THEN
C     **** ** write coor ** ****
                        WRITE(12,21) NAT(JJ),XX,YY,ZZ
 21                     FORMAT(I5,3F16.10)
                     ENDIF
                  ENDDO                  
               ENDDO !DO II=1,NTTC(IPIECE)
            ENDDO !K=
         ENDDO !J=
      ENDDO !I=
      
C     *********************************************************
C     ********     FRAMES --- FRAMES --- FRAMES      **********
C     *********************************************************
      IF (N.GT.0) WRITE(12,*) 'FRAMES'
      ADIS=VECSIZE3(DVEC,1)
      CDIS=VECSIZE3(DVEC,3)
      DO I=1,N-1
         DO J=I+1,N
            DIS=DIST(XC(I),YC(I),ZC(I),XC(J),YC(J),ZC(J))
            write(20,*) dis, adis, cdis
            IF(EQUAL(DIS,ADIS)) THEN
               write(21,*) dis,adis
               IF( (FLAG(I).EQ.1.AND.FLAG(J).EQ.3) .OR. 
     *              (FLAG(I).EQ.3.AND.FLAG(J).EQ.1) )
     *              WRITE(12,23) '1',XC(I),YC(I),ZC(I),XC(J),YC(J),ZC(J)
               IF( (FLAG(I).EQ.1.AND.FLAG(J).EQ.2) .OR.
     *              (FLAG(I).EQ.2.AND.FLAG(J).EQ.1) )
     *              WRITE(12,23) '2',XC(I),YC(I),ZC(I),XC(J),YC(J),ZC(J)
            ELSEIF(EQUAL(DIS,CDIS)) THEN
               write(21,*) dis,cdis
               IF( (FLAG(I).EQ.1.AND.FLAG(J).EQ.1) .OR.
     *              (FLAG(I).EQ.3.AND.FLAG(J).EQ.3) )
     *              WRITE(12,23) '1',XC(I),YC(I),ZC(I),XC(J),YC(J),ZC(J)
            ENDIF
         ENDDO
      ENDDO
 23   FORMAT(A1,2X,3F16.10,2X,3F16.10)
      RETURN
      END


