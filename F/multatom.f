      SUBROUTINE MULTATOM(CELLC,INCELL,IGROUP,NATOM,TPH,IDIM,IMODE2)
C     TPH...FLAG FOR TRIPLE-PRIMITIV-HEXAGONAL CELL
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
c     **********************************************************************
c     in order to allow compiling on non-g77 compilers; replace NATOM arrays
c     with INATOM parameter
c     **********************************************************************
      LOGICAL TPH
      REAL*8 X(NAC),Y(NAC),Z(NAC),DVEC(3,3),CELLC(3,INCELL),
     *     XC(6*INATOM),YC(6*INATOM),ZC(6*INATOM),
     *     XTPH(INATOM),YTPH(INATOM),ZTPH(INATOM),TPHFRAME(3,3),
     *     IDVC(3,3), FV(NAC,3), small_number
      INTEGER TPHFLAG(INATOM)
      INTEGER NAT(NAC) !I THINK 100 ATOMS PER CELL IS MORE THAN ENOUGH
      LOGICAL EQUAL,LPRINT,ISINSIDE,BOX
      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR
      COMMON/DIR/ NXDIR,NYDIR,NZDIR
      COMMON/FV/ FV,NFV
      PARAMETER (CON13=1.0d0/3.0d0, CON23=2.0d0/3.0d0, 
     $     small_number=1.0d-10)

      include 'mode.inc'

C     *************************************************************************
C     **** data for TPH frames ***
      DATA TPHFRAME/
     1     0.0d0, 0.0d0, 0.0d0,
     2     CON23, CON13, 0.0d0,
     3     CON13, CON23, 0.0d0/

c      print *,'MULTATOM> I am IN'
c      print *,'MULTATOM> NATR,NCELL',natr,ncell
c      print *,'INCELL,IGROUP',incell,igroup
c      print *,'NXDIR, NYDIR, NZDIR::',nxdir,nydir,nzdir
c      print *,'CELLC:: ',((CELLC(i,j),i=1,3),j=1,incell)
c      print *,'DVEC:: ',((DVEC(i,j),i=1,3),j=1,3)
c      print *,'NAT::',(NAT(i),i=1,natr)
c      print *,'X::',(X(i),i=1,natr)
c      print *,'Y::',(Y(i),i=1,natr)
c      print *,'Z::',(Z(i),i=1,natr)

      print *,'Entering multatom routine'
c     idv will be used for calculating fractional coorinates
      call InvertSqMat123(dvec, idvc, idim)

      NX=NXDIR
      NY=NYDIR
      NZ=NZDIR      

C     do we have M2_CELL or M2_TR_ASYM_UNIT???
C     if M2_TR_ASYM_UNIT then reduce (NXDIR,NYDIR,NZDIR) by 1
      IF (IMODE2.EQ.M2_TR_ASYM_UNIT)THEN
         NXDIR=NXDIR-1
         NYDIR=NYDIR-1
         NZDIR=NZDIR-1
      ENDIF
      
      IF(IDIM.LT.3) THEN
         NZDIR=0
         NZ   =0
      ENDIF
      IF(IDIM.LT.2) THEN
         NYDIR=0
         NY   =0
      ENDIF
      IF(IDIM.LT.1) THEN
         NXDIR=0
         NX   =0
      ENDIF
C     **** more data for TPH; crystal cells cut-out DATA
      CALL BOXES
      
C     ************************************************
C     ******           ATOMS SECTION            ******
C     ************************************************
C     WRITE "ATOMS" TO FILE2; MARK THE BEGINNING OF "ATOMS" SECTION
      WRITE(12,*) 'ATOMS'

C     THIS IS USED FOR CUT-OFF CRITERIA (LOOK BELOW)
      DNX=0.0d0
      DNY=0.0d0
      DNZ=0.0d0
      IF(NXDIR.GT.0 .and. IDIM.GT.0)DNX=DBLE(NXDIR)
      IF(NXDIR.EQ.0 .and. IDIM.GT.0)DNX=1.0d0
      IF(NYDIR.GT.0 .and. IDIM.GT.1)DNY=DBLE(NYDIR)
      IF(NYDIR.EQ.0 .and. IDIM.GT.1)DNY=1.0d0
      IF(NZDIR.GT.0 .and. IDIM.GT.2)DNZ=DBLE(NZDIR)
      IF(NZDIR.EQ.0 .and. IDIM.GT.2)DNZ=1.0d0

      DO I=0,NXDIR
         DO J=0,NYDIR
            DO K=0,NZDIR
               DI=DBLE(I)
               DJ=DBLE(J)
               DK=DBLE(K)
               DO II=1,NATR
                  rx=x(ii)*idvc(1,1) + y(ii)*idvc(2,1) + z(ii)*idvc(3,1)
                  ry=x(ii)*idvc(1,2) + y(ii)*idvc(2,2) + z(ii)*idvc(3,2)
                  rz=x(ii)*idvc(1,3) + y(ii)*idvc(2,3) + z(ii)*idvc(3,3)
                  DO IN=1,INCELL
C     lprint is logical parameter "for crystal-cut-out"
                     lprint=.true.  
c     *****************************************************
c     COORDINATES must be within [-0.5,0.5), check that !!!
c     *****************************************************

c     There was a roundoff problem for as simple structure as NaCl
c     (problem reported by W. YU). To avoid this round-off problem let's
c     add a small number
                     ACMP=CELLC(1,IN)+rx + small_number
                     BCMP=CELLC(2,IN)+ry + small_number
                     CCMP=CELLC(3,IN)+rz + small_number
                     print *,'Before: Xcmp:', nat(ii), acmp, bcmp, ccmp
                     IF(IMODE2.EQ.M2_CELL) then
                        if (idim.gt.0) ACMP=ACMP-DINT(ACMP)
                        if (idim.gt.1) BCMP=BCMP-DINT(BCMP)
                        if (idim.gt.2) CCMP=CCMP-DINT(CCMP)
                     endif
c     now we subtract that small number
                     ACMP=ACMP - small_number
                     BCMP=BCMP - small_number
                     CCMP=CCMP - small_number

                     print *,' after: Xcmp:', nat(ii), acmp, bcmp, ccmp
                     ACMP=ACMP+DI
                     BCMP=BCMP+DJ
                     CCMP=CCMP+DK
c                     print *,'MULTATOM> ACMP,BCMP,CCMP',ACMP,BCMP,CCMP
                     XX=ACMP*DVEC(1,1)+BCMP*DVEC(2,1)+CCMP*DVEC(3,1)
                     YY=ACMP*DVEC(1,2)+BCMP*DVEC(2,2)+CCMP*DVEC(3,2)
                     ZZ=ACMP*DVEC(1,3)+BCMP*DVEC(2,3)+CCMP*DVEC(3,3)
c                     XX=XPOS+X(II)
c                     YY=YPOS+Y(II)
c                     ZZ=ZPOS+Z(II)                     

C     if M2_TR_ASYM_UNIT then SKIP CUTTING
                     IF(IMODE2.EQ.M2_TR_ASYM_UNIT) GOTO 99

C     *************************************************************************
C     WE WANT THAT CRYSTALL IS CUT OUT NICELY ON THE OUT-BORDERS
C     IF .not."OUT OF RANGE" -> GO

C     *** CUTTING for .NOT.TPH ***
                     if(.not.tph.and.(
     1                    (i.eq.0.or.i.eq.nxdir).or.
     1                    (j.eq.0.or.j.eq.nydir).or.
     2                    (k.eq.0.or.k.eq.nzdir)) ) then

c     *** if we are dealing with SLAB, we don't perform AB cutting
                        if(idim.eq.3)then
c     *** *** lower AB plane *** ***
                           if(.not.isinside(dvec,1,2,3,0.0d0,0.0d0,
     1                       0.0d0,xx,yy,zz)) lprint=.false.
c     *** *** upper AB plane
                           if(.not.isinside(dvec,1,2,-3,
     1                          dnz*dvec(3,1),
     2                          dnz*dvec(3,2),
     3                          dnz*dvec(3,3),xx,yy,zz)) lprint=.false.
                        endif
c     *** *** lower AC plane
                        if(.not.isinside(dvec,1,3,2,0.0d0,0.0d0,0.0d0,
     1                       xx,yy,zz)) lprint=.false.
c     *** *** upper AC plane
                        if(.not.isinside(dvec,1,3,-2,
     1                       dny*dvec(2,1),
     2                       dny*dvec(2,2),
     3                       dny*dvec(2,3),xx,yy,zz)) lprint=.false.
c     *** *** lower BC plane
                        if(.not.isinside(dvec,2,3,1,0.0d0,0.0d0,0.0d0,
     1                       xx,yy,zz)) lprint=.false.
c     *** *** upper BC plane
                        if(.not.isinside(dvec,2,3,-1,
     1                       dnx*dvec(1,1),
     2                       dnx*dvec(1,2),
     3                       dnx*dvec(1,3),xx,yy,zz)) lprint=.false.
                     endif
c     *** *** if (lprint=false).or.(polymer(idim=1)).or.
c                 (imode2=m2_tr_asym_unit) -> print atom coordinates
 99                  continue
c     ********************************************************************
c     t.k: this is bug for (idim.eq.1), since it should be cutted accoring
c     to dvec(1,*) criteria
c     ********************************************************************

                     print *,'nat, xx, yy, zz, accepted:',
     $                    nat(ii),xx,yy,zz,lprint                        
                     
                     if((idim.eq.1)) then
                        if ( (XX.gt.-1.0d-6 .AND.
     $                       XX.lt.(dvec(1,1)*dnx+1.0d-6) .AND.
     $                       imode2.eq.m2_cell) .OR.
     $                       imode2.eq.m2_tr_asym_unit ) then
                           WRITE(12,21) NAT(II),XX,YY,ZZ,
     $                          (FV(II,ifv),ifv=1,nfv)
                        endif
                     else if((lprint.and.(.not.tph)).or.
     1                       (imode2.eq.m2_tr_asym_unit)) then
                        WRITE(12,21) NAT(II),XX,YY,ZZ,
     $                       (FV(II,ifv),ifv=1,nfv)
                     endif
C     *************************************************************************
C     ***         TPH CUTTING; TPH cutting is just for crystals !!!!        ***
                     IF(TPH.and.(idim.eq.3)) THEN                        
C     TAKE A SPECIAL CARE FOR TRIPLE-PRIMITIV HEXAGONAL CELL(S)
C     *** AB planes ***
c                        print *,'IF(TPH) THEN'
                        if(k.eq.0.or.k.eq.nzdir)then
c     *** *** lower AB plane *** ***
                           if(.not.isinside(dvec,1,2,3,
     1                          0.0d0,0.0d0,0.0d0,
     2                          xx,yy,zz)) lprint=.false.
c     *** *** upper AB plane
                           if(.not.isinside(dvec,1,2,-3,
     1                          dnz*dvec(3,1),
     2                          dnz*dvec(3,2),
     3                          dnz*dvec(3,3),xx,yy,zz)) lprint=.false. 
                        endif

C     *** AC & BC PLANES ***
                        IF(I.EQ.0.AND.J.EQ.0)THEN 
c                           print *,'BOX1'
                           IF(.NOT.BOX(1,0.0d0,0.0d0,0.0d0,
     1                          xx,yy,zz,dvec)) 
     *                          LPRINT=.FALSE.
                        ENDIF
                        IF(I.EQ.NXDIR.AND.J.EQ.NYDIR)THEN
c                           print *,'BOX5'
                           IF(.NOT.BOX(5,
     1                          dnx*dvec(1,1)+dny*dvec(2,1),
     2                          dnx*dvec(1,2)+dny*dvec(2,2),
     3                          dnx*dvec(1,3)+dny*dvec(2,3),
     4                          xx,yy,zz,dvec)) 
     *                       LPRINT=.FALSE.
                        ENDIF

                        IF(J.EQ.0.AND.(I.GT.0.AND.I.LT.NXDIR))THEN
c                           print *,'BOX2'
                           IF(.NOT.BOX(2,
     1                          di*dvec(1,1),di*dvec(1,2),di*dvec(1,3),
     2                          xx,yy,zz,dvec))
     *                          LPRINT=.FALSE.
                        ENDIF
                        IF(J.EQ.NYDIR.AND.(I.GT.0.AND.I.LT.NXDIR))THEN
c                           print *,'BOX-2'
                           IF(.NOT.BOX(-2,
     1                          di*dvec(1,1)+dny*dvec(2,1),
     2                          di*dvec(1,2)+dny*dvec(2,2),
     3                          di*dvec(1,3)+dny*dvec(2,3),
     4                          xx,yy,zz,dvec)) 
     *                          LPRINT=.FALSE.
                        ENDIF

                        IF(I.EQ.NXDIR.AND.J.EQ.0)THEN
C                           print *,'BOX3'
                           IF(.NOT.BOX(3,
     1                          dnx*dvec(1,1),dnx*dvec(1,2),
     2                          dnx*dvec(1,3),xx,yy,zz,dvec))
     *                          LPRINT=.FALSE.
                        ENDIF
                        IF(J.EQ.NYDIR.AND.I.EQ.0)THEN
c                           print *,'BOX6'
                           IF(.NOT.BOX(6,
     1                          dny*dvec(2,1),dny*dvec(2,2),
     2                          dny*dvec(2,3),xx,yy,zz,dvec))
     *                          LPRINT=.FALSE.
                        ENDIF

                        IF(I.EQ.NXDIR.AND.(J.GT.0.AND.J.LT.NYDIR))THEN
c                           print *,'BOX4'
                           IF(.NOT.BOX(4,
     1                          dnx*dvec(1,1)+dj*dvec(2,1),
     2                          dnx*dvec(1,2)+dj*dvec(2,2),
     3                          dnx*dvec(1,3)+dj*dvec(2,3),
     4                          xx,yy,zz,dvec))
     *                          LPRINT=.FALSE.
                        ENDIF
                        IF(I.EQ.0.AND.(J.GT.0.AND.J.LT.NYDIR))THEN
c                           print *,'BOX-4'
                           IF(.NOT.BOX(-4,
     1                          dj*dvec(2,1),dj*dvec(2,2),dj*dvec(2,3),
     2                          xx,yy,zz,dvec))
     *                          LPRINT=.FALSE.
                        ENDIF
C     if lprint=.true. write coordinates
                        IF(LPRINT) WRITE(12,21) NAT(II),XX,YY,ZZ,
     $                       (FV(II,ifv),ifv=1,nfv)
                     ENDIF
 21                  FORMAT(I5,6F16.10)                                       
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      
C     ************************************************
C     ******          FRAMES SECTION            ******
C     ************************************************
      N=0
      NTPH=0
      DO I=0,NX
         DO J=0,NY
            DO K=0,NZ               
               DI=DBLE(I)
               DJ=DBLE(J)
               DK=DBLE(K)
               IF(.NOT.TPH)THEN

c                  print *,'NATOM,N,',natom,n
c                  print *,'nxdir,nydir,nzdir::',i,j,k

C     *** THIS IS FOR CRYSTAL FRAMES; NOT FOR TPH CRYSTAL FRAMES ***
C     if i = nx than there is no frame in "A"-direction, because this
C     frame line would be from nx->nx+1 and thats out of range
                  IF(I.LT.NX) CALL MAKEFRAME(NATOM,N,DVEC,
     *                 XC,YC,ZC,DI,DJ,DK,DI+1,DJ,DK)
                  IF(J.LT.NY) CALL MAKEFRAME(NATOM,N,DVEC,
     *                 XC,YC,ZC,DI,DJ,DK,DI,DJ+1,DK)
                  IF(K.LT.NZ) CALL MAKEFRAME(NATOM,N,DVEC,
     *                 XC,YC,ZC,DI,DJ,DK,DI,DJ,DK+1)                  
               ELSE
C     ******************************************
C     *** NOW MAKE ALSO THP'S CRYSTAL FRAMES ***
C     ******************************************
                  DO IN=1,INCELL
C     *** frames are "written" only if this condition is .TRUE.
                     IF((.NOT.(IN.EQ.1.AND.(I.EQ.0.OR.J.EQ.0)))
     *                    .AND.(.NOT.(IN.EQ.2.AND.I.EQ.NXDIR.AND.
     *                    J.EQ.0))
     *                    .AND.(.NOT.(IN.EQ.3.AND.J.EQ.NYDIR.AND.
     *                    I.EQ.0)))THEN
                        NTPH=NTPH+1
                        ACMP=TPHFRAME(1,IN)+DI
                        BCMP=TPHFRAME(2,IN)+DJ
                        CCMP=TPHFRAME(3,IN)+DK
                        XPOS=ACMP*DVEC(1,1)+BCMP*DVEC(2,1)
     1                       +CCMP*DVEC(3,1)
                        YPOS=ACMP*DVEC(1,2)+BCMP*DVEC(2,2)
     1                       +CCMP*DVEC(3,2)
                        ZPOS=ACMP*DVEC(1,3)+BCMP*DVEC(2,3)
     1                       +CCMP*DVEC(3,3)
                        XTPH(NTPH)=XPOS
                        YTPH(NTPH)=YPOS
                        ZTPH(NTPH)=ZPOS
                        IF(IN.EQ.1) TPHFLAG(NTPH)=2
                        IF(IN.GT.1) TPHFLAG(NTPH)=1
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO

C     NOW MARK THE BEGINNING OF "FRAMES" SECTION
      IF(N.GT.0 .OR. NTPH.GT.0) WRITE(12,*) 'FRAMES'
      
      IF(.NOT.TPH)THEN
C     **********************************
C     *** FRAMES FOR NOT TPH'S CELLS ***
C     **********************************
         IFLAG=1    !this is just for now, but this flag is for frame's type
         DO I=1,N,2
            WRITE(12,22) IFLAG,XC(I),YC(I),ZC(I),XC(I+1),YC(I+1),ZC(I+1)
         ENDDO
 22      FORMAT(I1,2X,3F16.10,2X,3F16.10)
      ELSE
C     ******************************
C     *** FRAMES FOR TPH'S CELLS ***
C     ******************************
         ADIS=DIAGONAL(DVEC,1,2)/3
         CDIS=VECSIZE3(DVEC,3)
c         print *,'DVEC1=',(dvec(1,i),i=1,3)
c         print *,'DVEC2=',(dvec(2,i),i=1,3)
c         print *,'DVEC3=',(dvec(3,i),i=1,3)
c         print *,'CDIS=',CDIS
         DO I=1,NTPH-1
            DO J=I+1,NTPH
               DIS=DIST(XTPH(I),YTPH(I),ZTPH(I),XTPH(J),YTPH(J),ZTPH(J))
c               print *,'DIS=',DIS
               IF(EQUAL(DIS,ADIS)) THEN
                  IF(TPHFLAG(I).EQ.1.AND.TPHFLAG(J).EQ.1)
     *                 WRITE(12,23) '1',XTPH(I),YTPH(I),ZTPH(I),
     *                 XTPH(J),YTPH(J),ZTPH(J)
                  IF( (TPHFLAG(I).EQ.1.AND.TPHFLAG(J).EQ.2) .OR.
     *                 (TPHFLAG(I).EQ.2.AND.TPHFLAG(J).EQ.1) )
     *                 WRITE(12,23) '2',XTPH(I),YTPH(I),ZTPH(I),
     *                 XTPH(J),YTPH(J),ZTPH(J)
               ELSEIF(EQUAL(DIS,CDIS)) THEN
                  IF(TPHFLAG(I).EQ.1.AND.TPHFLAG(J).EQ.1)
     *                 WRITE(12,23) '1',XTPH(I),YTPH(I),ZTPH(I),
     *                 XTPH(J),YTPH(J),ZTPH(J)
 23               FORMAT(A1,2X,3F16.10,2X,3F16.10)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      RETURN
      END



C     *** THIS SUBROUTINE MAKE ONE FRAME ***
      SUBROUTINE MAKEFRAME(NATOM,N,DVEC,XC,YC,ZC,DI,DJ,DK,DI1,DJ1,DK1)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
c     **********************************************************************
c     in order to allow compiling on non-g77 compilers; replace NATOM arrays
c     with INATOM parameter
c     **********************************************************************
      REAL*8 XC(6*INATOM),YC(6*INATOM),ZC(6*INATOM),DVEC(3,3)
      
c      print *,'MAKEFRAME> n=',n

      N=N+1
      XC(N)=DI*DVEC(1,1)+DJ*DVEC(2,1)+DK*DVEC(3,1)
      YC(N)=DI*DVEC(1,2)+DJ*DVEC(2,2)+DK*DVEC(3,2)
      ZC(N)=DI*DVEC(1,3)+DJ*DVEC(2,3)+DK*DVEC(3,3)
      N=N+1
      XC(N)=DI1*DVEC(1,1)+DJ1*DVEC(2,1)+DK1*DVEC(3,1)
      YC(N)=DI1*DVEC(1,2)+DJ1*DVEC(2,2)+DK1*DVEC(3,2)
      ZC(N)=DI1*DVEC(1,3)+DJ1*DVEC(2,3)+DK1*DVEC(3,3)
      
      RETURN
      END


      
      REAL*8 FUNCTION DIAGONAL(DVEC,I,J)
      REAL*8 DVEC(3,3),NL,AI,AJ,PI,VECSIZE3,DSQRT
      PARAMETER (PI=3.14159265358979323844d0)

      NL=0.0d0
      AI=VECSIZE3(DVEC,I)
      AJ=VECSIZE3(DVEC,J)
      DIAGONAL=DSQRT(AI*AI+AJ*AJ-2*AI*AJ*DCOS(2*PI/3))
      RETURN
      END
