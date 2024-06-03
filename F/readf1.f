      SUBROUTINE READF1(ARG,IGROUP,IDIM)
C     THIS ROUTINE READS "GENGEOM" FILE
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      CHARACTER*(*) ARG
      CHARACTER*20 ARGM(12),CARG
      REAL*8 X(NAC),Y(NAC),Z(NAC),DVEC(3,3),IDVEC(3,3),
     1     TMPVEC(3,3),XC(NAC,4),YC(NAC,4),ZC(NAC,4),fv(NAC,3)
      INTEGER NAT(NAC) !I THINK NAC ATOMS PER CELL IS MORE THAN ENOUGH
      CHARACTER*256 LINE
      LOGICAL notfound, Lprimvec

      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR
      COMMON/FV/ FV,NFV
      COMMON/IVEC/ IDVEC
      COMMON/BOHRR/ BOHR,IMODE3
      COMMON/KEYWORD/ notfound, Lprimvec

      include 'mode.inc'
      PARAMETER (NARGM=12)

      DATA ARGM/'DIM-GROUP','MOLECULE','POLYMER','SLAB','CRYSTAL',
     $     'PRIMVEC','CONVVEC','PRIMCOORD',
     *     'CONVCOORD','RECIP-PRIMVEC','RECIP-CONVVEC','SYMMOP'/

      notfound = .false.
      IND=999

C     *********
      REWIND 11
C     *********

      LARG=INDEX(ARG,' ')-1
      IF(LARG.EQ.-1)LARG=LEN(ARG)

      DO I=1,NARGM
         CARG=ARGM(I)
         LCARG=INDEX(CARG,' ')-1
c     print *,'carg=',carg(1:5),'|'
         IF(CARG(1:LCARG).EQ.ARG(1:LARG)) THEN
            IND=I
            GOTO 1
         ENDIF
      ENDDO
c     **********************************************************
      WRITE(6,*) 'READF1> KEYWORD ', ARG, ' NOT KNOWN!!! <STOP>'
      STOP
c     **********************************************************


c     is desired keyword present in the file
 1    CONTINUE
   
c     DEBUG_BEGIN
c      READ(UNIT=11,FMT=*,END=99) CARG
      READ(UNIT=11,FMT='(a20)',END=99) CARG
      i = string_length(carg)
c     DEBUG_END

      LCARG=INDEX(CARG,' ')-1
c      print *,'|',carg(1:lcarg),'|'
      IF(CARG(1:LCARG).NE.ARG(1:LARG)) GOTO 1

      GOTO(5,6,7,8,9,10,15,20,25,35,45,999) IND

 99   continue
c     ****************************
c     this is nasty, but necessary
c     ****************************
      if ( arg(1:7) .eq. 'CONVVEC' ) then
         WRITE(12,'(1x,a)') 'CONVVEC'
         WRITE(12,'(3(1x,f15.10)/3(1x,f15.10)/3(1x,f15.10))') 
     1        ((DVEC(J,I),I=1,3),J=1,3)
      endif
      notfound=.true.
 999  RETURN

C     READ DIMENSION & GROUP NUMBER
 5    CONTINUE
      READ(11,*) IDIM,IGROUP
      WRITE(12,'(1x,a)') 'DIM-GROUP'
      WRITE(12,*) IDIM,IGROUP
      RETURN

C     keyword MOLECULE
 6    CONTINUE
      IDIM=0
      IGROUP=1
      RETURN
      
C     keyword POLYMER
 7    CONTINUE
c      WRITE(12,'(1x,a)') 'POLYMER'
      IDIM=1
      IGROUP=1
      WRITE(12,'(1x,a)') 'DIM-GROUP'
      WRITE(12,*) IDIM,IGROUP
      RETURN

C     keyword SLAB
 8    CONTINUE
c      WRITE(12,'(1x,a)') 'SLAB'
      IDIM=2
      IGROUP=1
      WRITE(12,'(1x,a)') 'DIM-GROUP'
      WRITE(12,*) IDIM,IGROUP
      RETURN

C     keyword CRYSTAL
 9    CONTINUE
c      WRITE(12,'(1x,a)') 'CRYSTAL'
      IDIM=3
      IGROUP=1
      WRITE(12,'(1x,a)') 'DIM-GROUP'
      WRITE(12,*) IDIM,IGROUP
      RETURN

 10   CONTINUE
      WRITE(12,'(1x,a)') 'PRIMVEC'
      GOTO 16
 15   CONTINUE
      WRITE(12,'(1x,a)') 'CONVVEC'
 16   CONTINUE
C     READ CELL VECTORS
      READ(11,*) ((DVEC(J,I),I=1,3),J=1,3)
      IF(IMODE3.EQ.M3_BOHR)THEN  !convert from borh to angstroms
         DO I=1,3
            DO J=1,3
               DVEC(J,I) = BOHR * DVEC(J,I)
            ENDDO
         ENDDO
      ENDIF
c     make the VECTORS consistent with dimenison
      call DIMIFY_VEC(DVEC,IDIM)

      WRITE(12,'(3(1x,f15.10)/3(1x,f15.10)/3(1x,f15.10))') 
     1     ((DVEC(J,I),I=1,3),J=1,3)
C     *** sometimes we will need inverse matrix of DVEC ***
      call InvertSqMat123(dvec,idvec,3)
      RETURN

 20   CONTINUE
      WRITE(12,'(1x,a)') 'PRIMCOORD'
      READ(11,*) NATR,NCELL
      WRITE(12,*) NATR,NCELL
c      print *,'natr=',natr,';  ncell=',ncell
      DO I=1,NATR
C     ALL NAT(I,*) ARE THE SAME -> SO NAT(I)            
 777     continue
         read(11,'(a)') line
c     print *,'line=',line,'|'
         nf=iCountFields(line)
         if (nf.lt.1) goto 777
         NFV=nf-4
         read(line,*) NAT(I),X(I),Y(I),Z(I), (fv(i,j), j=1,NFV)

         IF(IMODE3.EQ.M3_BOHR) THEN
            X(I) = X(I) * BOHR
            Y(I) = Y(I) * BOHR
            Z(I) = Z(I) * BOHR
         ENDIF
         WRITE(12,'(i3,3x,6(f15.10,2x))') NAT(I),X(I),Y(I),Z(I),
     $        (fv(i,j), j=1,NFV)
      ENDDO
      RETURN

 25   CONTINUE

c     ********
c     CONVCOOR
c     ********
      if (.not.Lprimvec) WRITE(12,'(1x,a)') 'PRIMCOORD'
C     READ COORDINATES
      READ(11,*) NATR,NCELL
      if (.not.Lprimvec) WRITE(12,*) NATR, ' 1'
c      print *,'natr=',natr,';  ncell=',ncell
      DO I=1,NATR
         DO J=1,NCELL
C     ALL NAT(I,*) ARE THE SAME -> SO NAT(I)            
            READ(11,*) NAT(I),XC(I,J),YC(I,J),ZC(I,J) 
            IF(IMODE3.EQ.M3_BOHR) THEN
               XC(I,J) = XC(I,J) * BOHR
               YC(I,J) = YC(I,J) * BOHR
               ZC(I,J) = ZC(I,J) * BOHR
            ENDIF
            if (.not.Lprimvec) WRITE(12,'(i3,3x,3(f15.10,2x))')
     $           NAT(I),XC(I,J),YC(I,J),ZC(I,J)
         ENDDO
         if(.not.Lprimvec) then
            x(i) = xc(i,1)
            y(i) = yc(i,1)
            z(i) = zc(i,1)
         endif
      ENDDO
      RETURN

 35   CONTINUE
      WRITE(12,'(1x,a)') 'RECIP-PRIMVEC'
      GOTO 46
 45   CONTINUE
      WRITE(12,'(1x,a)') 'RECIP-CONVVEC'
 46   CONTINUE
C     READ RECIPROCAL VECTORS: they are written in normal form, 
C     but should be read in TRANSPOSE form; take care
      READ(11,*) ((IDVEC(J,I),J=1,3),I=1,3)
      IF(IMODE3.EQ.M3_BOHR)THEN !convert from borh to angstroms
         DO I=1,3
            DO J=1,3
               IDVEC(J,I) = BOHR * IDVEC(J,I)
            ENDDO
         ENDDO
      ENDIF
      WRITE(12,'(3(1x,f15.10)/3(1x,f15.10)/3(1x,f15.10))')
     $     ((IDVEC(J,I),J=1,3),I=1,3)
      RETURN
C     *** sometimes we will need inverse matrix of DVEC ***

      END
