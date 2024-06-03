      program SaveSFile
c     **********************************************************************
c     program CONVERT XSF format to WIEN STRUCT FILE format
c     **********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      INTEGER NAT(NAC)
      INTEGER C2I
      REAL*8 X(NAC),Y(NAC),Z(NAC),DVEC(3,3),IDVEC(3,3),TMPVEC(3,3)
      REAL*8 iRdvec(3,3)
      REAL*8 A(NAC),B(NAC),C(NAC)
      REAL*8 SYMMOP(1,3,3), SYMMTR(1,3)
      CHARACTER*80 FILE1, TITLE, MODE
      CHARACTER*4 GTYPE
      LOGICAL l_inv_center, l_inv
      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR
      COMMON/IVEC/ IDVEC
      COMMON/BOHRR/ BOHR,IMODE3 !this is not used, but is needed because 
                                !this program uses readf1 routine
      PARAMETER (RAD2DEG = 57.295779513084, TOLMIN=1.d-7)

      DATA SYMMOP/
     $     1.0d0, 0.0d0, 0.0d0,
     $     0.0d0, 1.0d0, 0.0d0,
     $     0.0d0, 0.0d0, 1.0d0/
      DATA SYMMTR/
     $     0.0d0, 0.0d0, 0.0d0/

      narg=IArgc()
      if (narg.ne.1) then         
         stop 'Usage: savestruct file1'
      endif

      Call GetArg(1,FILE1)
      Open(unit=11,file=file1,status='old')

      Call ReadF1('DIM-GROUP', igroup, idim)
      if (idim.ne.3) stop 'SaveStruct:: so far I am converting to
     $     WIEN97 STRUCT FILE just from crystals (dim=3) !!!'

      if (igroup .eq. 1) GTYPE='P   '
      if (igroup .eq. 2) GTYPE='CYZ ' ! A type
      if (igroup .eq. 3) GTYPE='CXZ ' ! B type
      if (igroup .eq. 4) GTYPE='CXY ' ! C type
      if (igroup .eq. 5) GTYPE='F   ' 
      if (igroup .eq. 6) GTYPE='B   '
      if (igroup .eq. 7) GTYPE='R   '
      if (igroup .eq. 8) GTYPE='H   '
      if (igroup .eq. 9) GTYPE='H   '
      if (igroup .gt. 10) GTYPE='P   ' !just in case

c     this should works for H,P; but for F, B, C R, CONVVEC/PRIMCOORD 
c     should probably be used  !!!!!!!!!!!!!1
      Call ReadF1('PRIMVEC', 0, idim) !dvec & idvec are assigned
      Call ReadF1('PRIMCOORD', 0, idim)
            
c     get the inverse matrix
      call Invert3x3(dvec,iRdvec)
      
      if ((igroup .gt. 1) .and. (igroup .lt. 8)) then
         Call ReadF1('CONVVEC', 0, idim) !dvec & idvec are assigned
         if (igroup .eq. 7) then
c     for rhombohedral hexagonal vectors and rhobohedral coordinates 
c     must be specified
            do i=1,3
               do j=1,3
                  idvec(j,i)=iRdvec(j,i)
               enddo
            enddo
         endif
      endif

      Close(11)
      
      Call SaveStructFile(gtype, dvec, idvec, natr, nat, x, y, z)         
      END


      subroutine SaveStructFile(gtype, dvec, idvec, natr, nat, x, y, z)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      INTEGER NAT(NAC), INN(NAC)
      INTEGER MAXNAT            ! maximum atomic number in structure
      INTEGER NPT, ISPLIT       ! used for Wien97
      REAL*8 X(NAC),Y(NAC),Z(NAC),DVEC(3,3),IDVEC(3,3),TMPVEC(3,3)
      REAL*8 A(NAC),B(NAC),C(NAC), MIN_NNDIST(NAC), RMT(NAC), R0(NAC)
      REAL*8 SYMMOP(1,3,3), SYMMTR(1,3)
      REAL*8 V1(3), V2(3), V3(3)
      CHARACTER*10 ATOM(100), ATOM_UPPER(100)
      CHARACTER*4 gtype
      REAL*8 FRMT(NAC)      

      PARAMETER (RAD2DEG = 57.295779513084, TOLMIN=1.d-7)

      include 'atoms.inc'

C     CONVERSION FROM BOHRS TO ANGS.
      BOHR=0.529177d0

c     ****************************************************************
c     now get the fractional coordinates of atoms and 
c     translate them if necessary ( in WIEN97 fractional coorinates 
c     are in range [0,1) )
      maxnat = 1
      do i=1,natr
         nat(i) = nat(i) - int( nat(i) / 100 ) * 100 ! nat from 0 to 99
         if(nat(i).gt.maxnat) maxnat = nat(i)
         a(i) = x(i)*idvec(1,1) + y(i)*idvec(2,1) + z(i)*idvec(3,1)
         b(i) = x(i)*idvec(1,2) + y(i)*idvec(2,2) + z(i)*idvec(3,2)
         c(i) = x(i)*idvec(1,3) + y(i)*idvec(2,3) + z(i)*idvec(3,3)
         if(a(i).lt.0.0d0) a(i) = a(i) + 1.0d0
         if(b(i).lt.0.0d0) b(i) = b(i) + 1.0d0
         if(c(i).lt.0.0d0) c(i) = c(i) + 1.0d0
c     Write(6,21) nat(i), a(i), b(i), c(i)
      enddo
      
c     get the nearest neighbor distances
      do i=1,natr
         min_nndist(i) = Dist(dvec(1,1), dvec(1,2), dvec(1,3),
     $        dvec(2,1), dvec(2,2), dvec(2,3)) !this always greater then nndist
         do j=1,natr
            if (i.ne.j) then
               a2 = a(j) + nint(a(i)-a(j))
               b2 = b(j) + nint(b(i)-b(j))
               c2 = c(j) + nint(c(i)-c(j))
               
               x1 = a(i)*dvec(1,1) + b(i)*dvec(2,1) + c(i)*dvec(3,1)
               y1 = a(i)*dvec(1,2) + b(i)*dvec(2,2) + c(i)*dvec(3,2)
               z1 = a(i)*dvec(1,3) + b(i)*dvec(2,3) + c(i)*dvec(3,3)
               
               x2 = a2*dvec(1,1) + b2*dvec(2,1) + c2*dvec(3,1)
               y2 = a2*dvec(1,2) + b2*dvec(2,2) + c2*dvec(3,2)
               z2 = a2*dvec(1,3) + b2*dvec(2,3) + c2*dvec(3,3)
               
               d = Dist(x1,y1,z1, x2,y2,z2)
               if (d.lt.min_nndist(i)) then
                  min_nndist(i) = d
                  inn(i) = j
               endif
            endif
         enddo
      enddo
      
c     *************************
c     get a,b,c,alpha,beta,gamma
      ap = VecSize3(dvec,1) / BOHR
      bp = VecSize3(dvec,2) / BOHR
      cp = VecSize3(dvec,3) / BOHR
      do i=1,3
         v1(i) = dvec(1,i)
         v2(i) = dvec(2,i)
         v3(i) = dvec(3,i)
      enddo
      alpha = dacos( ScalarProduct(v2, v3) /
     $     (VecSize(v2) * VecSize(v3)) ) * RAD2DEG
      beta  = dacos( ScalarProduct(v3,v1) /
     $     (VecSize(v3) * VecSize(v1)) ) * RAD2DEG
      gamma = dacos( ScalarProduct(v1,v2) /
     $     (VecSize(v1) * VecSize(v2)) ) * RAD2DEG
      
c     *******************************************
c     get some aproximation for muffin tin radius
c     *******************************************
c     let the smallest offset between RMT spheres be 0.1 a.u.
      OFFSET = 0.1d0/2.d0       ! each sphre gets half of offset
      
      do i=1,natr
         if (nat(i).eq.1) frmt(i) = 0.50d0
         if (nat(i).gt.1 .and. nat(i).lt.21)  frmt(i) = 1.00d0
         if (nat(i).ge.21 .and. nat(i).lt.39) frmt(i) = 1.10d0
         if (nat(i).ge.39 .and. nat(i).lt.57) frmt(i) = 1.15d0
         if (nat(i).ge.57) frmt(i) = 1.20d0
         
         r0(i) = 0.000005d0
         if (nat(i).le.71) r0(i) = 0.00001d0
         if (nat(i).le.36) r0(i) = 0.00005d0
         if (nat(i).le.18) r0(i) = 0.0001d0
      enddo
      
      do i=1,natr
         rmt(i) = frmt(i) * min_nndist(i) /
     $        ( BOHR * (frmt(i) + frmt(inn(i))) ) - OFFSET
      enddo
      isplit = 8
      npt    = 781
      
c     ****************************
c     now write WIEN97 struct file
c     ****************************
      write(*,'(a)')
     $     'Wien97 struct file generated by XCrySDen program'
      write(*,'(a4,a23,i3)') gtype, 'LATTICE.NONEQUIV. ATOMS',natr
      if (maxnat.lt.40) then
         write(*,'(a13,a4)') 'MODE OF CALC=', 'NREL'
      else
         write(*,'(a13,a4)') 'MODE OF CALC=', 'RELA'
      endif
      write(*,'(6f10.6)') ap, bp, cp, alpha, beta, gamma
      
      do i=1,natr
c     i should be positive for cubic symmetry and negative for non-cubic
c     *** ** * CHANGE THAT * ** ***
         write(*,'(a4,i4,a4,f10.8,a3,f10.8,a3,f10.8)')
     $        'ATOM', -i,': X=', a(i), ' Y=', b(i), ' Z=', c(i)
         write(*,'(a15,i2,a17,i2)') 'MULT=', 1, 'ISPLIT=', isplit
         write(*,'(a10,a5,i5,a5,f10.8,a5,f10.4,a5,f5.1)')
     $        atom(nat(i)),'NPT=',npt,'R0=',r0(i),
     $        'RMT=',rmt(i),'Z:', real(nat(i))
         write(*,'(a20,3f10.7)') 'LOCAL ROT MATRIX:   ', 1.0,0.0,0.0
         write(*,'(20x,3f10.7)') 0.0, 1.0, 0.0
         write(*,'(20x,3f10.7)') 0.0, 0.0, 1.0
      enddo
      write(*,'(i4,A)') 1,'      NUMBER OF SYMMETRY OPERATIONS'
      write(*,'(3i2,f10.7)') 1, 0, 0, 0.0
      write(*,'(3i2,f10.7)') 0, 1, 0, 0.0
      write(*,'(3i2,f10.7)') 0, 0, 1, 0.0
      write(*,'(i8)') 1
 21   FORMAT(I5,3F16.10)
      
      END
      

 
