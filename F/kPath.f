      program kpath
c     ************************************
c     Usage: kpath <XSFfile> <k-list-file>
c     ************************************
      
c     k-list-dat should look like:
c     ----------------------------
c     KP, NK, IDV, Emin, Emax ... # of special k-points, # of-points, 
c                                 devider (IDV), Emin, Emax
c     IX IY IZ  LABEL ... ix, iy, iz -> k-coordinate: (ix,iy,iz)/idv
c                         REPEAT KP times
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      CHARACTER*128 file1,file2
      INTEGER icoor(3,NAC), nkp_seg(NAC), NAT(NAC)
      REAL*8 idvec(3,3), dis(NAC), frcoor(3,NAC), coor(3,NAC)
      REAL*8 prim(3), conv(3), inv(3,3)
      REAL*8 ivm(3,3)
      REAL*8 dis_sum, rpv(3,3), rcv(3,3), mat(3,3), dvec(3,3)
      REAL*8 X(NAC), Y(NAC), Z(NAC)
      CHARACTER*10 label(NAC)
      CHARACTER*4 input_vec
      LOGICAL ldiv
      COMMON/IVEC/ IDVEC
      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR

      PARAMETER ( TOLMIN=1.0d-8 )

      include 'get.inc'

      if ( iargc() .lt. 2 ) then
         write(*,*) 'Usage: kpath <XSFfile> <k-list-file>'
         STOP
      endif

c     read XSF file
      Call GetArg(1,FILE1)
      len=index(file1,' ')-1
      Open(11,FILE=FILE1(1:len),STATUS='OLD')
      
c     reciprocal vectors are in form:
c     ---------------------------------
c             | x1 y1 z1 |   | 11 21 31 |
c     IDVEC = | x2 y2 z2 | = | 12 22 32 |
c             | x3 y3 z3 |   | 13 23 33 |
      
c     calculate the distance between k-coordinates and estimate number 
c     of points for each segment

c     read k-list-file

      Call GetArg(2,FILE2)
      len=index(file2,' ')-1
      Open(21,FILE=FILE2(1:len),STATUS='OLD')      
      read(21,*,END=99) KP, NK, IDV, Emin, Emax
      read(21,'(A4)') input_vec
      if (kp.le.1) then
         write(*,*) 'less then two special k-points selected'
         stop
      endif

c     supportInfo-1.kpath file
      open(33,file='supportInfo-1.kpath',status='unknown')
c     supportInfo-2.kpath file
      open(34,file='supportInfo-2.kpath',status='unknown')

      Call ReadF1('DIM-GROUP',IGROUP,IDIM)
c     read both sets of reciprocal vectors

      if ( input_vec(1:4) .eq. 'prim' ) then
         Call ReadF1('RECIP-CONVVEC', IDUM1, IDUM2)
         do i=1,3
            do j=1,3
               rcv(j,i) = idvec(j,i)
            enddo
         enddo
         Call ReadF1('RECIP-PRIMVEC', IDUM1, IDUM2)
         do i=1,3
            do j=1,3
               rpv(j,i) = idvec(j,i)
            enddo
         enddo
      else
         Call ReadF1('RECIP-PRIMVEC', IDUM1, IDUM2)
         do i=1,3
            do j=1,3
               rpv(j,i) = idvec(j,i)
            enddo
         enddo
         Call ReadF1('RECIP-CONVVEC', IDUM1, IDUM2)
         do i=1,3
            do j=1,3
               rcv(j,i) = idvec(j,i)
            enddo
         enddo
      endif

c     **************************************************************
c     get the transpose of idvec; this is the matrix for getting the
c     Cartesian coordinates.
      call MatTranspose(idvec, ivm, 3, 3)      
c     **************************************************************

      write(33,'(''   Selected k-points in CARTESIAN coordinates:'')')

      dis_sum=0.0d0
      do i=1,KP
         read(21,*,END=99) (icoor(j,i),j=1,3), label(i)
         do j=1,3
            frcoor(j,i) = dble(icoor(j,i)) / dble(idv)
         enddo
c     transform coor from fractional to cartesian COORDINATES
         coor(1,i) = frcoor(1,i)*ivm(1,1) + 
     1        frcoor(2,i)*ivm(2,1) + frcoor(3,i)*ivm(3,1)
         
         coor(2,i) = frcoor(1,i)*ivm(1,2) + 
     1        frcoor(2,i)*ivm(2,2) + frcoor(3,i)*ivm(3,2)
         
         coor(3,i) = frcoor(1,i)*ivm(1,3) + 
     1        frcoor(2,i)*ivm(2,3) + frcoor(3,i)*ivm(3,3)
         
         write(33,'(3x,3(f10.5,1x),3x,a5)')
     $        coor(1,i), coor(2,i), coor(3,i), label(i)
         
         if (i .gt. 1) then
            dis(i-1) = Dist( coor(1,i), coor(2,i), coor(3,i),
     1           coor(1,i-1), coor(2,i-1), coor(3,i-1) )
            dis_sum = dis_sum + dis(i-1)
         endif
      enddo

      call ZeroMat(mat,3,3)

c     if needed convert frcoor to other basis (prim/conv)

      if (IConvCoor(igroup,input_vec,rcv) .eq. IGET_PRIM) then
         Call ReadF1('PRIMVEC', idum1, idum2)
c     bug resolved by Florent Boucher
         call MatMult(rcv, 3, 3, dvec, 3, mat)
         call GetIntCoor(frcoor, icoor, mat, idv, kp)
      else if (IConvCoor(igroup,input_vec,rcv) .eq. IGET_CONV) then
c     convert from primitive to conventional frcoor
         Call ReadF1('CONVVEC', idum1, idum2)
c     bug resolved by Florent Boucher
         call MatMult(rpv, 3, 3, dvec, 3, mat)
         call GetIntCoor(frcoor, icoor, mat, idv, kp)
      endif

c     now write the K-list

      do i=1,kp-1
         nkp_seg(i) = nint( dble(nk) * dis(i) / dis_sum ) 
         if (nkp_seg(i) .eq. 0) then
            WRITE(*,*) 
     1           'ERROR: IDV*NK to big! Decrease number of k-points'
            STOP
         endif
cfb
         new_idv    = nkp_seg(i)*idv
         mult       = nkp_seg(i)
cfb	 
         ix = icoor(1,i)*mult
         iy = icoor(2,i)*mult
         iz = icoor(3,i)*mult
         if (i.eq.1) then
            write(*,'(A10,4I5,F5.1,2F5.2,2x,A30)') 
     1           label(i), ix, iy, iz, new_idv, 2.0, Emin, Emax,
     $           'k-list generated by XCrySDen'
         else
            write(*,'(A10,4I5,F5.1)') label(i), ix, iy, iz, new_idv, 2.0
         endif

         do j=1,nkp_seg(i)-1
            ix = icoor(1,i)*mult + 
     1           j*mult*(icoor(1,i+1) - icoor(1,i))/nkp_seg(i)
            iy = icoor(2,i)*mult +                   
     1           j*mult*(icoor(2,i+1) - icoor(2,i))/nkp_seg(i)
            iz = icoor(3,i)*mult +                   
     1           j*mult*(icoor(3,i+1) - icoor(3,i))/nkp_seg(i)
            write(*,'(10X,4I5,F5.1)') ix, iy, iz, new_idv, 2.0
         enddo

c     supporting info         
         write(34,'(''K-point '',a10,'':'', 3i5,2x,''/'',1x,i5,
     $        1x,''-->'',1x,3f8.4)')
     $        label(i),
     $        icoor(1,i)*mult, icoor(2,i)*mult, icoor(3,i)*mult,
     $        new_idv,
     $        dble(icoor(1,i)*mult)/dble(new_idv),
     $        dble(icoor(2,i)*mult)/dble(new_idv),
     $        dble(icoor(3,i)*mult)/dble(new_idv)
c     /
      enddo
      ix = icoor(1,kp)*mult
      iy = icoor(2,kp)*mult
      iz = icoor(3,kp)*mult
      write(*,'(A10,4I5,F5.1)') label(kp), ix, iy, iz, new_idv, 2.0
      write(*,'(A3)') 'END'
c     supporting info         
      write(34,'(''K-point '',a10,'':'', 3i5,2x,''/'',1x,i5,
     $     1x,''-->'',1x,3f8.4)')
     $     label(kp), ix, iy, iz, new_idv,
     $     dble(ix)/dble(new_idv),
     $     dble(iy)/dble(new_idv),
     $     dble(iz)/dble(new_idv)
c     /
      GOTO 1000

 99   continue
      STOP 'unexpected end of k-list-file'

 1000 continue

c     
c     POSTSCRIPT: write some supporting information
c     

c     RECIPROCAL-PRIMITIVE
      write(33,'('' '')')
      write(33,'(''   Selected k-points in crystal coordinates:'')')
      write(33,
     $     '(''   (with respect to RECIPROCAL-PRIMITIVE vectors)'')')
      
      call Invert3x3(rpv,inv)
      do i=1,KP
         do j=1,3
c            prim(j) =
c     $           coor(1,i)*inv(1,j) +
c     $           coor(2,i)*inv(2,j) +
c     $           coor(3,i)*inv(3,j)
            prim(j) =
     $           coor(1,i)*inv(j,1) +
     $           coor(2,i)*inv(j,2) +
     $           coor(3,i)*inv(j,3)
         enddo
         write(33,'(3x,3(f10.5,1x),3x,a5)')
     $        prim(1), prim(2), prim(3), label(i)
      enddo

c     RECIPROCAL-CONVENTIONAL
      write(33,'('' '')')
      write(33,'(''   Selected k-points in crystal coordinates:'')')
      write(33,
     $     '(''   (with respect to RECIPROCAL-CONVENTIONAL vectors)'')')
      
      call Invert3x3(rcv,inv)
      do i=1,KP
         do j=1,3
c            conv(j) =
c     $           coor(1,i)*inv(1,j) +
c     $           coor(2,i)*inv(2,j) +
c     $           coor(3,i)*inv(3,j)
            conv(j) =
     $           coor(1,i)*inv(j,1) +
     $           coor(2,i)*inv(j,2) +
     $           coor(3,i)*inv(j,3)
         enddo
         write(33,'(3x,3(f10.5,1x),3x,a5)')
     $        conv(1), conv(2), conv(3), label(i)
      enddo

      END

