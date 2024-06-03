c     ------------------------------------------------------------------------
      program kpath
c     Usage: kpath <XSFfile> <k-list-file>
c     ------------------------------------------------------------------------
      
c     k-list-dat should look like:
c     ----------------------------
c     KP, NK            ... No. of special k-points, No. interp. k-points
c     crystal_coor_type ... prim or conv
c     KX KY KZ label    ... special k-point coordinates, label

      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'

      CHARACTER*128 file1,file2
      INTEGER nkp_seg(NAC), NAT(NAC), nk, kp, nkp_sum, i, j
      REAL*8 idvec(3,3), dis(NAC), frcoor(3,NAC), coor(3,NAC), fr(3,NAC)
      REAL*8 prim(3), conv(3), inv(3,3)
      REAL*8 ivm(3,3)
      REAL*8 dis_sum, rpv(3,3), rcv(3,3), m(3,3), dvec(3,3)
      REAL*8 X(NAC), Y(NAC), Z(NAC), fx, fy, fz, f
      CHARACTER*10 label(NAC)
      CHARACTER*4 input_vec
      LOGICAL ldiv
      COMMON/IVEC/ IDVEC
      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR

      PARAMETER ( TOLMIN=1.0d-8 )

c     reciprocal vectors are in form:
c     ---------------------------------
c             | x1 y1 z1 |   | 11 21 31 |
c     IDVEC = | x2 y2 z2 | = | 12 22 32 |
c             | x3 y3 z3 |   | 13 23 33 |


      if ( iargc() .lt. 2 ) then
         write(*,*) 'Usage: kpath <XSFfile> <k-list-file>'
         STOP
      endif

c     read XSF file

      call GetArg(1,FILE1)
      len=index(file1,' ')-1
      open(11,FILE=FILE1(1:len),STATUS='OLD')
            
c     read k-list-file

      call GetArg(2,FILE2)
      len=index(file2,' ')-1
      open(21,FILE=FILE2(1:len),STATUS='OLD')      
      read(21,*,END=99) KP, NK
      read(21,'(A4)') input_vec

      if (kp.le.1) then
         write(*,*) 'less then two special k-points selected'
         stop
      endif

c     supportInfo-1.kpath and supportInfo-2.kpath files

      open(33,file='supportInfo-1.kpath',status='unknown')
c      open(34,file='supportInfo-2.kpath',status='unknown')

      call ReadF1('DIM-GROUP',IGROUP,IDIM)

c     read both sets of reciprocal vectors

      if ( input_vec(1:4) .eq. 'prim' ) then
         call ReadF1('RECIP-CONVVEC', IDUM1, IDUM2)
         call MatCopy(idvec,rcv,3,3)
         
         call ReadF1('RECIP-PRIMVEC', IDUM1, IDUM2)
         call MatCopy(idvec,rpv,3,3)
      else
         call ReadF1('RECIP-PRIMVEC', IDUM1, IDUM2)
         call MatCopy(idvec,rpv,3,3)

         call ReadF1('RECIP-CONVVEC', IDUM1, IDUM2)
         call MatCopy(idvec,rcv,3,3)
      endif

c     **************************************************************
c     get the transpose of idvec; this is the matrix for getting the
c     Cartesian coordinates.
      call MatTranspose(idvec, ivm, 3, 3)      
c     **************************************************************

      write(33,'(''   Selected k-points in CARTESIAN coordinates:'')')

c     calculate the distance between k-coordinates and estimate number 
c     of points for each segment
      
      dis_sum=0.0d0
      do i=1,KP
         read(21,*,END=99) (frcoor(j,i),j=1,3), label(i)

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
            dis(i-1) = dist( coor(1,i), coor(2,i), coor(3,i),
     1           coor(1,i-1), coor(2,i-1), coor(3,i-1) )
            dis_sum = dis_sum + dis(i-1)
         endif
      enddo

      call ZeroMat(m,3,3)

      if ( input_vec(1:4) .eq. 'conv' ) then
c     convert frcoor to primitive basis
         call Invert3x3(rpv,inv)
         do i=1,kp
            do j=1,3
               frcoor(j,i) =
     $              coor(1,i)*inv(j,1) +
     $              coor(2,i)*inv(j,2) +
     $              coor(3,i)*inv(j,3)
            enddo
         enddo
      endif
      
c     now write the K-list

      nkp_sum=0
      do i=1,kp-1
         nkp_seg(i) = nint( dble(nk-1) * dis(i) / dis_sum ) 
         nkp_sum = nkp_sum + nkp_seg(i)
         if (nkp_seg(i) .eq. 0) then
            write(*,*) 
     1           'ERROR: NK too small!!! Increase number of k-points'
            STOP
         endif
      enddo

      write(*,*) 'K_POINTS crystal'
      write(*,*) nkp_sum+1

      do i=1,kp-1
         do j=1,nkp_seg(i)
            f = dble(j-1)/dble(nkp_seg(i)) 
            fx = frcoor(1,i) + f*(frcoor(1,i+1) - frcoor(1,i))
            fy = frcoor(2,i) + f*(frcoor(2,i+1) - frcoor(2,i))
            fz = frcoor(3,i) + f*(frcoor(3,i+1) - frcoor(3,i))
            write(*, '(5x,3(f15.10,1x),3x,''1.0'')') fx, fy, fz
c     supporting info         
c            write(34,'(5x,3(f15.10,1x),3x,''1.0'')') fx, fy, fz
         enddo
      enddo
      fx = frcoor(1,kp)
      fy = frcoor(2,kp)
      fz = frcoor(3,kp)
      write(*, '(5x,3(f15.10,1x),3x,''1.0'')') fx, fy, fz
c     supporting info         
c      write(34,'(5x,3(f15.10,1x),3x,''1.0'')') fx, fy, fz
      
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

