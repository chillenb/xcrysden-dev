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
         
c         print *,'i, coor::', i,  coor(1,i), coor(2,i), coor(3,i)
         if (i .gt. 1) then
            dis(i-1) = Dist( coor(1,i), coor(2,i), coor(3,i),
     1           coor(1,i-1), coor(2,i-1), coor(3,i-1) )
c            print *, 'dis::', i-1, dis(i-1)
            dis_sum = dis_sum + dis(i-1)
         endif
      enddo

      call ZeroMat(mat,3,3)
c     if needed convert frcoor to other basis (prim/conv)
      if (IConvCoor(igroup,input_vec,rcv) .eq. IGET_PRIM) then
         Call ReadF1('PRIMVEC', idum1, idum2)
c         call MatMult(dvec, 3, 3, rcv, 3, mat) 
c     bug resolved thanks to Florent Boucher
         call MatMult(rcv, 3, 3, dvec, 3, mat)
         call GetIntCoor(frcoor, icoor, mat, idv, kp)
      else if (IConvCoor(igroup,input_vec,rcv) .eq. IGET_CONV) then
c     convert from primitive to conventional frcoor
         Call ReadF1('CONVVEC', idum1, idum2)
c         call MatMult(dvec, 3, 3, rpv, 3, mat)
c     bug resolved thanks to Florent Boucher
         call MatMult(rpv, 3, 3, dvec, 3, mat)
         call GetIntCoor(frcoor, icoor, mat, idv, kp)
      endif

c     najmanjsi skupni delitelj bo PRODUCT(nkp_seg(i),i=1,kp-1)*idv,
c     zato morajo biti nkp_seg zelo lepo deljivi, mogoce z deset ali 
c     kaj takega. KP, ki jo je podal uporabnik naj bo samo za informacijo!!!
      k_aver = nint( dble(nk) / dble(kp-1))
      nn = 10
      if (k_aver .lt. 40) nn = 5
      if (k_aver .lt. 30) nn = 4
      if (k_aver .lt. 20) nn = 3
      if (k_aver .lt. 10) nn = 2

 100  continue
c      print *,'nn::', nn

c     estimate number of k-points for each segment
      mult=1
      nkp_max = 1
      do i=1,KP-1
c     here nk_seg(i) is nn times smaler then it should be
         nkp_seg(i) = nint( dble(nk) * dis(i) / (dis_sum * nn) ) 
c         print *, 'ratio:', dis(i) / dis_sum,
c     1        ' nkp_seg(',i,')', nkp_seg(i)
         if (nkp_seg(i) .gt. nkp_max) nkp_max = nkp_seg(i)
         mult = mult * nkp_seg(i)
      enddo      

c      print *, 'mult, nkp_max::', mult, nkp_max
      if (kp.gt.3) then
         do i=nkp_max,2,-1
            ldiv=.true.
            do j=1,kp-1
               if ( Mod(nkp_seg(j),i) .ne. 0 ) ldiv=.false.
c               print *, 'i, j, Mod(nkp_seg(j),i)', i, j,
c     $              mod(nkp_seg(j),i)
            enddo
            if (ldiv) then
c               print *, 'mult & i:', i
               mult = mult / i
            endif
         enddo
      endif
      mult=mult*nn

      if ( idv*mult .gt. 9999) then
c     idv_common*mult is to large to fit to i4 format 
c     (i4 because some (ix,iy,iz) might be negative, so we need one 
c     place for minus); ENLARGE the nn
         if(nn.ge.10) nn=2*nn
         if(nn.eq.5)  nn=10
         if(nn.eq.4)  nn=5
         if(nn.eq.3)  nn=4
         if(nn.eq.2)  nn=3
         goto 100
      endif
c      print *, 'mult #2', mult

c     *************************************************
c     *** get the right fractional coordinates for WIEN
c     *************************************************
      new_idv = mult*idv 
c     now write the K-list
      do i=1,kp-1
         nkp_nn     = nkp_seg(i)*nn
         if (nkp_seg(i) .eq. 0) then
            WRITE(*,*) 
     1           'ERROR: IDV*NK to big! Decrease number of k-points'
            STOP
         endif
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
c         print *,'i,nkp_seg::', i, nkp_seg(i)
         do j=1,nkp_nn-1
            ix = icoor(1,i)*mult + 
     1           j*mult*(icoor(1,i+1) - icoor(1,i))/nkp_nn
            iy = icoor(2,i)*mult +                   
     1           j*mult*(icoor(2,i+1) - icoor(2,i))/nkp_nn
            iz = icoor(3,i)*mult +                   
     1           j*mult*(icoor(3,i+1) - icoor(3,i))/nkp_nn
            write(*,'(10X,4I5,F5.1)') ix, iy, iz, new_idv, 2.0
         enddo
      enddo
      ix = icoor(1,kp)*mult
      iy = icoor(2,kp)*mult
      iz = icoor(3,kp)*mult
      write(*,'(A10,4I5,F5.1)') label(kp), ix, iy, iz, new_idv, 2.0
      write(*,'(A3)') 'END'
      GOTO 1000
            
 99   continue
      STOP 'unexpected end of k-list-file'

 1000 continue
      END
