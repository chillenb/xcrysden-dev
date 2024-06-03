      Integer Function IConvCoor(igroup, input_vec, rcv)
      CHARACTER*4 input_vec      
      REAL*8 rcv(3,3), scalar, tolmin
      include 'get.inc'
      tolmin=1.0d-6

C--------------------
C     GROUP NUMBERS::
C--------------------
C     N.  symbol
C     1.....P   
C     2.....A (=CXY)  
C     3.....B (=CXZ)  
C     4.....C (=CYZ)
C     5.....F   
C     6.....I   (body centered)
C     7.....R   
C     8.....H   
C     9.....TRIGONAL_NOT_RHOMOHEDRAL (this is basically hexagonal)

c     ------------------------------------------------------------------------
c     **** **** WIENXX definition::
c     ------------------------------------------------------------------------
c     H, R, P, C(monoclinic) --> correspond to PRIMITIVE reciprocal vectors
c              ^^^^^^^^^^^^  is this really true !!!
c     F, B, C(orthorhombic)  --> correspond to CONVENTIONAL reciprocal vectors

      if ( (igroup.eq.2 .or. igroup.eq.4 .or.
     $     igroup.eq.5 .or. igroup.eq.6) ) then
c
c     CASE: A, C, F, I(=body-centered)
c     
c     bug resolved by Florent Boucher
         if ( input_vec(1:4).eq.'prim' ) then
            IConvCoor=IGET_CONV
         else
            IConvCoor=IDO_NOTHING
         endif
         RETURN

      elseif ( igroup.eq.3 ) then
c
c     CASE: B, i.e, CXZ centered
c
c     check if it is orthorombic or monoclinic !!!
c     if gamma.eq.90.0d0 -> scalar product is zero !!!
         scalar = rcv(1,1)*rcv(2,1) + rcv(1,2)*rcv(2,2) +
     $        rcv(1,3)*rcv(2,3)
         if ( abs(scalar).lt.tolmin .and.
     $        input_vec(1:4).eq.'prim') then
c     it is orthorombic
            IConvCoor=IGET_CONV
         else if ( abs(scalar).ge.tolmin .and.
     $           input_vec(1:4).eq.'conv' ) then            
c     it is monoclinic
            IConvCoor=IGET_PRIM
         endif
         RETURN

      else
c
c     CASE: P, H, R
c
         if ( input_vec(1:4).eq.'conv' ) then
            IConvCoor=IGET_PRIM
            RETURN
         endif
      endif      

      IConvCoor=IDO_NOTHING
      RETURN
      END
      

      subroutine GetIntCoor(fr, icoo, mat, idv, kp)
      include 'param.inc'
      INTEGER icoo(3,NAC)
      REAL*8  fr(3,NAC), coo(3,NAC), mat(3,3)

      do i=1,kp
         coo(1,i)=fr(1,i)*mat(1,1) + fr(2,i)*mat(1,2) + fr(3,i)*mat(1,3)
         coo(2,i)=fr(1,i)*mat(2,1) + fr(2,i)*mat(2,2) + fr(3,i)*mat(2,3)
         coo(3,i)=fr(1,i)*mat(3,1) + fr(2,i)*mat(3,2) + fr(3,i)*mat(3,3)
c     print *, 'FR:: ',fr(1,i),fr(2,i),fr(3,i)
c     print *, 'COO::   ',coo(1,i),coo(2,i),coo(3,i)
      enddo

c     determine new IDV
      idv = IGetIDV(coo, kp)
c      print *,'IDV::',idv
c     get new icoo
      do i=1,kp
         do j=1,3
            icoo(j,i) = nint( dble(idv) * coo(j,i) )            
         enddo
c         print *, 'ICOO:: ',icoo(1,i),icoo(2,i),icoo(3,i)
      enddo
      END
      

      Integer Function IGetIDV(coor, kp)
      include 'param.inc'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 coor(3,NAC), min
      INTEGER imen(300)
      PARAMETER (TOLMIN = 1.0d-5)

      do i=1,kp
         n=(i-1)*3
         do 100 j=1,3
            imen(n+j) = 0
            imin      = 1
            min       = 1.0d0
c     if coor > 1.0 -> subtract integer part!!!
            p         = abs(coor(j,i)) - dint(coor(j,i)) 
c            print *, 'POINT(',i,',',j,'):: ',p            
            if ( p.lt.TOLMIN .or. abs(p-1.0d0).lt.TOLMIN ) then
               imen(n+j) = 1
               goto 100
            endif
            do ii=1,100    !100 should be more then enough for special points
               do jj=1,ii
                  f = abs( p - dble(jj)/dble(ii) )
c                  print *,'F::',f
                  if ( f .lt. min ) then
                     min  = f
                     imin = ii
                  endif
                  if ( f .lt. TOLMIN ) then
                     imen(n+j) = ii
                     goto 100
                  endif
               enddo
            enddo
            if ( imen(n+j) .eq. 0 ) imen(n+j)=imin
c            print *, 'IMEN:: ',imen(n+j),imin
c            print *,' '
 100     continue
      enddo

c     *** now find the minimum common multiplier
 110  continue
      iss = imen(1)
      idv = iss
      n   = 3*kp
      do i=2,n
         iss = idv
         idv = IGetMinMult( iss, imen(i) )
c         print *, 'IGETMULT:: ', idv
      enddo      
      IGetIDV = idv
      RETURN
      END

      Integer Function IGetMinMult( ia, ib)
      if ( ia.eq.0 .and. ib.eq.0 ) then
         IGetMinMult=1
         RETURN
      endif
c     suppose that: ia>ib --> check if it is
      if ( ib .gt. ia ) then
         ic = ia
         ia = ib
         ib = ic
      endif

      if ( ib .eq. 0 ) then
         IGetMinMult=ia
         RETURN
      endif

      ic = 1
c     now: ia is greater than ib */
      do i=1,ib
         if ( mod(ia*ib,ib) .eq. 0 ) then
            ic = i
            goto 120
         endif
      enddo      
 120  CONTINUE
      IGetMinMult=ia*ic
      RETURN
      END
