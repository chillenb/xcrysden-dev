!------------------------------------------------------------------------
      INTEGER FUNCTION ibrav2igroup(ibrav)
!-----------------------------------------------------------------------
      implicit NONE
      integer ibrav

      ibrav2igroup = 1

      if (abs(ibrav) .eq. 1)  ibrav2igroup = 1
      if (abs(ibrav) .eq. 2)  ibrav2igroup = 5
      if (abs(ibrav) .eq. 3)  ibrav2igroup = 6
      if (abs(ibrav) .eq. 4)  ibrav2igroup = 9
      if (abs(ibrav) .eq. 5)  ibrav2igroup = 7
      if (abs(ibrav) .eq. 6)  ibrav2igroup = 1
      if (abs(ibrav) .eq. 7)  ibrav2igroup = 6
      if (abs(ibrav) .eq. 8)  ibrav2igroup = 1
      if (abs(ibrav) .eq. 9)  ibrav2igroup = 4
      if (abs(ibrav) .eq. 91) ibrav2igroup = 4
      if (abs(ibrav) .eq. 10) ibrav2igroup = 5
      if (abs(ibrav) .eq. 11) ibrav2igroup = 6
      if (abs(ibrav) .eq. 12) ibrav2igroup = 1
      if (abs(ibrav) .eq. 13) ibrav2igroup = 4
      if (abs(ibrav) .eq. 14) ibrav2igroup = 1
      
      RETURN
      END      
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      LOGICAL FUNCTION zero_latvec(v)
!-----------------------------------------------------------------------
      implicit NONE
      real*8 v(3,3), VecSize
      
      zero_latvec = (
     &     (VecSize(v(1,1)) .lt. 1d-6) .and.
     &     (VecSize(v(1,2)) .lt. 1d-6) .and.
     &     (VecSize(v(1,3)) .lt. 1d-6) )
      return
      END
!-----------------------------------------------------------------------

      
!-----------------------------------------------------------------------
      SUBROUTINE set_convvec(ibrav, celldm, pv, cv)
!-----------------------------------------------------------------------
!     This routine is used for pw.x I/O files and generates the
!     conventional lattice vectors, cv(3,3), based on ibrav, celldm(),
!     and pv() data, where pv() are primitive lattice vectors
!
!     N.B:
!     * pv() must be in atomic Bohr units
!     * in output the cv() are also in atomic Bohr units
!------------------------------------------------------------------------      
      
      IMPLICIT NONE
      integer ibrav, i, j
      real*8 celldm(6), pv(3,3), cv(3,3)
      
!      CALL ZeroMat(cv, 3, 3)

!
!     only centered lattices have cv() different from pv()
!
      IF ( (ibrav .eq. 2) .or. (abs(ibrav) .eq. 3) ) THEN
!     
!     fcc (2) and bcc (3 or -3)
!     
         cv(1,1) = 1.0d0
         cv(2,2) = 1.0d0
         cv(3,3) = 1.0d0

      ELSEIF ( ibrav .eq. 7 ) THEN
!
!     tetragonal I (bct)
!
         cv(1,1) = 1.0d0
         cv(2,2) = 1.0d0
         cv(3,3) = celldm(3) ! C/A

      ELSEIF ( (abs(ibrav) .eq. 9)
     &        .or. (ibrav .eq. 91)
     &        .or. (ibrav .eq. 10)
     &        .or. (ibrav .eq. 11) ) THEN
!
!     Orthorhombic lattices:
!     * ibrav == 9 or -9:   1-face (C) centered orthorhombic
!     * ibrav == 91:        1-face (A) centered orthorhombic
!     * ibrav == 10:        face-centered 
!     * ibrav == 11:        body-centered
         cv(1,1) = 1.0d0
         cv(2,2) = celldm(2)     ! B/A
         cv(3,3) = celldm(3)     ! C/A

      ELSEIF ( ibrav .eq. 13 ) THEN
!
!     Monoclinic base-centered (unique axis c) 
!
!     cv(1,1) -- 1st vector
!     cv(1,2) -- 2nd vector
!     cv(1,3) -- 3rd vector
!     
         cv(1,1) = 1.0d0
         cv(1,2) = pv(1,2)/celldm(1)
         cv(2,2) = pv(2,2)/celldm(1)
         cv(3,3) = celldm(3)    ! C/A
         
      ELSEIF ( ibrav .eq. -13 ) THEN
!
!     Monoclinic base-centered (unique axis b)
!
!     cv(1,1) -- 1st vector
!     cv(1,2) -- 2nd vector
!     cv(1,3) -- 3rd vector
!     
         cv(1,1) = 1.0d0
         cv(2,2) = celldm(2)    ! B/A
         cv(1,3) = pv(1,3)/celldm(1)
         cv(3,3) = pv(3,3)/celldm(1)
      ENDIF

!
!     scale cv() with celldm(1) as to convert them to Bohr atomic units
!
      DO i=1,3
         DO j=1,3
            cv(j,i) = celldm(1) * cv(j,i)
         ENDDO
      ENDDO
      
      RETURN
      END
!-----------------------------------------------------------------


c     -----------------------------------------------------------------------
      LOGICAL FUNCTION matches (str1, str2)  
c     .true. if str1 is contained in str2, .false. otherwise
c     This function is taken from PWscf package (www.pwscf.org).
c     -----------------------------------------------------------------------
      implicit none  
      character str1*(*), str2*(*)  
      integer len1, len2, l  
      
      len1 = len(str1)  
      len2 = len(str2)  
      do l = 1, len2 - len1 + 1  
         if ( str1(1:len1) .eq. str2(l:l + len1 - 1) ) then  
            matches = .true.
            return
         endif
      enddo
      matches = .false.  
      return  
      end
      
