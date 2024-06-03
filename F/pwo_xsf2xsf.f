!-----------------------------------------------------------------------
      PROGRAM pwo_xsf2xsf
!-----------------------------------------------------------------------
!
!     Usage: pwo_xsf2xsf ibrav 
!     
!     This program reads the XSF file produced by pwo2xsf.sh filter and
!     tries to assign CONVVEC on the basis of provided "ibrav" index,
!     which is specified as command-line argument.
!
!     The resulting XSF file is written to standard output.
!-----------------------------------------------------------------------
      IMPLICIT NONE
      integer narg, ibrav, natoms, nsteps, i, j, ia
      integer igroup, ibrav2igroup
      
      character*256 line, abrav
      character atm*4, keyword*32
      real*8 pv(3,3), cv(3,3), celldm(6), tau(3), VecSize
      logical zero_latvec, matches
      
      call ZeroMat(cv,3,3)
      
      narg = iargc()
      if (narg .ne. 1)
     $     stop 'Usage: pwo_xsf2xsf ibrav'

      call getarg(1,abrav)
      read(abrav,*) ibrav
      
      igroup = ibrav2igroup(ibrav)
      
      nsteps = 1
      
!     ANIMSTEPS or CRYSTAL|SLAB|POLYMER|MOLECULE
      call read_line(line)
         
      IF ( matches('ANIMSTEPS', line) ) THEN
         write (*,*)  trim(line)
         read(line,*) keyword, nsteps
         
!     CRYSTAL|SLAB|POLYMER|MOLECULE
         call read_line(line)
      ENDIF

      IF ( matches('CRYSTAL', line) ) THEN
         write(*,'('' DIM-GROUP'')')
         write(*,'('' 3 '',i3)') igroup
      ELSE
         write (*,*)  trim(line)
      ENDIF
      
      DO ia=1,nsteps
         
!     PRIMVEC or PRIMCOOR
         call readwrite_line(line)
!
         IF ( matches('PRIMVEC', line) ) THEN
            read (*,*) ((pv(i,j),i=1,3),j=1,3)
            write(*,'(3(f15.10,2x,f15.10,2x,f15.10,/))')
     $           ((pv(i,j),i=1,3),j=1,3)
            
!     try to assign CONVVEC
            call at2celldm
     $           (ibrav, pv(1,1), pv(1,2), pv(1,3), celldm)
            call set_convvec(ibrav, celldm, pv, cv)
            
!     write CONVENTIONAL VECTORS is they are not zero
            IF ( .not. zero_latvec(cv) ) THEN
!     
               write(*,'('' CONVVEC'')')
               write(*,'(3(f15.10,2x,f15.10,2x,f15.10,/))')
     $              ((cv(i,j),i=1,3),j=1,3)
            ENDIF
      
!     read PRIMCOOR line
            call readwrite_line(line)
         ENDIF
         
         IF ( matches('PRIMCOOR', line) ) THEN
            call readwrite_line(line) ! natoms, 1
            read (line,*) natoms

            do i=1,natoms
               call readwrite_line(line)
            enddo
            write(*,*) ''
         ENDIF
      ENDDO
      
      END
!-----------------------------------------------------------------------

      
      SUBROUTINE readwrite_line(line)
      character*256 line
      
 10   continue
      read (*,'(a256)') line
      line = adjustl(line)

!     check for empty line
      if ( len_trim(line) .eq. 0) goto 10

      write (*,*) trim(line)

      return
      END
     
      SUBROUTINE read_line(line)
      character*256 line
      
 10   continue
      read (*,'(a256)') line
      line = adjustl(line)

!     check for empty line
      if ( len_trim(line) .eq. 0) goto 10

      return
      END
 
