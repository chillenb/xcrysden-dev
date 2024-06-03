c     ------------------------------------------------------------------------
      program pwi2xsf
c     Reads preprocesed (with pwi2xsf.sh) PW-input file
c     and converts to XSF format file
c
c     This program reads the NEWLY formated preprocesed-PW.X input
c     
c     Usage: pwi2xsf.sh < PW-preprocessed file
c     ------------------------------------------------------------------------

      implicit none

      integer
     $     maxtyp,              ! maximum number of different types of atoms
     $     maxatom,             ! maximum number fo atoms
     $     maximage,            ! maximum number of images
     $     NOT_SPECIFIED_UNIT,
     $     ALAT_UNIT,
     $     BOHR_UNIT,
     $     ANGSTROM_UNIT,
     $     CRYSTAL_UNIT

      real*8
     $     bohr

      parameter (
     $     maxtyp  = 100,
     $     maxatom = 10000,
     $     maximage= 50,
     $     bohr    = 0.52917720859d0,
     $
     $     NOT_SPECIFIED_UNIT = 0,
     $     ALAT_UNIT     = 1,
     $     BOHR_UNIT     = 2,
     $     ANGSTROM_UNIT = 3,
     $     CRYSTAL_UNIT  = 4 )          

      integer
     $     ibrav,               ! label for Bravais lattice
     $     igroup,              ! the XSF's group label of Bravais lattice
     $     nat,                 ! number of atoms
     $     ntyp,                ! number of pseudopotentials
     $     num_of_images,       ! number of NEB images
     $     inp_num_of_images,   ! number of NEB images in the input
     $     cell_unit,           ! length-unit of cell-parameters
     $     atomic_posunit(maximage) ! length-unit of atomic positions

      real*8
     $     celldm(6),           ! cell parameters
     $     alat,                ! lattice parameter
     $     a, b, c, cosab, cosac, cosbc, ! lattice parameters
     $     omega ! unit cell volume (not used)

      character
     $     calculation*80,      ! type of calculation
     $     line*120             ! line of input
      character*3
     $     atm(maxatom,maximage) ! atomic symbols

      integer
     $     ityp,                ! type of PP
     $     ounit,               ! output unit     
     $     i, j, ipol,          ! dummies
     $     inat, iim, iim_old,  ! counters
     $     length, string_length ! length of string

      real*8
     $     x,y,z,               ! Cartesian coordinates & weights
     $     w1,w2,               ! linear interpolation weights
     $     dx, dy, dz,          ! auxiliary
     $     tau(3,maxatom,maximage), ! atomic coordinates
     +     pv( 3,3 ),           ! lattice vectors (PRIMITIVE)
     +     cv( 3,3 ),           ! lattice vectors (CONVENTIONAL)
     $     old_total_dist, old_dist(maximage), ! old(=input) inter-image distances
     $     new_total_dist, new_dist, ! new(=output) inter-image distances
     $     scale                ! auxiliary

      logical
     $     ltaucry, matches, last_image

      data cv/
     $     0.0d0, 0.0d0, 0.0d0,
     $     0.0d0, 0.0d0, 0.0d0,
     $     0.0d0, 0.0d0, 0.0d0/

      namelist/system/
     $     ibrav, nat, celldm, a, b, c, cosab, cosac, cosbc,
     $     calculation, num_of_images
      
      ounit = 6
      iim = 0
      last_image = .false.

c     set default values
      calculation   = 'scf'
      num_of_images = 1
      nat       = 0
      ibrav     = 0
      celldm(1) = 0.0d0
      a         = 0.0D0 
      b         = 0.0D0
      c         = 0.0D0
      cosab     = 0.0D0
      cosac     = 0.0D0
      cosbc     = 0.0D0

      
c-----------------------------------------------------------------------     
c     read namelist system
c-----------------------------------------------------------------------
      read (5,system)
c-----------------------------------------------------------------------
      
      if ( nat.eq.0 ) then
         call errore('pwi2xsf', 'error while reading pw.x input',1)
      endif

c     was lattice specified in terms of A, B, C,...
      if ( celldm(1) .eq. 0.0D0 .AND. a .ne. 0.0D0 ) THEN
         CALL abc2celldm ( ibrav, a,b,c,cosab,cosac,cosbc, celldm )
      else if ( celldm(1) .ne. 0.0D0 .AND. a .ne. 0.0D0 ) THEN
         !
         call errore('pwi2xsf', 'both celldm() and a,b,c specified', 2)
         !
      endif
      
c
c     read the rest of the input
c
 990  continue
      read(5,'(a120)',end=999) line
      line = adjustl(line)
      length = string_length( line )

c     
c     CELL_PARAMETERS
c     
      if ( line(1:15) .eq. 'CELL_PARAMETERS' ) then

c     find out the length-unit ...
         line(1:length-15) = adjustl(line(16:length))
         length = string_length (line)

         cell_unit = NOT_SPECIFIED_UNIT

         if ( celldm(1).eq.0.0d0 ) then
            celldm(1) = 1.0d0
         endif      

         if (length.gt.0 ) then
            if ( matches('ALAT',line) ) then
               cell_unit = ALAT_UNIT
            elseif ( matches('BOHR',line) ) then
               cell_unit = BOHR_UNIT
            elseif ( matches('ANGSTROM',line) ) then
               cell_unit = ANGSTROM_UNIT
            elseif ( matches('CUBIC',line) .or.
     $               matches('HEXA', line) ) then 
               ! old style syntax
               cell_unit = ALAT_UNIT
!            else
!               print *, 'ERROR: wrong unit in CELL_PARAMETERS'
!               print *, '       must be alat, bohr, or angstrom'
!               STOP
            endif
         endif
         
         read (5,*) ((pv(i,j),i=1,3),j=1,3)

         scale = 1.0d0

         if ( (cell_unit .eq. ALAT_UNIT) .or.
     $        (cell_unit .eq. NOT_SPECIFIED_UNIT) ) then
            
            scale = celldm(1)   ! lattice is specified in alat units
            
         elseif ( cell_unit .eq. ANGSTROM_UNIT ) then
            
            scale = 1.0d0/bohr
            
         elseif ( cell_unit .eq. BOHR_UNIT ) then
            
            scale = 1.0d0
            
         endif

c     convert pv() lattice vectors to atomic units
         
         do i=1,3
            pv(1,i) = pv(1,i)*scale
            pv(2,i) = pv(2,i)*scale
            pv(3,i) = pv(3,i)*scale
         enddo
         !
      elseif ( line(1:16) .eq. 'ATOMIC_POSITIONS' ) then
c     
c     ATOMIC_POSITIONS
c
c     find out the length-unit
         line(1:length-16) = adjustl(line(17:length))
         length = string_length (line)
         !
         iim = iim + 1
         !
         atomic_posunit(iim) = ALAT_UNIT
         !
         if (length.gt.0 ) then
            !
            if ( matches('CRYSTAL_SG',line) ) then
               call errore('pwi2xsf',
     &              'crystal_sg currently not supported', 3)
               !
            elseif ( matches('ALAT',line) ) then
               atomic_posunit(iim) = ALAT_UNIT
            elseif ( matches('BOHR',line) ) then
               atomic_posunit(iim) = BOHR_UNIT            
            elseif ( matches('CRYSTAL',line) ) then
               atomic_posunit(iim) = CRYSTAL_UNIT
            elseif ( matches('ANGSTROM',line) ) then
               atomic_posunit(iim) = ANGSTROM_UNIT
            endif
         endif

c     
c     READ ATOMS
c     
         IF ( (calculation(1:3) .ne. 'NEB') .and.
     $        (calculation(1:3) .ne. 'SMD') ) THEN

            if ( (iim .gt. 1) .and.
     $           (calculation(1:4) .ne. 'PATH') ) then
               call errore('pwi2xsf',
     &              'several ATOMIC_POSITIONS cards 
     &for non PATH calculation', 4)
            endif
            !
            call read_atoms(nat,atm(1,iim),tau(1,1,iim))
            !
         ELSE
c     
c     old format path calculation (NEB or SMD)
c     
            last_image = .false.
            iim = iim - 1
            do while(.not.last_image)
               !
               iim = iim + 1
               !
               read (5,'(a120)') line ! line: first_image
               if ( matches('LAST_IMAGE',line) ) last_image = .true.            
               !
               if (iim .gt. 1)
     $              atomic_posunit(iim) = atomic_posunit(iim-1)
               !
               call read_atoms(nat,atm(1,iim),tau(1,1,iim))
            enddo
         ENDIF
      endif
      !
      inp_num_of_images = iim
      goto 990
      !
 999  continue
      
      if ( ibrav.ne.0 ) then
         
!     get primitive-lattice vectors in atomic units         
         call latgen(ibrav,celldm, pv(1,1), pv(1,2), pv(1,3), omega)
         
!     get conventional-lattice vectors in atomic units         
         call set_convvec(ibrav, celldm, pv, cv)

      endif

      call write_XSF_header
     &     (num_of_images, ibrav, bohr, pv, cv, nat, ounit)
      
c     coordinates to ANGSTROMs

      alat = bohr*celldm(1)     ! alat in angstrom units

      do iim=1,inp_num_of_images
         do inat=1,nat
            if     ( atomic_posunit(iim) .eq. BOHR_UNIT ) then            
               tau(1,inat,iim) = bohr * tau(1,inat,iim)
               tau(2,inat,iim) = bohr * tau(2,inat,iim)
               tau(3,inat,iim) = bohr * tau(3,inat,iim)
               
            elseif ( atomic_posunit(iim) .eq. ALAT_UNIT ) then
               tau(1,inat,iim) = alat * tau(1,inat,iim)
               tau(2,inat,iim) = alat * tau(2,inat,iim)
               tau(3,inat,iim) = alat * tau(3,inat,iim)
               
            elseif ( atomic_posunit(iim) .eq. CRYSTAL_UNIT ) then
               call pwCryst_to_cart(1, tau(1,inat,iim), pv, 1)
!     pwCryst_to_cart returns coordinates in atomic units
               tau(1,inat,iim) = bohr * tau(1,inat,iim)
               tau(2,inat,iim) = bohr * tau(2,inat,iim)
               tau(3,inat,iim) = bohr * tau(3,inat,iim)
            endif
         enddo
      enddo

      if ( num_of_images .lt. 2 ) then
         
c     write atoms for non-PATH calculation
         do inat=1,nat
            write(ounit,'(a3,2x,3f15.10)') atm(inat,1),
     $           tau(1,inat,1), tau(2,inat,1), tau(3,inat,1)
         enddo
         
      else

c     calculate intermediate images for PATH calculation

         old_total_dist = 0.0d0
         old_dist(1)    = 0.0d0
         
         do iim = 2, inp_num_of_images

            old_dist(iim) = 0.0
            
            do inat=1,nat
               dx = tau(1,inat,iim) - tau(1,inat,iim-1)
               dy = tau(2,inat,iim) - tau(2,inat,iim-1)
               dz = tau(3,inat,iim) - tau(3,inat,iim-1)
               old_dist(iim) = old_dist(iim) + dx*dx + dy*dy + dz*dz
            enddo

            old_dist(iim) = sqrt( old_dist(iim) )
            
            old_total_dist = old_total_dist + old_dist(iim)
            old_dist(iim)  = old_total_dist
         enddo
         
         new_dist = old_total_dist / dble(num_of_images-1)
     
c     --------------------------------------------------
c     perform INTERPOLATION
c     --------------------------------------------------
         
         new_total_dist = 0.0
         do iim=1,num_of_images-1
            do iim_old=1,inp_num_of_images-1            
               if ( new_total_dist .ge. old_dist(iim_old)
     $              .and.
     $              new_total_dist .lt. old_dist(iim_old+1) + 1d-10 )
     $              then
               
                  w1 = ( old_dist(iim_old+1) - new_total_dist )
     $                 /
     $                 ( old_dist(iim_old+1) - old_dist(iim_old) )
                  w2 = 1.0d0 - w1
                  
                  write(ounit,'('' PRIMCOORD '',i5)') iim
                  write(ounit,*) nat, 1

                  do inat=1,nat
                     x = w1*tau(1,inat,iim_old)+w2*tau(1,inat,iim_old+1)
                     y = w1*tau(2,inat,iim_old)+w2*tau(2,inat,iim_old+1)
                     z = w1*tau(3,inat,iim_old)+w2*tau(3,inat,iim_old+1)
                     write(ounit,'(a3,2x,3f15.10)')
     $                    atm(inat,iim_old), x, y, z
                  enddo
                  goto 11
               endif
            enddo
 11         continue
            new_total_dist = new_total_dist + new_dist
         enddo
      
c     print last image
         write(ounit,'('' PRIMCOORD '',i5)') iim
         write(ounit,*) nat, 1
         
         do inat=1,nat
            x = tau(1,inat,inp_num_of_images)
            y = tau(2,inat,inp_num_of_images)
            z = tau(3,inat,inp_num_of_images)
            write(ounit,'(a3,2x,3f15.10)')
     $           atm(inat,inp_num_of_images), x, y, z
         enddo
      endif
      END


c     ------------------------------------------------------------------------
      subroutine read_atoms(nat,atm,coor)
c     read atomic coordinates
c     ------------------------------------------------------------------------
      implicit none
      integer
     $     nat,                 ! number of atoms
     $     ipol,inat,           ! counters
     $     length, string_length ! length of a string
      character
     $     line*120             ! line of input
      character*3
     $     atm(*)               ! atomic symbols
      real*8
     $     coor(3,*)
      
      do inat=1,nat
 10      continue
         read (5,'(a120)') line
         length = string_length( line )
         if (length.eq.0) then
c     an empty line, read again
            goto 10
         endif
         read (line,*) atm(inat),(coor(ipol,inat),ipol=1,3)
      enddo
      return
      end


c     ------------------------------------------------------------------------
      SUBROUTINE write_XSF_header
     &     (num_of_images, ibrav, bohr, p, c, nat, ounit)
c     writes the header for XSF structure file
c     ------------------------------------------------------------------------
      implicit none
      
      real*8
     $     bohr,                ! Bohr radius in ANGSTROMS unit
     $     p(3,3), c(3,3),      ! lattive vectors in Bohr units (PRIMITIVE & CONVENTIONAL)
     $     p1(3,3), c1(3,3)     ! lattive vectors in ANGSTROMS unit
      integer
     $     ibrav,               ! label for Bravais lattice
     $     igroup,              ! group of DIM-GROUP
     $     ibrav2igroup,        ! function for transforming ibrav to igroup
     $     nat,                 ! number of atoms
     $     num_of_images,       ! number of NEB images
     $     ounit                ! output unit     
      integer
     $     i, j                 ! dummies
      logical zero_latvec

      do i=1,3
         do j=1,3
            p1(i,j) = bohr*p(i,j)
            c1(i,j) = bohr*c(i,j)
         enddo
      enddo

      igroup = ibrav2igroup(ibrav)
      
      if (num_of_images .gt. 1)
     $     write(ounit,'('' ANIMSTEPS '',i5)') num_of_images

      write(ounit,'('' DIM-GROUP'')')
      write(ounit,'('' 3 '',i3)') igroup
      write(ounit,'(/,'' PRIMVEC'')')
      write(ounit,'(3(f15.10,2x,f15.10,2x,f15.10,/))')
     $     ((p1(i,j),i=1,3),j=1,3)

!     write CONVENTIONAL VECTORS is they are not zero
      IF ( .not. zero_latvec(c1) ) THEN
!     
         write(ounit,'('' CONVVEC'')')
         write(ounit,'(3(f15.10,2x,f15.10,2x,f15.10,/))')
     $        ((c1(i,j),i=1,3),j=1,3)
      ENDIF
!
      if (num_of_images .eq. 1) then
         write(ounit,'('' PRIMCOORD'')')
         write(ounit,*) nat, 1
      endif
      return
      end
      
      
c-----------------------------------------------------------------------
      subroutine pwCryst_to_cart ( nvec, vec, trmat, iflag )
c-----------------------------------------------------------------------
c
c     This routine transforms the atomic positions or the k-point 
c     components from crystallographic to carthesian coordinates ( iflag=1)
c     and viceversa ( iflag=-1 ).
c     Output carth. coordinates are stored in the input ('vec') array. 
c
c
      implicit none
c
c     first the dummy variables
c
      integer 
     +        nvec,        ! input: number of vectors (atom. pos. or k-points)
     +                     !        to be transf. from cryst. to carth. axes 
     +        iflag        ! input: gives the sense of the transformation
      real*8 
     +        vec(3,nvec), ! input/output: cryst./carth. coord. of the vectors 
     +                     !               (atom. pos. or k-points)
     +        trmat(3,3)   ! input: transformation matrix
                           ! if iflag=1:
                           !    trmat = at ,  basis of the real-space latt. 
                           !                  for atoms   or 
                           !          = bg ,  basis of the rec.-space latt. 
                           !                  for k-points 
                           ! if iflag=-1: the opposite
c
c    here the local variables
c
      integer 
     +        nv,          ! counter on vectors 
     +        kpol         ! counter on polarizations

      real*8 
     +       vau(3)        ! auxil. vector (containing the temp. transf. coord.)
c
c     Compute the carth. coordinates of each vectors 
c     (atomic positions or k-points components)
c
      do nv = 1, nvec
         if ( iflag.eq.1 ) then
            do kpol = 1,3
               vau(kpol) = trmat(kpol,1)*vec(1,nv) + 
     +                     trmat(kpol,2)*vec(2,nv) +
     +                     trmat(kpol,3)*vec(3,nv)
            enddo
         else
            do kpol = 1,3
               vau(kpol) = trmat(1,kpol)*vec(1,nv) + 
     +                     trmat(2,kpol)*vec(2,nv) +
     +                     trmat(3,kpol)*vec(3,nv)
            enddo
         endif
         do kpol = 1,3
            vec(kpol,nv) = vau(kpol) 
         enddo
      enddo
c
      return
      end


