c     ------------------------------------------------------------------------
      program pwi2xsf_old
c     reads preprocesed (with pwi2xsf.sh) PW-input file (OLD versions ~ 1.0)
c     and converts to XSF format file
c     
c     Usage: pwi2xsf.sh < PW-preprocessed file
c     ------------------------------------------------------------------------

      implicit none

      integer maxtyp
      real*8
     $     bohr
      parameter (maxtyp = 100, bohr = 0.529177d0)

      integer
     $     ibrav,               ! label for Bravais lattice
     $     nat,                 ! number of atoms
     $     ntyp                 ! number of pseudopotentials

      real*8
     $     celldm(6),           ! cell parameters
     $     alat                 ! lattice parameter

      character
     $     dummy*80

      integer
     $     ityp,                ! type of PP
     $     atn(maxtyp),          ! nuclear charge
     $     ounit,               ! output unit
     $     i, j                 ! dummies

      real*8
     $     tau( 3 ),
     +     p( 3,3 ),            ! lattice vectors (PRIMITIVE)
     +     c( 3,3 )             ! lattice vectors (CONVETIONAL)

      logical
     $     ltaucry

      namelist/input/ ibrav, nat, celldm, ltaucry
      
      ltaucry = .false.

      ounit=6

c     set default values
      nat       = 0
      ibrav     = 0
      celldm(1) = 0.0d0
      do i=1,maxtyp
         atn(i) = 0
      enddo

      open(1, file='nuclei.charges', status='old')
      read (1,*) ntyp
      do i=1,ntyp
         read (1,*) j,atn(j)
      enddo
      close(1)

      read (5,input)

      if ( nat.eq.0 .or. celldm(1).eq.0.0d0 ) then
         print *,'ERROR reading INPUT.    STOPPING !!!'
         STOP
      endif

      if ( ibrav.eq.0 ) then
c     read custom lattice
         read (5,*) ((p(i,j),i=1,3),j=1,3)
         read (5,*) dummy      
         call dcopy(p,c,9)
      else
         call pwLatgen( ibrav, celldm,
     $        p(1,1), p(1,2), p(1,3), c(1,1), c(1,2), c(1,3) )
      endif

      alat = bohr*celldm(1)
      call write_XSF_header (alat, p, c, nat, ounit)

      do i=1,nat
         read(5,*) tau(1), tau(2), tau(3), ityp
         if (ltaucry) call pwCryst_to_cart(1, tau, p, 1)
         write(ounit,'(i3,2x,3f15.10)') atn(ityp),
     $        alat*tau(1), alat*tau(2), alat*tau(3)
      enddo
      END

c     ------------------------------------------------------------------------
      subroutine write_XSF_header (alat, p, c, nat, ounit)
c     writes the header for XSF structure file
c     ------------------------------------------------------------------------
      real*8
     $     alat,                ! lattice parameter
     $     p(3,3), c(3,3),      ! lattive vectors (PRIMITIVE & CONVETIONAL)
     $     p1(3,3), c1(3,3)     ! lattive vectors in ANGSTROMS unit
      integer
     $     nat,                 ! number of atoms
     $     ounit                ! output unit     
      integer
     $     i, j                 ! dummies

      do i=1,3
         do j=1,3
            p1(i,j) = alat*p(i,j)
            c1(i,j) = alat*c(i,j)
         enddo
      enddo

      write(ounit,'('' DIM-GROUP'')')
      write(ounit,*) 3, 1
      write(ounit,'(/,'' PRIMVEC'')')
      write(ounit,'(3(f15.10,2x,f15.10,2x,f15.10,/))')
     $     ((p1(i,j),i=1,3),j=1,3)
      write(ounit,'('' CONVVEC'')')
      write(ounit,'(3(f15.10,2x,f15.10,2x,f15.10,/))')
     $     ((c1(i,j),i=1,3),j=1,3)
      write(ounit,'('' PRIMCOORD'')')
      write(ounit,*) nat, 1
      return
      end

c     ------------------------------------------------------------------------
      subroutine dcopy(src,dst,n)
c     ------------------------------------------------------------------------
      real*8
     $     src(*), dst(*)
      integer
     $     n, i
      do i=1,n
         dst(i) = src(i)
      enddo
      return
      end

      
c This is basically a rewrite of Quantum ESPRESSO's latgen routine,
c rewritten to suite current needs.
c
c This file is distributed under the terms of the GNU General Public
c License. See the file `License' in the root directory of the present
c distribution, or http://www.gnu.org/copyleft/gpl.txt .
c
c
c---------------------------------------------------------------------
      subroutine pwLatgen( ibrav, celldm, p1, p2, p3, c1, c2, c3 )
c---------------------------------------------------------------------
c
c   This routine sets up the crystallographic vectors p1, p2, and p3.
c   The a's are expressed in units of celldm(1) ( a0 for
c   cubic lattices). This version contains all 14 Bravais lattices.
c
      implicit none
c
c     First the input variables
c
      real*8 
     +     celldm( 6 ),         ! input: the dimensions of the lattice
     +     p1( 3 ),             ! output: the first lattice vector (PRIMITIVE)
     +     p2( 3 ),             ! output: the second lattice vector
     +     p3( 3 ),             ! output: the third lattice vector
     +     c1( 3 ),             ! output: the first lattice vector(CONVETIONAL)
     +     c2( 3 ),             ! output: the second lattice vector
     +     c3( 3 )              ! output: the third lattice vector
      integer
     +       ibrav          ! input: the index of the Bravais lattice
c
c    Here the local variables required by the routine
c
      real*8
     +        sr2,          ! the square root of 2.0 
     +        sr3           ! the square root of 3.0
      parameter (sr2 = 1.4142 13562 373d0, sr3 = 1.7320 50807 569d0)
      real*8 
     +        term,         !\
     +        term1,        ! \ 
     +        term2,        !   Auxiliary variables
     +        cbya,         ! /
     +        singam,       !/ 
     +        sin,          ! the sine function
     +        sen
      integer 
     +        ipol          ! counter on the coordinates 

      if ( celldm( 1 ) .le. 0.d0 ) 
     +     call pwError( 'latgen', 'wrong celldm', max(1,ibrav) )

      if ( ibrav .eq. 0 ) then
c
c     user supplied lattice
c
         do  ipol = 1, 3
            p1( ipol ) = 0.d0
            p2( ipol ) = 0.d0
            p3( ipol ) = 0.d0
            
            c1( ipol ) = 0.d0
            c2( ipol ) = 0.d0
            c3( ipol ) = 0.d0
         end do
      endif

      if ( ibrav .eq. 1 ) then
c
c     simple cubic lattice
c

         p1( 1 ) = 1.0d0
         p2( 2 ) = 1.0d0
         p3( 3 ) = 1.0d0

      elseif ( ibrav .eq. 2 ) then
c
c     fcc lattice
c
         term = 1.0d0 / 2.d0

         p1( 1 ) = - term
         p1( 3 ) = term
         p2( 2 ) = term
         p2( 3 ) = term
         p3( 1 ) = - term
         p3( 2 ) = term

         c1( 1 ) = 1.0d0
         c2( 2 ) = 1.0d0
         c3( 3 ) = 1.0d0
c     **********
         return

      elseif ( ibrav .eq. 3 ) then
c
c     bcc lattice
c
         term = 1.0d0 / 2.d0

         do ipol = 1, 3
            p1( ipol ) = term
            p2( ipol ) = term
            p3( ipol ) = term
         enddo
         p2( 1 ) = - term
         p3( 1 ) = - term
         p3( 2 ) = - term

         c1( 1 ) = 1.0d0
         c2( 2 ) = 1.0d0
         c3( 3 ) = 1.0d0
c     **********
         return

      elseif ( ibrav .eq. 4 ) then
c
c     hexagonal lattice
c
         if ( celldm( 3 ) .le. 0.d0 ) 
     +        call pwError( 'latgen', 'wrong celldm', 4 )

         cbya = celldm( 3 )
         p1( 1 ) = 1.0d0
         p2( 1 ) = -1.0d0 / 2.d0
         p2( 2 ) = sr3 / 2.d0
         p3( 3 ) = cbya

      elseif ( abs(ibrav) .eq. 5 ) then
c
c     trigonal lattice
c
         if ( celldm( 4 ) .le. -0.5d0  .or.  celldm (4) .ge. 1.0d0 ) 
     +      call pwError( 'latgen', 'wrong celldm', 5 )

         term1=sqrt(1.d0 + 2.d0*celldm(4))
         term2=sqrt(1.d0 - celldm(4))

         if ( ibrav == 5 ) then
!     threefold axis along c (001)
            p2(2) = sr2*term2/sr3
            p2(3) = term1/sr3
            p1(1) = term2/sr2
            p1(2) =-p1(1)/sr3
            p1(3) = p2(3)
            p3(1) =-p1(1)
            p3(2) = p1(2)
            p3(3) = p2(3)
         elseif ( ibrav == -5 ) then
!     threefold axis along (111)
            p1(1) = (term1 - 2.0d0*term2) / 3.0d0
            p1(2) = (term1 + term2) / 3.0d0
            p1(3) = p1(2)
            p2(1) = p1(3)
            p2(2) = p1(1)
            p2(3) = p1(2)
            p3(1) = p1(2)
            p3(2) = p1(3)
            p3(3) = p1(1)
         endif

      elseif ( ibrav .eq. 6 ) then
c
c     tetragonal lattice
c
         if ( celldm( 3 ) .le. 0.d0 ) 
     +      call pwError( 'latgen', 'wrong celldm', 6 )

         cbya = celldm( 3 )

         p1( 1 ) = 1.0d0
         p2( 2 ) = 1.0d0
         p3( 3 ) = cbya

      elseif ( ibrav .eq. 7 ) then
c
c     body centered tetragonal lattice
c
         if ( celldm( 3 ) .le. 0.d0 ) 
     +      call pwError( 'latgen', 'wrong celldm', 7 )

         cbya=celldm(3)

         p2(1)=1.0d0/2.d0
         p2(2)=p2(1)
         p2(3)=cbya/2.d0
         p1(1)= p2(1)
         p1(2)=-p2(1)
         p1(3)= p2(3)
         p3(1)=-p2(1)
         p3(2)=-p2(1)
         p3(3)= p2(3)

         c1(1) = 1.0d0
         c2(2) = 1.0d0
         c3(3) = cbya
c     **********
         return

      elseif ( ibrav .eq. 8 ) then
c
c     Simple orthorombic lattice
c
         if ( celldm( 2 ) .le. 0.d0 .or. celldm( 3 ) .le. 0.d0 ) 
     $        call pwError( 'latgen', 'wrong celldm', 8 )
         p1( 1 ) = 1.0d0
         p2( 2 ) = celldm( 2 )
         p3( 3 ) = celldm( 3 )

      elseif ( abs(ibrav) .eq. 9 ) then
c
c     One face centered orthorombic lattice
c
         if ( celldm( 1 ) .le. 0.d0 .or. celldm( 2 ) .le. 0.d0
     +        .or. celldm( 3 ) .le. 0.d0 ) 
     +        call pwError( 'latgen', 'wrong celldm', 9 )
         
         if ( ibrav .eq. 9 ) then
!     old PWscf description
            
            p1( 1 ) = 0.5d0
            p1( 2 ) = p1(1) * celldm( 2 )
            p2( 1 ) = - p1( 1 )
            p2( 2 ) = p1( 2 )
            
         else
!     alternate description
            p1(1) = 0.5d0
            p1(2) =-p1(1) * celldm(2)
            p2(1) = p1(1)
            p2(2) =-p1(2)            
         endif
         p3( 3 ) = celldm( 3 )  

         c1(1) = 1.0d0
         c2(2) = celldm(2)
         c3(3) = celldm(3)
c     ********** 
         return

      elseif ( ibrav .eq. 10 ) then
c
c     All face centered orthorombic lattice
c
         if ( celldm( 2 ) .le. 0.d0 .or. celldm( 3 ) .le. 0.d0 ) 
     +        call pwError( 'latgen', 'wrong celldm', 10 )

         p2( 1 ) = 0.5d0
         p2( 2 ) = 0.5d0 * celldm( 2 )
         p1( 1 ) = p2( 1 )
         p1( 3 ) = 0.5d0 * celldm( 3 )
         p3( 2 ) = 0.5d0 * celldm( 2 )
         p3( 3 ) = p1( 3 )

         c1(1) = 1.0d0
         c2(2) = celldm(2)
         c3(3) = celldm(3)
c     **********
         return

      elseif ( ibrav .eq. 11 ) then
c
c     Body centered orthorombic lattice *** HERE
c
         if ( celldm( 2 ) .le. 0.d0 .or. celldm( 3 ) .le. 0.d0 ) 
     +        call pwError( 'latgen', 'wrong celldm', 11 )

         p1( 1 ) = 0.5d0
         p1( 2 ) = 0.5d0 * celldm( 2 )
         p1( 3 ) = 0.5d0 * celldm( 3 )
         p2( 1 ) = - p1( 1 )
         p2( 2 ) = p1( 2 )
         p2( 3 ) = p1( 3 )
         p3( 1 ) = - p1( 1 )
         p3( 2 ) = - p1( 2 )
         p3( 3 ) = p1( 3 )

         c1(1) = 1.0d0
         c2(2) = celldm(2)
         c3(3) = celldm(3)
c     **********
         return

      elseif ( ibrav .eq. 12 ) then
c
c     Simple monoclinic lattice, unique (i.e. orthogonal to a) axis: c
c
         if ( celldm( 2 ) .le. 0.d0 .or. celldm( 3 ) .le. 0.d0
     $        .or. abs(celldm(4)) .ge. 1.d0 )
     $        call pwError( 'latgen', 'wrong celldm', 12 )

         sen = sqrt( 1.d0 - celldm( 4 ) ** 2 )

         p1( 1 ) = 1.0d0
         p2( 1 ) = celldm( 2 ) * celldm( 4 )
         p2( 2 ) = celldm( 2 ) * sen
         p3( 3 ) = celldm( 3 )

      else if (ibrav ==-12) then
!     
!        Simple monoclinic lattice, unique axis: b (more common)
!        
         if (celldm(2) .le. 0.d0)
     $        call pwError ('latgen', 'wrong celldm(2)', ibrav)
         if (celldm(3) .le. 0.d0)
     $        call pwError ('latgen', 'wrong celldm(3)', ibrav)
         if (abs(celldm(5)).ge.1.d0)
     $        call pwError ('latgen', 'wrong celldm(5)', ibrav)      
         
         sen=sqrt(1.d0-celldm(5)**2)
         
         p1(1) = 1.0d0
         p2(2) = celldm(2)
         p3(1) = celldm(3)*celldm(5)
         p3(3) = celldm(3)*sen     

      elseif ( ibrav .eq. 13 ) then
c
c     One face centered monoclinic lattice
c
         if ( celldm( 2 ) .le. 0.d0
     +        .or. celldm( 3 ) .le. 0.d0 
     +        .or. abs ( celldm( 4 ) ) .ge. 1.d0 ) 
     +        call pwError( 'latgen', 'wrong celldm', 13 )

         sen = sqrt( 1.d0 - celldm(4) ** 2 )

         p1(1) = 0.5 
         p1(3) =-p1(1) * celldm(3)
         p2(1) = celldm(2) * celldm(4)
         p2(2) = celldm(2) * sen
         p3(1) = p1(1)
         p3(3) =-p1(3)
         
         c1(1) = 1.0d0
         c2(1) = p2(1)
         c2(2) = p2(2)
         c3(3) = celldm(3)

         return

      elseif ( ibrav .eq. 14 ) then
c
c     Triclinic lattice
c
         if (      celldm( 2 ) .le. 0.d0
     +        .or. celldm( 3 ) .le. 0.d0
     +        .or. abs ( celldm( 4 ) ) .ge. 1.d0
     +        .or. abs ( celldm( 5 ) ) .ge. 1.d0
     +        .or. abs ( celldm( 6 ) ) .ge. 1.d0 )
     +      call pwError( 'latgen', 'wrong celldm', 14 )

         singam = sqrt( 1.d0 - celldm( 6 ) ** 2 )

         term = 1.d0
     $        +  2.d0 * celldm( 4 ) * celldm( 5 ) * celldm( 6 )
     $        - celldm( 4 )**2 
     +        - celldm( 5 )**2
     $        - celldm( 6 )**2 

         if (term < 0.d0) call pwError('latgen',
     $        'celldm do not make sense, check your data', ibrav)

         term = sqrt(term / ( 1.d0 - celldm( 6 ) ** 2 ))
         
         p1( 1 ) = 1.0d0
         p2( 1 ) = celldm( 2 ) * celldm( 6 )
         p2( 2 ) = celldm( 2 ) * singam
         p3( 1 ) = celldm( 3 ) * celldm( 5 )
         p3( 2 ) = celldm( 3 ) *
     $        (celldm(4) - celldm(5)*celldm(6)) / singam
         p3( 3 ) = celldm( 3 ) * term

      else

         call pwError('latgen', ' wrong ibrav ', ibrav )

      endif

c     **********
c     if we came here then just copy p vectors to c vectors !!!
c     **********

      do  ipol = 1, 3
         c1( ipol ) = p1( ipol )
         c2( ipol ) = p2( ipol )
         c3( ipol ) = p3( ipol )
      enddo

      return
      end



c
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


!
! This is taken from an old version of Quantum ESPRESSO's errore
! routine.
!
! This file is distributed under the terms of the GNU General Public
! License. See the file `License' in the root directory of the present
! distribution, or http://www.gnu.org/copyleft/gpl.txt .
!
c----------------------------------------------------------------------
      subroutine pwError( routin, messag, ierr )
c----------------------------------------------------------------------
c
c    This is a simple routine which writes an error message to
c    output. If ierr = 0 it does nothing, if ierr < 0 it
c    writes the message but does not stop, if ierr > 0 stops.
c
      implicit none


      character*(*) 
     +          routin,       ! the name of the calling routine
     +          messag        ! the output message
      integer 
     +          ierr          ! the error flag

      if(ierr.eq.0) return
      write(6,*) ' '
      write(6,'(1x,78(''%''))')
      write(*,'(5x,''from '',a,'' : error #'',i10)') routin, ierr
      write(*,'(5x,a)') messag
      write(6,'(1x,78(''%''))')

      if(ierr.gt.0) then
         write(*,'(''     stopping ...'')')
         stop
      else
         write(6,*) ' '
         return
      end if
      end
      
