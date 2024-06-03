      subroutine WriteFTN34(
     $     title,        ! tittle (arbitrary)
     $     idim,         ! dimensionality of the system
     $     igroup,       ! group of the system
     $     vec,          ! primitive vectors
     $     nsymm,        ! n of symmetry operations
     $     symmop,       ! symmetry operators matrices
     $     symmtr,       ! translation of symmop
     $     natr,         ! number of atoms in primitive unit cell
     $     nat, x, y, z) ! atomic number & cartesian coordinates
      
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 vec(3,3), symmop(nsymm,3,3), symmtr(nsymm,3),
     $     x(natr), y(natr), z(natr)
      INTEGER nat(natr)
      CHARACTER*80 title

      Write(34, '(A)') title
      Write(34, *) idim, igroup
      do i=1,3
         Write(34,*) (vec(i,j), j=1,3)
      enddo
      Write(34,*) nsymm
      do i=1,nsymm
         do j=1,3
            Write(34,*) (symmop(i,j,k), k=1,3)
         enddo
         Write(34,*) (symmtr(i,j), j=1,3)
      enddo
      Write(34,*) natr
      do i=1,natr
         Write(34,*) nat(i), x(i), y(i), z(i)
      enddo
      Write(34,'(A)') 'END'      
      END

