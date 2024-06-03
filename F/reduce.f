!------------------------------------------------------------------------
!           
!     below are routines for vector reduction that are needed for proper
!     construction of voronoi polyhdra of unreduced lattices
!
!------------------------------------------------------------------------
      
      SUBROUTINE ReduceBasis3d(src,dst)
!        
!     BEWARE: xcrysden usage of vectors in this directory is:
!
!           |11   21   31|           |x1  x2  x3|
!     VEC = |12   22   32| ==> VEC = |y1  y2  y3|
!           |13   23   33|           |z1  z2  z3|
!
! but the "reduce" routines uses vectors in standard sense, i.e.:
!
!     vecA = src(*,1), vecB = src(*,2), vecC = src(*,3)
!
! which is transposed wrt xcrysden usage, hence we will need
! to call MatTranspose 2x
!
      implicit none
      REAL*8 src(3,3), dst(3,3), red(3,3)
      LOGICAL is_basis_reduced
      
!     transpose to vecA = src(*,1), ...
      call MatTranspose(src,red,3,3) ! vec = src^T
      
      if ( .not. is_basis_reduced(red) ) then
         call reduce_basis_doit(red)
      endif
      
!     transpose back to vecA = src(1,*), ...
      call MatTranspose(red,dst,3,3) ! dst = red^T
      
      return
      END

      SUBROUTINE ReduceBasis2d(src,dst)
!
!     BEWARE: xcrysden usage of vectors in this directory is: vecA = src(1,*)
!     but the "reduce" routines uses vectors in standard sense, i.e.: vecA = src(*,1), ...
!     hence we need to call MatTranspose 2x
      !
      implicit none
      REAL*8 src(3,3), dst(3,3), red(3,3)

!     transpose to vecA = src(*,1), ...
      call MatTranspose(src,red,3,3) ! vec = src^T
      
      call gauss_reduce(red(1,1), red(1,2)) ! vecs A, B
      
!     transpose back to vecA = src(1,*), ...
      call MatTranspose(red,dst,3,3) ! dst = red^T
      
      return
      END
      
      SUBROUTINE reduce_basis_doit(vec)
!     BEWARE: vec needs to be in the following form:
!     A = src(*,1), B = src(*,2), C = src(*,3)
      implicit none
      REAL*8 vec(3,3)
      LOGICAL reduced, are_dot_products_reduced
      INTEGER ic, i, j

      ic = 0
      reduced = .false.
      
      do while ( .not. reduced)
         call gauss_reduce(vec(1,2), vec(1,3)) ! vecs B, C
         call gauss_reduce(vec(1,1), vec(1,3)) ! vecs A, C
         call gauss_reduce(vec(1,1), vec(1,2)) ! vecs A, B

         reduced = are_dot_products_reduced(vec)

         ic = ic + 1
         if (ic .gt. 1000) goto 111 
      enddo

 111  continue      
      return
      END

      SUBROUTINE gauss_reduce(v1,v2)
      implicit none
      real*8 v1(3), v2(3), swap(3), VecSize, ScalarProduct, r
      integer i, ir
      logical reduced
      !
      reduced = .false.
      !
      do while ( .not. reduced )
         
         if ( VecSize(v1) - VecSize(v2) .gt. 1d-6 )
     $        call anti_swap_v3(v1,v2) ! v1 must be shorter, anti-swap
             
         r = ScalarProduct(v1,v2) / ScalarProduct(v1,v1)
         ir = nint(r)
         
         if (ir .ne. 0) then
            v2(1) = v2(1) - ir*v1(1)
            v2(2) = v2(2) - ir*v1(2)
            v2(3) = v2(3) - ir*v1(3)
         endif

         reduced = ( VecSize(v2) - VecSize(v1) .gt. -1d-6 )
      enddo
      
      return
      END

      SUBROUTINE anti_swap_v3(v1,v2)
      implicit none
      real*8 v1(3), v2(3), swap(3)
      integer i
      do i=1,3
         swap(i) = -v1(i)
         v1(i)   = v2(i)
         v2(i)   = swap(i)
      enddo
      return
      END
      
      SUBROUTINE sort_basis_vecs(v)
      implicit none
      REAL*8 v(3,3), VecSize
!     BEWARE: "v" needs to be in the following form:
!     A = v(*,1), B = v(*,2), C = v(*,3)

      if ( VecSize(v(1,1)) - VecSize(v(1,2)) .gt. 1d-6 )
     $     call anti_swap_v3(v(1,1), v(1,2))

      if ( VecSize(v(1,1)) - VecSize(v(1,3)) .gt. 1d-6 )
     $     call anti_swap_v3(v(1,1), v(1,3))
      
      if ( VecSize(v(1,2)) - VecSize(v(1,3)) .gt. 1d-6 )
     $     call anti_swap_v3(v(1,2), v(1,3))      
     
      return
      END
       
      logical function is_basis_reduced(v)
      implicit none
!     BEWARE: "v" needs to be in the following form:
!     A = v(*,1), B = v(*,2), C = v(*,3)
      real*8 v(3,3), VecSize
      logical are_dot_products_reduced
            
!     check that |a| <= |b|
      if ( (VecSize(v(1,1)) - VecSize(v(1,2))) .gt. 1d-6 ) then
         is_basis_reduced=.false.
         return
      endif

!     dot-products seem to be a sufficient condition ...
      is_basis_reduced = are_dot_products_reduced(v)
      
      return
      end

      logical function are_dot_products_reduced(v)
      implicit none
      real*8 v(3,3), ScalarProduct
!     BEWARE: "v" needs to be in the following form:
!     A = v(*,1), B = v(*,2), C = v(*,3)
    
      logical ab, ac, bc

      ab = ( dabs(ScalarProduct(v(1,1),v(1,2))) - 0.5d0*
     $     ScalarProduct(v(1,1),v(1,1)) .le. 1d-6 )

      ac = ( dabs(ScalarProduct(v(1,1),v(1,3))) - 0.5d0*
     $     ScalarProduct(v(1,1),v(1,1)) .le. 1d-6 )

      bc = ( dabs(ScalarProduct(v(1,2),v(1,3))) - 0.5d0*
     $     ScalarProduct(v(1,2),v(1,2)) .le. 1d-6 )

      are_dot_products_reduced = ( ab .and. ac .and. bc )
      
      return
      END
