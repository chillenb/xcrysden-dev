      SUBROUTINE GetCenteredPoints(dv, dvc, idim, cp, npoi)
      implicit real*8 (a-h,o-z)
      real*8 dv(3,3), dvc(3,3), cp(3,4)
      real*8 xyz2con(3,3),res_dv(3),res_dvc(3),poi(3),poiC(3)
      real*8 parallelepipedDiagonal
      integer SHORTEST_DIAGONAL, LONGEST_DIAGONAL
      parameter (tolmin= 1.0d-7)
      include 'p2c.inc'

c     dv   ... primitive direct vectors
c     dvc  ... conventional reciprocal vectors
c     cp   ... special centered points (calculated & returned)
c     npoi ... # of special centered points (calculated & returned)

c     ********************************************************************
c     calculate the inverse matrix of conventional vectors matrix in order
c     to convert from XYZ coor to conventional fractional coor
c     ********************************************************************
c      call MatCopy(dvc, xyz2con, idim, idim)
      call InvertSqMat123(dvc, xyz2con, idim)


c      print *,'GCP: xyz2con:'
c      do 10 i=1,3
c 10         print *,xyz2con(1,i),xyz2con(2,i),xyz2con(3,i)
c
c      print *,'GCP: dvc:'
c      do 11 i=1,3
c 11      print *,dvc(1,i),dvc(2,i),dvc(3,i)
c
c      print *,'GPC: dv:',dv


c     this is just in case: (polymers & slabs)
c      do i=1,3
c         do j=1,3
c            if ((i.gt.idim) .or. (j.gt.idim)) then
c               if (i.eq.j) then
c                  xyz2con(i,j) = 1.0d0
c               else
c                  xyz2con(i,j) = 0.0d0
c               endif
c            endif
c         enddo
c      enddo

c     *************************************************************************
c     determine sphere size of primitive and conventional vectors and its ratio
c     *************************************************************************
      do i=1,3
c     shortest possible diagonal
         res_dv(i)  = parallelepipedDiagonal ( dv, SHORTEST_DIAGONAL )

c         print *,'res_dv:',res_dv

c     longest possible diagonal
         res_dvc(i) = parallelepipedDiagonal ( dvc, LONGEST_DIAGONAL )
      enddo
      rp = VecSize( res_dv )

c      print *,'rp:',rp
      
c     due to above calculation res_dv can be zero vector and rp can be zero;
c     if that happens set it to half Angstroms
      if (rp .lt. tolmin) rp=1.5d0
      rc = VecSize( res_dvc )
      ir = int( rc/rp )
      if ( rc/rp - dble(int(rc/rp)) .gt. tolmin ) ir = ir + 1

      print *, 'GCP: ', rp,rc,ir

c     **********************
c     get the special points
c     **********************
      irx = ir
      iry = ir
      irz = ir
      if (idim.lt.3) irz=0
      if (idim.lt.2) iry=0
      if (idim.lt.1) then
         npoi    = 1
         cp(1,1) = 0.0d0
         cp(2,1) = 0.0d0
         cp(3,1) = 0.0d0
         RETURN                 ! it is a MOLECULE; return
      endif

      npoi=0      
      do i=-irx,irx
         do j=-iry,iry
            do k=-irx,irx
               di = dble(i)
               dj = dble(j)
               dk = dble(k)
               do l=1,3
                  poi(l) = di*dv(1,l) + dj*dv(2,l) + dk*dv(3,l)
               enddo

c     multiply p with conversion matrix to get conventinal fractional coor
c               call MatMult(xyz2con, 3, 3, poi, 1, poiC)
               do ii=1,3
                  poiC(ii) = xyz2con(1,ii)*poi(1) +
     $                 xyz2con(2,ii)*poi(2) + xyz2con(3,ii)*poi(3)
               enddo

c     print *, 'GCP: poi=', poi(1), poi(2), poi(3)
c     print *, 'GCP: fract_poiC=', poiC(1), poiC(2), poiC(3)

c     check if poiC is inside the cell: condition 0<=poiC<1
               if ( poiC(1).gt.-tolmin .and. poiC(1).lt.(1.d0-tolmin)
     $              .and.
     $              poiC(2).gt.-tolmin .and. poiC(2).lt.(1.d0-tolmin)
     $              .and.
     $              poiC(3).gt.-tolmin .and. poiC(3).lt.(1.d0-tolmin) )
     $              then
                  npoi = npoi + 1
                  do l=1,3
                     cp(l,npoi) = poiC(l)
                  enddo
c                  print *, 'Centered point: ',
c     $                 cp(1,npoi), cp(2,npoi), cp(3,npoi)
c                  print *, 'Centered point: ',
c     $                 poiC(1), poiC(2), poiC(3)
               endif

            enddo
         enddo
      enddo

      RETURN
      END


      REAL*8 FUNCTION parallelepipedDiagonal(vec,type)
      implicit none
      real*8  vec(3,3), sq_diag(4), diag_d(3,4), min, max
      integer i, type, SHORTEST_DIAGONAL, LONGEST_DIAGONAL

      include 'p2c.inc'

c      diag_d(i,1) = +vec(i,1) + vec(i,2) + vec(i,3)
c      diag_d(i,3) = -vec(i,1) + vec(i,2) + vec(i,3) 
c      diag_d(i,4) = +vec(i,1) - vec(i,2) + vec(i,3)
c      diag_d(i,2) = -vec(i,1) - vec(i,2) + vec(i,3) 

      do i=1,3
c     diagonal-1: (0,0,0)-->(1,1,1)
         diag_d(i,1) = vec(i,1) + vec(i,2) + vec(i,3)
         
c     diagonal-2: (1,1,0)-->(0,0,1)
         diag_d(i,2) = vec(i,3) - vec(i,1) - vec(i,2) 
         
c     diagonal-3: (1,0,0)-->(0,1,1)
         diag_d(i,3) = vec(i,2) + vec(i,3) - vec(i,1)
         
c     diagonal-4: (0,1,0)-->(1,0,1)
         diag_d(i,4) = vec(i,1) + vec(i,3) - vec(i,2)
      enddo

      do i=1,4
         sq_diag(i) =
     $        diag_d(1,i)*diag_d(1,i) +
     $        diag_d(2,i)*diag_d(2,i) + 
     $        diag_d(3,i)*diag_d(3,i)
      enddo

      if (type .eq. SHORTEST_DIAGONAL) then
         min = sq_diag(1)
         if (sq_diag(2) .lt. min) min = sq_diag(2)
         if (sq_diag(3) .lt. min) min = sq_diag(3)
         if (sq_diag(4) .lt. min) min = sq_diag(4)

         parallelepipedDiagonal = sqrt(min)
      else
         max = sq_diag(1)
         if (sq_diag(2) .gt. max) max = sq_diag(2)
         if (sq_diag(3) .gt. max) max = sq_diag(3)
         if (sq_diag(4) .gt. max) max = sq_diag(4)
         
         parallelepipedDiagonal = sqrt(max)
      endif

      return
      end
