      program recv
      character*10 idum
      real*8 a(3,3),ar(3,3)
      read(34,*) IDUM
      read(34,*) l,k
      read(34,*) ((a(i,j),j=1,3),i=1,3)
      write(*,'(3(3f12.5/))') ((a(i,j),j=1,3),i=1,3)
      call recvec(a,ar)
      write(*,'(3(3f12.5/))') ((ar(i,j),j=1,3),i=1,3)
      end


C     VECTOR1 GOES TO FIRST COLOUMN, VETCOR2 TO SECOND & VECTOR3 TO 
C     THIRD COLOUMN;;; DV(COL,ROW)
C           |11   21   31|          |x1  x2  x3|
C     VEC = |12   22   32|  => DV = |y1  y2  y3|
C           |13   23   33|          |z1  z2  z3|
      SUBROUTINE RecVec(a,ar)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 a(3,3),ar(3,3)
      REAL*8 ab(3),bc(3),ca(3),PI2

      PI2=8.0d0*DATAN(1.0d0) 

      do i=1,3
         do j=1,3
            a(i,j)=a(i,j)/0.529d0
c            print *,a(i,j)
         enddo
      enddo
      CALL VecProduct(a(1,2),a(1,3),bc)
      CALL VecProduct(a(1,3),a(1,1),ca)
      CALL VecProduct(a(1,1),a(1,2),ab)

c      write(*,'(3(3f12.5/))') (a(i),i=1,3),(bv(i),i=1,3),(cv(i),i=1,3)

      vol=a(1,1)*bc(1) + a(2,1)*bc(2) + a(3,1)*bc(3)
      do i=1,3
         ar(1,i)=PI2*bc(i)/vol
         ar(2,i)=PI2*ca(i)/vol
         ar(3,i)=PI2*ab(i)/vol
      enddo
      RETURN
      END
