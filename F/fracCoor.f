c     ------------------------------
      program fracCoor
c     Usage: fracCoor XSF-file x y z
c     ------------------------------
      implicit real*8 (a-h, o-z)
      include 'param.inc'

      CHARACTER filen*120, ax*32, ay*32, az*32
      REAL*8 xx(3), x(NAC), y(NAC), z(NAC), dvec(3,3), rvec(3,3), v(3)
      INTEGER nat(NAC)

      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR

      narg=iargc()
      if ( narg .ne. 4 ) then
         print *,'Usage: fracCoor XSF-file x y z'
         STOP
      endif

      call getarg(1,filen)
      open(unit=11, file=filen, status='old')

      call getarg(2,ax)
      call getarg(3,ay)
      call getarg(4,az)
      read(ax,*) xx(1)
      read(ay,*) xx(2)
      read(az,*) xx(3)
      
      call readf1('DIM-GROUP', igroup, idim)
      call readf1('PRIMVEC', idum1, idim)

      if (idim .eq. 0) then
         write(6,'(3f15.8,2x)') xx(1), xx(2), xx(3)
         goto 999
      else
         call InvertSqMat123(dvec, rvec, idim)
c         print *,dvec
c         print *,rvec
         if (idim .eq. 1) then
            write(6,'(3f15.8,2x)') xx(1)*rvec(1,1), xx(2), xx(3)
         else if (idim .eq. 2) then
            a = rvec(1,1)*xx(1) + rvec(2,1)*xx(2)
            b = rvec(1,2)*xx(1) + rvec(2,2)*xx(2)
            write(6,'(3f15.8,2x)') a, b, xx(3)
         else if (idim .eq. 3) then
            a = rvec(1,1)*xx(1) + rvec(2,1)*xx(2) + rvec(3,1)*xx(3)
            b = rvec(1,2)*xx(1) + rvec(2,2)*xx(2) + rvec(3,2)*xx(3)
            c = rvec(1,3)*xx(1) + rvec(2,3)*xx(2) + rvec(3,3)*xx(3)
            write(6,'(3f15.8,2x)') a, b, c
         endif
      endif
      
 999  continue
      END

c      subroutine ZeroVec(v,ind)
c      REAL*8 v(3,3)
c      do i=1,3
c         v(ind,i)=0.0d0
c      enddo
c      return
c      END
c
c      subroutine UnitVec(v,ind)
c      REAL*8 v(3,3)
c      do i=1,3
c         v(ind,i)=0.0d0
c      enddo
c      v(ind,ind)=1.0d0
c      return
c      END
