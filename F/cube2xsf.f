c     *****************************************************************
      program cube2xsf      
c     Usage: cube2xsf mode cubefile
c
c     MODE:
c          1 ... scalar property (like density)
c          4 ... scalar+vector   (like density+gradient)
c          5 ... scalar+vector+scalar (like density+gradient+laplacian)
c
c     BEWARE: onle MODE.eq.1 is implemented !!!
c
c     AUTHOR:  Anton Kokalj (tone.kokalj@ijs.si)
c     DATE:    Sat Feb 22 15:36:05 CET 2003
c     *****************************************************************
      implicit none
      integer M_SCA, M_SCA_VEC, M_SCA_VEC_SCA, MAXATOM, MAXMO, MAXGRID

      real*8 B2A

      PARAMETER (
     $     M_SCA=1,
     $     M_SCA_VEC=4,
     $     M_SCA_VEC_SCA=5,
     $     B2A=0.529177d0,
     $     MAXATOM = 5000,
     $     MAXMO   = 1000,
     $     MAXGRID = 1000
     $     )

      character filename*256, title*80, mode*1, fileout*256
      integer n(3),             !increments
     $     mo_index(MAXMO),      !indices of MOs
     $     nat(MAXATOM),
     $     natoms, narg, i, i1, i2, i3, j, imo, nmo, imode

      logical lmo               !do we have MO plot

      real*8 c,                 ! charge
     $     xo(3),               !origin
     $     dx(3,3),             !increments
     $     v(3,3),              !spanning vectors
     $     val(MAXGRID,MAXMO),  !values in cubefile
     $     x(MAXATOM), y(MAXATOM), z(MAXATOM)
      
      lmo=.false.

      narg = iargc()
      if (narg.ne.2)
     $     stop 'Usage: cube2xsf mode cubefile'

      call getarg(1,mode)
      call getarg(2,filename)
      read(mode,*) imode
      if (imode.ne.1) then
         write(*,*) 'imode.gt.1 not implemented'
         STOP
      endif
      open(unit=1, file=filename(1:len(filename)), status='old')

      read(1,'(a80)') title
      read(1,'(a80)') title
      read(1,*) natoms, (xo(i),i=1,3)
      do i=1,3
         read(1,*) n(i), (dx(i,j),j=1,3)
      enddo

      if ( n(3) .gt. MAXGRID )
     $     STOP 'NZ > MAXGRID, increase MAXGRID parameter'

c     calculate the spanning vectors
      do i=1,3
         do j=1,3
c     matrix array is as: v(#,xyz)
            v(i,j) = dble(n(i)-1)*dx(i,j)
         enddo
      enddo

c     **
c     ** do we have MO plot ???
c     **
      if ( natoms .lt. 0 ) then
         lmo=.true.
         natoms=-1*natoms
      endif

c     **
c     ** fastest runing direction in G98 is 1,2,3 and the corrsponding in
c     ** XSF is 3,2,1; CORRECT THAT
      do i=1,natoms
         read(1,*) nat(i), c, x(i), y(i), z(i)
         x(i) = x(i) * B2A
         y(i) = y(i) * B2A
         z(i) = z(i) * B2A
c         write(*,'(i3,2x,3f12.6)') nat, x*B2A, y*B2A, z*B2A
      enddo

      nmo=1
      if (lmo) then
         read(1,'(10I5:)') nmo, (mo_index(i),i=1,nmo)
      else
         mo_index(1)=1
      endif

      if ( nmo .gt. MAXMO )
     $     STOP 'NMO > MAXMO, increase MAXMO parameter'

c     **
c     ** write XSF DATAGRID header 
c     ** 
      do imo=1,nmo
         if (lmo) then
            write(fileout,'(''mo-'',i3.3,''.xsf'')') mo_index(imo)
         else
            write(fileout,'(''denpot.xsf'')')
         endif

         open(unit=9+imo, file=fileout, status='unknown')
         write(9+imo,'(1x,a5)') 'ATOMS'
         do i=1,natoms
            write(9+imo,'(i3,2x,3f12.6)') nat(i), x(i), y(i), z(i)
         enddo
         
         write(9+imo,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
         if(lmo)then 
            write(9+imo,'(a16,i3.3)') 'g98_3D_ORBITAL_#',mo_index(imo)
         else
            write(9+imo,'(a)') 'g98_3D_unknown'
         endif
         write(9+imo,'(a)') 'DATAGRID_3D_g98Cube'
         write(9+imo,*) n(3),n(2),n(1)
         write(9+imo,'(3e20.12)') (xo(i)*B2A,i=1,3)
         do i=3,1,-1
            write(9+imo,'(3e20.12)') (v(i,j)*B2A,j=1,3)
         enddo
      enddo

c     write values to XSF files

      do i1=1,n(1)
         do i2=1,n(2)
c            read(1,'(6(E13.5):)') ((val(i3,j),j=1,nmo),i3=1,n(3))

c
c     Michael Rutter requested for the free format reading of cube
c     file. Let's see how this will work.
c

            read(1,err=99,fmt=*) ((val(i3,j),j=1,nmo),i3=1,n(3))

            do imo=1,nmo
               write(9+imo,'(6(E13.5):)') (val(i3,imo),i3=1,n(3))
            enddo
         enddo
      enddo

 100  continue
c     write XFS footer

      do imo=1,nmo
         write(9+imo,'(a)') 'END_DATAGRID_3D'
         write(9+imo,'(a)') 'END_BLOCK_DATAGRID_3D'
      enddo
      goto 999
      
 99   continue
c     error occured during reading
      write(*,*) 'Error reading CUBE datagrid: i1 =', i1, '; i2 =', i2
      write(*,*) '             current i3 and imo = ', i3, imo
      goto 100
      
 999  continue
      END



