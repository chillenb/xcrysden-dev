c     *****************************************************************
      program cube2xsf_fabris
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

      implicit real*8 (a-h, o-z)

      PARAMETER (
     $     M_SCA=1,
     $     M_SCA_VEC=4,
     $     M_SCA_VEC_SCA=5,
     $     B2A=0.529177d0,
     $     MAXATOM = 5000,
     $     MAXMO   = 150, ! maximum number of orbitals
     $     MAXGRID = 75   ! maximum size of the grid, i.e., (MAGRID x MAXGRID x MAXGRID)
     $     )

      character filename*256, title*80, mode*1, fileout*256
      integer n(3),             !increments
     $     mo_index(MAXMO),      !indices of MOs
     $     nat(MAXATOM)

      logical lmo               !do we have MO plot

      real*8
     $     xo(3),               !origin
     $     dx(3,3),             !increments
     $     v(3,3),              !spanning vectors
     $     val(MAXGRID,MAXGRID,MAXGRID,MAXMO),  !values in cubefile
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

      read(1,*) title
      read(1,*) title
      read(1,*) natoms, (xo(i),i=1,3)
      do i=1,3
         read(1,*) n(i), (dx(i,j),j=1,3)
      enddo

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

c     read all the values and store in memory
      do i1=1,n(1)
         do i2=1,n(2)
            read(1,'(6(E13.5):)') ((val(i1,i2,i3,j),j=1,nmo),i3=1,n(3))
         enddo
      enddo

c
c     loop over all the orbitals
c
      do imo=1,nmo

         if (lmo) then
            write(fileout,'(''mo-'',i3.3,''.xsf'')') mo_index(imo)
         else
            write(fileout,'(''denpot.xsf'')')
         endif

         open(unit=9, file=fileout, status='unknown')
                  
         write(9,'(1x,a5)') 'ATOMS'
         do i=1,natoms
            write(9,'(i3,2x,3f12.6)') nat(i), x(i), y(i), z(i)
         enddo
         
         write(9,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
         if(lmo)then 
            write(9,'(a16,i3.3)') 'g98_3D_ORBITAL_#',mo_index(imo)
         else
            write(9,'(a)') 'g98_3D_unknown'
         endif
         write(9,'(a)') 'DATAGRID_3D_g98Cube'
         write(9,*) n(3),n(2),n(1)
         write(9,*) (xo(i)*B2A,i=1,3)
         do i=3,1,-1
            write(9,*) (v(i,j)*B2A,j=1,3)
         enddo

c     write values to XSF files
         
         do i1=1,n(1)
            do i2=1,n(2)
               write(9,'(6(E13.5):)') (val(i1,i2,i3,imo),i3=1,n(3))
            enddo
         enddo

c     write XFS footer
         
         write(9,'(a)') 'END_DATAGRID_3D'
         write(9,'(a)') 'END_BLOCK_DATAGRID_3D'
      enddo

      END



