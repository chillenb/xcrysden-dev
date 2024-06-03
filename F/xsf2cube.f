c     *****************************************************************
      program xsf2cube
c     Usage: xsf2cube xsffile
c
c     BEWARE: the cube file should contain solely a single scalar field !!!!
c
c     AUTHOR:  Anton Kokalj (tone.kokalj@ijs.si)
c     DATE:    Tue Jan 2011
c     *****************************************************************

      implicit none

      integer
     $     maxat,
     $     nxx,
     $     nyy,
     $     nzz
      real*8 B2A

      parameter (
     $     maxat = 10000,
     $     nxx   = 500,
     $     nyy   = 500,
     $     nzz   = 800,
     $     B2A=0.529177d0
     $     )

      integer i, j, k, natoms, nx, ny, nz, atn(maxat)
      integer xsf_unit, out_unit, count
      character dummy*256, xsfFile*256
      character*4 nat(maxat)
      real value(nxx,nyy,nzz), charge(maxat), val
      real*8 vec(3,3), tau(3,maxat), origin(3), span(3,3)
      character*2 atom(100), atom_upper(100), aval*25
      include 'atoms.inc'

      call getarg(1,xsfFile)

      out_unit=6

      if ( (len_trim(xsfFile).eq.1) .and. xsfFile(1:1).eq.'-' ) then
c     read XSF from standard
         xsf_unit=5
      else
         xsf_unit=20
         open(unit=xsf_unit,
     $        file=xsfFile(1:len_trim(xsfFile)), status='old')
      endif

c     CRYSTAL|SLAB|POLYMER|MOLECULE
      read  (xsf_unit,'(a80)') dummy
      

c     PRIMVEC
      read  (xsf_unit,'(a80)') dummy
      read (xsf_unit,*) ((vec(i,j),i=1,3),j=1,3)
      
c     PRIMCOORD
      read  (xsf_unit,'(a80)') dummy
      read (xsf_unit,*) natoms
      do i=1,natoms
         read(xsf_unit,*) nat(i), (tau(j,i),j=1,3)
         do j=1,100
            if ( (nat(i)(1:min(len_trim(nat(i)),2)).eq.
     $           atom(j)(1:len_trim(atom(j)))) .or. 
     $           (nat(i)(1:min(len_trim(nat(i)),2)).eq.
     $           atom_upper(j)(1:len_trim(atom(j)))) ) then
               atn(i)=j
               charge(i)=real(j)
            endif
         enddo
      enddo

c     DATAGRID header
      read (xsf_unit,*) dummy
      read (xsf_unit,*) dummy
      read (xsf_unit,*) dummy

      read (xsf_unit,*) nx, ny, nz
      read (xsf_unit,*) (origin(i),i=1,3)
      read (xsf_unit,*) ((span(i,j),i=1,3),j=1,3)

      if ( nx.gt.nxx) stop 'NX > NXX, increase the NXX parameter'
      if ( ny.gt.nyy) stop 'NY > NYY, increase the NYY parameter'
      if ( nz.gt.nzz) stop 'NZ > NZZ, increase the NZZ parameter'

c     read the datagrid

      read (xsf_unit,*)
     $     (((value(i,j,k),i=1,nx),j=1,ny),k=1,nz)

c     Cube header

      write(out_unit,*) 'scalar_field'
      write(out_unit,*) 'converted from XSF file'
      write(out_unit,'(I5,3F12.6)') natoms, (origin(i)/B2A,i=1,3)

c     N.B.: to make a proper cube file for Bader program, the cube file
c     should not contain the edge points, which are periorid replicas,
c     hence the actual number of points is (nx-1,ny-1,nz-1) !!!
      
      write(out_unit,'(I5,3F12.6)')
     $     nx-1, (span(i,1)/B2A/dble(nx-1), i=1,3)
      write(out_unit,'(I5,3F12.6)')
     $     ny-1, (span(i,2)/B2A/dble(ny-1), i=1,3)
      write(out_unit,'(I5,3F12.6)')
     $     nz-1, (span(i,3)/B2A/dble(nz-1), i=1,3)
      do i=1,natoms
         write(out_unit,'(I5,4F12.6)')
     $        atn(i), charge(i), (tau(j,i)/B2A,j=1,3)
      enddo

c     write the datagrid

      do i=1,nx-1
         do j=1,ny-1
            write (out_unit,'(6(E13.5):)') (value(i,j,k),k=1,nz-1)
         enddo
      enddo

      END



      
