c     *******************************************
      program InpIni2XCR
      implicit real*8 (a-h,o-z)
c     Usage: inpini2xcr inpini_file [get_species]
c
c     if get_species argument is present (it can be what-ever) 
c     then just return the list of species
c     *******************************************
      include 'param.inc'
      PARAMETER( !type of execution
     $     IGET_FTN34=0,   ! produce ftn34
     $     IGET_SPECIES=1, ! get species names
     $     BOHR=0.529177d0 ! bohr
     $     )

      character*256 file
      character*20
     $     file_out,
     $     name(NAC),
     $     coor_name(NAC)
      integer
     $     nat(NAC),
     $     coor_index(NAC)
      logical
     $     tford(NAC)      
      real*8
     $     a(3,3),
     $     coor(3,NAC)
      
      narg=iargc()
      if (narg.lt.1 .or. narg.gt.2)
     $     stop 'Usage: fhi_coord2xcr <fhi_coord.out file> [getspecies]'
      iexec_type = IGET_FTN34
      if (narg.eq.2) iexec_type = IGET_SPECIES

      call getarg(1,file)
      len = index(file,' ') - 1
      open(unit=1, file=file(1:len), status='old')
      file_out = 'fhi_coord.xcr'
      open(unit=2, file=file_out, status='unknown')

      do i=1,3
         read(1,*) (a(j,i), j=1,3)
         do j=1,3
            a(j,i) = BOHR*a(j,i)
         enddo
      enddo
      read(1,*) n_all_species      
      iat = 0
      do i=1,n_all_species
         read(1,*) n_i_species
         read(1,*) name(i)
         do ii=1,n_i_species
            iat = iat + 1
            read(1,*) (coor(j,iat), j=1,3), tford(iat)
            coor_name(iat)  = name(i)
            coor_index(iat) = i
            do j=1,3
               coor(j,iat) = BOHR * coor(j,iat)
            enddo
         enddo
      enddo

c     ***
c     *** RETURN list of species
c     ***
      if (iexec_type .eq. IGET_SPECIES)then
         open(unit=3,file='fhi_species_name.list', status='unknown')
         write(3,*) n_all_species
         write(3,'(a10)') (name(i),i=1,n_all_species)
      else
c     ***
c     *** WRITE XSF file
c     ***
         open(unit=4, file='fhi_species_nat.list', status='old')
         read(4,*) n_all_species
         read(4,*) (nat(i),i=1,n_all_species)

         write(2,*) 'DIM-GROUP'
         write(2,*) '3 1'
         write(2,*) 'PRIMVEC'
         write(2,1000) a
         write(2,*) 'PRIMCOORD'
         write(2,*) iat, 1
         do i=1,iat
            write(2,1001) nat(coor_index(i)), (coor(j,i), j=1,3)
         enddo
      endif
 1000 format(2(3(F15.9,2X),/),3(F15.9,2X))
 1001 format(I5,3(F15.9,2X))

      END
