c     ********************************************************
      program InpIni2XCR
c     Usage: inpini2xcr inpini_file [get_species]
c
c     if get_species argument is present (it can be what-ever) 
c     then just return the list of species
c     ********************************************************
      include 'param.inc'
      PARAMETER( !type of execution
     $     IGET_FTN34=0,        ! produce ftn34
     $     IGET_SPECIES=1       ! get species names
     $     )                    

      implicit real*8 (a-h,o-z)
      character*256 file
      character*10
     $     name(NAC)
      integer
     $     nrx(3),
     $     ineq_pos(3),
     $     symop(48,3,3),
     $     coor_index(NAC),
     $     nat(NAC)
      real*8
     $     nel,
     $     xk(3,10000),
     $     wkpt(10000),
     $     a(3,3),
     $     b(3,3),
     $     valence_charge(NAC),
     $     ion_fac,
     $     ion_damp,
     $     coor(3,NAC),
     $     vau0(3,NAC),
     $     coor_valcharge(NAC)
      logical
     $     tmetal, tdegen,
     $     tmold, tband,
     $     tpsmesh, coordwave,
     $     t_init_basis_s, t_init_basis_p, t_init_basis_d,
     $     coorflag(3,NAC), tford(NAC)
      
      narg = iargc()
      if (narg.lt.1 .or. narg.gt.2)
     $     stop 'Usage: fhi_inpini2xcr <inpini_file> [get_species]'
      iexec_type = IGET_FTN34
      if (narg.eq.2) iexec_type = IGET_SPECIES

      call getarg(1,file)
      len = index(file,' ') - 1
      open(unit=1, file=file(1:len), status='old')

      read(1,*) nsx,nax,nx,ngwx1,ngwx88
      read(1,*) ngwix,nx_init,nrx(1),nrx(2),nrx(3),nschltz
      read(1,*) nx_basis,max_basis_n,nlmax_init
      read(1,*) nnrx,nkpt,nlmax,mmaxx,n_fft_store
      read(1,*) minpes, ngrpx, nrpes
      read(1,*) ibrav, pgind
      read(1,*) nel,tmetal,ekt,tdegen
      read(1,*) ecut,ecuti
      read(1,*) tmold,tband,nrho
      read(1,*) npos, nthm, nseed
      read(1,*) T_ion, T_init, Q, nfi_rescale
      read(1,*) nsp,tpsmesh,coordwave
      read(1,*) nkpt
      do i=1,nkpt
         read(1,*) (xk(j,i), j=1,3), wkpt(i)
      enddo
      do i=1,3
         read(1,*) (a(j,i), j=1,3)
      enddo
      do i=1,3
         read(1,*) (b(j,i), j=1,3)
      enddo
      read(1,*) alat,omega      

      iat=0
      do i=1,nsx
         read(1,*) name(i),number,valence_charge(i), ion_fac
         read(1,*) ion_damp,rgauss,l_max,l_loc
         read(1,*) t_init_basis_s, t_init_basis_p, t_init_basis_d
         do ii=1,number
            iat = iat + 1
            read(1,*) (coor(j,iat), j=1,3),
     $           (coorflag(j,iat), j=1,3), tford(iat)
            if (npos.eq.1 .or. npos.eq.4) read(1,*) (vau0(j,iat), j=1,3)
            coor_index(iat)     = i
            coor_valcharge(iat) = valence_charge(i)
         enddo
         read(1,*) (ineq_pos(j),j=1,3)
      enddo
c     **********************************************************************
c     we already have enough information, but lets read the rest of the file
c     **********************************************************************
      read(1,*) nsymop
      do ii=1,nsymop
         read(1,'(i3)') isymop
         do i=1,3
            read(1,*) (symop(isymop,j,i),j=1,3)
         enddo
      enddo
c     after SYMOP there are two lines - I don't know the meaning - IGNORE

c     ***
c     *** RETURN list of species
c     ***
      if (iexec_type .eq. IGET_SPECIES)then
         open(unit=2,file='fhi_species_name.list', status='unknown')
         write(2,*) nsx
         write(2,'(a10)') (name(i),i=1,nsx)
      else
c     ***
c     *** WRITE fort.34 file
c     ***
         open(unit=3, file='fhi_species_nat.list', status='old')
         read(3,*) nsx
         read(3,*) (nat(i),i=1,nsx)
         igroup=1
         if ( ibrav.eq.0 ) igroup=1 !user-supplied --> primitive (SO FAR!!!)
         if ( ibrav.eq.1 ) igroup=1 !simple-cubic  --> primitive
         if ( ibrav.eq.2 ) igroup=5 !fcc           --> fcc
         if ( ibrav.eq.3 ) igroup=6 !bcc           --> bcc
         if ( ibrav.eq.4 ) igroup=8 !haxagonal     --> hexagonal
         if ( ibrav.eq.8 ) igroup=1 !orthorhombic  --> primitive
         if ( ibrav.eq.10) igroup=7 !rhombohedral  --> rhombohedral
         if ( ibrav.eq.12) igroup=1 !bc orthorombic--> primitive
         write(34,*) 'EXTERNAL'
         write(34,*) 3,igroup
         write(34,1000) a
         write(34,*) '1'
         write(34,*) '1. 0. 0.'
         write(34,*) '0. 1. 0.'
         write(34,*) '0. 0. 1.'
         write(34,*) '0. 0. 0.'
         write(34,*) iat
         do i=1,iat
            write(34,*) nat(coor_index(i)), (coor(j,i), j=1,3)
         enddo
         write(34,*) 'END'
      endif
 1000 format(2(3(F15.9,2X),/),3(F15.9,2X))

      END
