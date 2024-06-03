      program MultiSlab

c     *************************************************************
c     program reads SLAB structure from xc_struct.$$ and generate a 
c     multislab
c     *************************************************************
c
c     USAGE: multislab <mode> <file1> [< <vacuum_thickness_file>]
c
c     mode:  SYMMINFO              
c            MULTISLAB_NO_INVERSION
c            MULTISLAB_INVERSION   
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'param.inc'
      INTEGER NAT(NAC), INN(NAC)
      INTEGER MAXNAT ! maximum atomic number in structure
      INTEGER NPT, ISPLIT ! used for Wien97
      REAL*8 X(NAC),Y(NAC),Z(NAC),DVEC(3,3),IDVEC(3,3),TMPVEC(3,3)
      REAL*8 A(NAC),B(NAC),C(NAC), MIN_NNDIST(NAC), RMT(NAC), R0(NAC)
      REAL*8 SYMMOP(1,3,3), SYMMTR(1,3)
      REAL*8 V1(3), V2(3), V3(3)
      CHARACTER*80 FILE1, TITLE, MODE
      CHARACTER*10 ATOM(100), ATOM_UPPER(100)
      CHARACTER*4 type
      REAL*8 FRMT(NAC)
      LOGICAL l_inv_center, l_inv, lzero3
      COMMON/MULTAT/ X,Y,Z,NAT,DVEC,NATR
      COMMON/IVEC/ IDVEC
      COMMON/BOHRR/ BOHR,IMODE3 !this is not used, but is needed because 
                                !this program uses readf1 routine
      PARAMETER (RAD2DEG = 57.295779513084d0, TOLMIN=1.d-7)
      PARAMETER (
     $     M1_SYMMINFO = 0,
     $     M1_MULTISLAB_NO_INVERSION = 1,
     $     M1_MULTISLAB_INVERSION = 2 )

      DATA SYMMOP/
     $     1.0d0, 0.0d0, 0.0d0,
     $     0.0d0, 1.0d0, 0.0d0,
     $     0.0d0, 0.0d0, 1.0d0/
      DATA SYMMTR/
     $     0.0d0, 0.0d0, 0.0d0/
     
      include 'atoms.inc'

C     CONVERSION FROM BOHRS TO ANGS.
      BOHR=0.529177d0

      narg=IArgc()
      if (narg.ne.2) then
         print *,
     $'Usage: multislab <mode> <file1> [< <vacuum_thickness_file>]'
         stop
      endif

      Call GetArg(1,MODE)
      len = index(mode, ' ') - 1
      if (mode(1:len) .eq. 'SYMMINFO') imode=M1_SYMMINFO
      if (mode(1:len) .eq. 'MULTISLAB_NO_INVERSION')
     $     imode=M1_MULTISLAB_NO_INVERSION
      if (mode(1:len) .eq. 'MULTISLAB_INVERSION')
     $     imode=M1_MULTISLAB_INVERSION
      
      Call GetArg(2,FILE1)
      Open(unit=11,file=file1,status='old')

      Call ReadF1('DIM-GROUP', igroup, idim)
      if (idim.lt.2) stop 'MultiSlab:: so far I am converting to
     $     multislab just from slabs and multi-slabs (dim=2,3) !!!'

      Call ReadF1('PRIMVEC', 0, idim)
      Call ReadF1('PRIMCOORD', 0, idim)
      Close(11)

      if(imode.gt.M1_SYMMINFO) then
         Read(*,*) vac_thickness
      else
         vac_thickness = 0.0d0
      endif

c     find an atom with lower Z and set Z -> 0.0
      zmin=z(1)
      zmax=z(1)
      do i=2,natr
         if (z(i).lt.zmin) zmin=z(i)
         if (z(i).gt.zmax) zmax=z(i)
      enddo

      do i=1,natr
         z(i)=z(i)-zmin
      enddo

c     assign the third vector
      dvec(3,1)=0.0d0
      dvec(3,2)=0.0d0
      dvec(3,3)=zmax-zmin
      lzero3=.false.
      if (abs(dvec(3,3)).lt.tolmin) then
         dvec(3,3) = 1.0d0
         lzero3=.true.
      endif

c     get the inverse matrix
      call Invert3x3(dvec,idvec)

c     ****************************************************************
c     now get the fractional coordinates of atoms and 
c     translate them if necessary ( in WIEN97 fractional coorinates 
c     are in range [0,1) )
      maxnat = 1
      do i=1,natr
         nat(i) = nat(i) - int( nat(i) / 100 ) * 100 ! nat from 0 to 99
         if(nat(i).gt.maxnat) maxnat = nat(i)
         a(i) = x(i)*idvec(1,1) + y(i)*idvec(2,1) + z(i)*idvec(3,1)
         b(i) = x(i)*idvec(1,2) + y(i)*idvec(2,2) + z(i)*idvec(3,2)
         c(i) = x(i)*idvec(1,3) + y(i)*idvec(2,3) + z(i)*idvec(3,3)
         if(lzero3) c(i)=0.0d0
         if(a(i).lt.0.0d0) a(i) = a(i) + 1.0d0
         if(b(i).lt.0.0d0) b(i) = b(i) + 1.0d0
         if(c(i).lt.0.0d0) c(i) = c(i) + 1.0d0
c        Write(6,21) nat(i), a(i), b(i), c(i)
      enddo

c    ////////////////////////////////////
c     do we have center of inversion !!! 
c
c     THIS DOESN'T WORK AT ALL --> fix in the future !!!
c    ////////////////////////////////////  
      l_inv_center=.true.
      i=1
      do while (i.lt.natr .and. l_inv_center)
         l_inv=.false.
         j=i+1
         do while (j.le.natr .and. .not.l_inv)
            if ( abs(1.0d0-a(i)-a(j)) .lt. tolmin .and.
     $           abs(1.0d0-b(i)-b(j)) .lt. tolmin .and.
     $           abs(1.0d0-c(i)-c(j)) .lt. tolmin ) l_inv=.true.
            j=j+1
         enddo
         if (.not. l_inv) l_inv_center=.false.
         i=i+1
      enddo

      if(imode.eq.M1_SYMMINFO) then
         if(l_inv_center) then
            write(*,*) 1
         else
            write(*,*) 0
         endif
      else 
         if(imode.eq.M1_MULTISLAB_INVERSION .and. l_inv_center) then
            vact2 = vac_thickness / 2.0d0
            do i=1,natr
               z(i)=z(i) + vact2
            enddo
         endif
c     re-assign the 4th component of third vector
         dvec(3,3)=zmax-zmin + vac_thickness

c     write the C95's unit 34; this will not be used in future
         title = 'EXTERNAL'
         Call WriteFTN34(title, 3, 1, dvec, 1, symmop, symmtr,
     $        natr, nat, x, y, z)
         
c     get the inverse matrix
      call Invert3x3(dvec,idvec)
         
c     ****************************************************************
c     now get the fractional coordinates of atoms and 
c     translate them if necessary ( in WIEN97 fractional coorinates 
c     are in range [0,1) )
         maxnat = 1
         do i=1,natr
            nat(i) = nat(i) - int( nat(i) / 100 ) * 100 ! nat from 0 to 99
            if(nat(i).gt.maxnat) maxnat = nat(i)
            a(i) = x(i)*idvec(1,1) + y(i)*idvec(2,1) + z(i)*idvec(3,1)
            b(i) = x(i)*idvec(1,2) + y(i)*idvec(2,2) + z(i)*idvec(3,2)
            c(i) = x(i)*idvec(1,3) + y(i)*idvec(2,3) + z(i)*idvec(3,3)
            if(a(i).lt.0.0d0) a(i) = a(i) + 1.0d0
            if(b(i).lt.0.0d0) b(i) = b(i) + 1.0d0
            if(c(i).lt.0.0d0) c(i) = c(i) + 1.0d0
c     Write(6,21) nat(i), a(i), b(i), c(i)
         enddo
      
c     get the nearest neighbor distances
         do i=1,natr
            min_nndist(i) = 99.9d10

            do j=1,natr
               if (i.ne.j) then
                  a2 = a(j) + nint(a(i)-a(j))
                  b2 = b(j) + nint(b(i)-b(j))
                  c2 = c(j) + nint(c(i)-c(j))
                  
                  x1 = a(i)*dvec(1,1) + b(i)*dvec(2,1) + c(i)*dvec(3,1)
                  y1 = a(i)*dvec(1,2) + b(i)*dvec(2,2) + c(i)*dvec(3,2)
                  z1 = a(i)*dvec(1,3) + b(i)*dvec(2,3) + c(i)*dvec(3,3)
                  
                  x2 = a2*dvec(1,1) + b2*dvec(2,1) + c2*dvec(3,1)
                  y2 = a2*dvec(1,2) + b2*dvec(2,2) + c2*dvec(3,2)
                  z2 = a2*dvec(1,3) + b2*dvec(2,3) + c2*dvec(3,3)
                  
                  d = Dist(x1,y1,z1, x2,y2,z2)
                  if (d.lt.min_nndist(i)) then
                     min_nndist(i) = d
                     inn(i) = j
                  endif
               endif
            enddo
         enddo
         
c     *************************
c     get a,b,c,alpha,beta,gamma
         ap = VecSize3(dvec,1) / BOHR
         bp = VecSize3(dvec,2) / BOHR
         cp = VecSize3(dvec,3) / BOHR
         do i=1,3
            v1(i) = dvec(1,i)
            v2(i) = dvec(2,i)
            v3(i) = dvec(3,i)
         enddo
         alpha = dacos( ScalarProduct(v2, v3) /
     $        (VecSize(v2) * VecSize(v3)) ) * RAD2DEG
         beta  = dacos( ScalarProduct(v3,v1) /
     $        (VecSize(v3) * VecSize(v1)) ) * RAD2DEG
         gamma = dacos( ScalarProduct(v1,v2) /
     $        (VecSize(v1) * VecSize(v2)) ) * RAD2DEG

         type='P   '
c     if gamma is very close to 120.d0 make it 120.0d0 and set lattice
c     type to hexagonal
         if ( abs(gamma-120.d0) .lt. tolmin ) then
            gamma=120.0d0
            type='H   '
         endif

c     *******************************************
c     get some aproximation for muffin tin radius
c     *******************************************
c     let the smallest offset between RMT spheres be 0.1 a.u.
         OFFSET = 0.1d0/2.d0    ! each sphre gets half of offset
         
         do i=1,natr
            if (nat(i).eq.1) frmt(i) = 0.50d0
            if (nat(i).gt.1  .and. nat(i).lt.21)  frmt(i) = 1.00d0
            if (nat(i).ge.21 .and. nat(i).lt.39) frmt(i) = 1.10d0
            if (nat(i).ge.39 .and. nat(i).lt.57) frmt(i) = 1.15d0
            if (nat(i).ge.57) frmt(i) = 1.20d0
            
            r0(i) = 0.000005d0
            if (nat(i).le.71) r0(i) = 0.00001d0
            if (nat(i).le.36) r0(i) = 0.00005d0
            if (nat(i).le.18) r0(i) = 0.0001d0
         enddo
         
         do i=1,natr
            rmt(i) = frmt(i) * min_nndist(i) /
     $           ( BOHR * (frmt(i) + frmt(inn(i))) ) - OFFSET
         enddo
         isplit = 8
         npt    = 781
         
c     ****************************
c     now write WIEN97 struct file
c     ****************************
         write(*,'(a)')
     $        'Wien97 struct file generated by XCrySDen program'
         write(*,'(a4,a24,i2)') type, 'LATTICE.NONEQUIV. ATOMS:',natr
         if (maxnat.lt.40) then
            write(*,'(a13,a4)') 'MODE OF CALC=', 'NREL'
         else
            write(*,'(a13,a4)') 'MODE OF CALC=', 'RELA'
         endif
         write(*,'(6f10.6)') ap, bp, cp, alpha, beta, gamma
         
         do i=1,natr
            write(*,'(a5,i3,a4,f10.8,a3,f10.8,a3,f10.8)')
     $           'ATOM=', -i,': X=', a(i), ' Y=', b(i), ' Z=', c(i)
            write(*,'(a15,i2,a17,i2)') 'MULT=', 1, 'ISPLIT=', isplit
            write(*,'(a10,a5,i5,a5,f10.8,a5,f10.4,a5,f5.1)')
     $           atom(nat(i)),'NPT=',npt,'R0=',r0(i),
     $           'RMT=',rmt(i),'Z:', real(nat(i))
            write(*,'(a20,3f10.7)') 'LOCAL ROT MATRIX:   ', 1.0,0.0,0.0
            write(*,'(20x,3f10.7)') 0.0, 1.0, 0.0
            write(*,'(20x,3f10.7)') 0.0, 0.0, 1.0
         enddo
         write(*,'(i4,A)') 1,'      NUMBER OF SYMMETRY OPERATIONS'
         write(*,'(3i2,f10.7)') 1, 0, 0, 0.0
         write(*,'(3i2,f10.7)') 0, 1, 0, 0.0
         write(*,'(3i2,f10.7)') 0, 0, 1, 0.0
         write(*,'(i8)') 1
 21      FORMAT(I5,3F16.10)                                       
      endif

      END
      

 
