
C     *********************************************
C     *** program based on spag.f program of WIENxx
C     *********************************************


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Program modified for the present purposes by:                         C
C     ------                                                                C
C     Anton Kokalj                               Email: Tone.Kokalj@ijs.si  C
C     Dept. of Physical and Organic Chemistry    Phone: x 386 1 477 3520    C
C     Jozef Stefan Institute                       Fax: x 386 1 477 3811    C
C     Jamova 39, SI-1000 Ljubljana                                          C
C     SLOVENIA                                                              C
C                                                                           C
C     Source: $TESTS/wnReadBands.f                                          C
C     ------                                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     case spaghetti:
c     set exe = spaghetti  
c     cat << theend > $def
c     5    , '$file.insp',                   'old','formatted',0
c     6    , '$file.outputsp$updn',      'unknown','formatted',0
c     9    , '$file.qtl$updn',           'unknown','formatted',0
c     10   ,'$file.spaghetti${updn}_ene','unknown','formatted',0
c     11   ,'$file.spaghetti${updn}_ps', 'unknown','formatted',0
c     20   ,'$file.struct',                  'old','formatted',0
c     7    ,'$file.output1${updn}',          'old','formatted',0 
c     theend

      program ReadBands
      IMPLICIT REAL*8 (A-H,O-Z)

ccc      include 'SRC_spaghetti/param.inc'
      PARAMETER  (MINE = 1, MAXE = 2)
c     
      character  aline*80,fname*80,ram_file*80
      CHARACTER*11      STATUS,FORM
      character  k_pat*7,ei_pat*28
      character  dummy*1,label1*12
      character  k_name*12,symbol*12,label*12
      character  title*80,lattice*4
c     
ccc      dimension  jatom_list(NATO)

ccc      dimension  eigen(NEVL,NKP),char(NEVL,NKP)
ccc      dimension  n_ene(NKP)
ccc      dimension  vk(3,NKP),k_name(NKP)
ccc      dimension  lines(NKP),xval(NKP)
ccc      dimension  band_w(NEVL,2)

      real*8,allocatable  ::  eigen(:,:)    
      integer,allocatable ::  n_ene(:)
      dimension  vk(3)     
      real*8,allocatable  ::  band_w(:,:)
      real*8 tmp_eigen(9999)

      dimension  qtl(13),icomma(15)
ccc      logical    break(nkp)
c     
      data  k_pat  /'     K='/
      data  ei_pat /'EIGENVALUES BELOW THE ENERGY'/

c     
      call getarg(2,fname)
      if(fname.eq.'      ') call getarg(1,fname)
      OPEN(1,FILE=fname,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,
     *     ERR=8002)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING BAND.DEF !!!!'
      STOP 'BAND.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS,
     *     '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
      label='            '
      label1='            '
      do i=1,12
         if(fname(i:i).ne.'.') then
            label(i:i)=fname(i:i)
         else
            goto 9
         endif
      enddo

c     
c.....READ K-VECTORS
c
 9    n_kpt=0
      nn_ene_max=0
  10   CONTINUE
      read(7,'(a80)',end=20) aline
      if  (aline(1:7).eq.k_pat)  then
c     THIS IS THE BEGINNING OF AN EIGENVALUE SECTION (output1)
         n_kpt=n_kpt + 1
         nn_ene=0
         read(7,'(a1)') dummy
         read(7,'(a1)') dummy
c     READ EIGENVALUES OF K-VECTOR
 11     read(7,'(a80)') aline
         if (aline(15:42).eq.ei_pat)  then
            goto 10
         else
            call get_ei(aline, tmp_eigen, nn_ene)
            if(nn_ene_max.lt.nn_ene) nn_ene_max=nn_ene
            goto 11
         endif
      else
         goto 10
      endif
 20   continue
      allocate(band_w(nn_ene_max,2),eigen(nn_ene_max,n_kpt))
      allocate(n_ene(n_kpt))

      rewind 7
      n_kpt=0
 100  CONTINUE
      read(7,'(a80)',end=200) aline
      if  (aline(1:7).eq.k_pat)  then
c     THIS IS THE BEGINNING OF AN EIGENVALUE SECTION (output1)
         n_kpt=n_kpt + 1
cc         if(n_kpt.gt.nkp) then
cc            write(6,*) ' parameter NKP too small:',nkp
cc            stop 'NKP TOO SMALL'
cc         end if
         n_ene(n_kpt)=0
         call get_k (aline,vk(1),vk(2),
     &        vk(3),k_name)
         read(7,'(a1)') dummy
         read(7,'(a1)') dummy
c     READ EIGENVALUES OF K-VECTOR
 110     read(7,'(a80)') aline
         if (aline(15:42).eq.ei_pat)  then
c     write(6,*) 'k=',k_name(n_kpt),'  done'
c     write(6,*) 'numb e(k) found=',n_ene(n_kpt)
            goto 100
         else
cc            if(n_ene(n_kpt).ge.nevl-4) then
cc               write(6,*) 'parameter nevl too small:',nevl
cc               stop 'NEVL TOO SMALL'
cc            endif
            call get_ei(aline, eigen(1,n_kpt), n_ene(n_kpt))
            goto 110
         endif
      else
         goto 100
      endif
c     
c.....ALL K-VECTORS HAVE BEEN READ; SEARCH FOR K-POINT WITH SMALLEST
c     NUMBER OF EIGENVALUES
c     
 200  continue
c      write(*,*) 'number of k-points read=',n_kpt
cc      nu_min=NEVL+1
      nu_min=9999999
      do 205 j=1,n_kpt
         if (n_ene(j).lt.nu_min)  then
            nu_min=n_ene(j)
            k_min=j
         endif
 205  continue
c      write(6,*) 'smallest number eigenvalues at k=',k_min,' (',
c     *     k_name(k_min),')'
c      write(6,*) '         =',nu_min

c     tk
c     DETERMINE BAND WIDTHS in terms of (Emin,Emax)
c
c     insert code here !!!
      emin=+99999.9
      emax=-99999.9
      do i=1,nu_min
         band_w(i,MINE)=+99999.9
         band_w(i,MAXE)=-99999.9
         do j=1,n_kpt
            if(eigen(i,j).lt.band_w(i,MINE)) band_w(i,MINE)=eigen(i,j)
            if(eigen(i,j).gt.band_w(i,MAXE)) band_w(i,MAXE)=eigen(i,j)
         enddo
         if(band_w(i,MINE).lt.emin) emin=band_w(i,MINE)
         if(band_w(i,MAXE).gt.emax) emax=band_w(i,MAXE)
      enddo

c
c     write IRREDUBIBLE SET of K-POINTS and ENERGIES
c   
      write(6,*) 'IRREDUCIBLE-KPOINT-SET'
      write(6,'(2F12.6)') emin, emax
      write(6,*) nu_min, n_kpt
      write(8,*) nu_min
      do i=1,nu_min
         write(6,'(10x,I6,10x,2F12.6)')
     $        i, band_w(i,MINE), band_w(i,MAXE)
         write(8,'(10x,I6,10x,2F12.6)')
     $        i, band_w(i,MINE), band_w(i,MAXE)
         do j=1,n_kpt
c            write(*,'(i6,2x,3F10.6,5x,F12.6)')
c     $           j, vk(1,j), vk(2,j), vk(3,j), eigen(i,j)
            write(6,'(i6,2x,F12.6)') j, eigen(i,j)
         enddo
      enddo
      write(6,*) 'END-IRREDUCIBLE-KPOINT-SET'

      END
