      program ReadBaKgen
      implicit real*8 (a-h,o-z)
!     =========================================
!     Usage: $0 OUTPUTBAND-file OUTPUTKGEN-file
!     =========================================
      PARAMETER ( MINE=1,MAXE=2)
!    $     MAX_NKPT    = 500000,
!    $     MAX_IRNKPT  = 7000,
!    $     MAX_NBANDS  = 700,
!    $     MINE        = 1,
!    $     MAXE        = 2)

      character*6   match_line1
      character*13  match_line2
      character*256 filename
      character*80  line
      character*12  label, label1
      CHARACTER*11  status, form
      logical       not_match
      integer       n(3) 
      integer       index(6),itk
!      integer       coor(3,MAX_NKPT), k_ired_ind(MAX_NKPT)      
!      integer       fileseq_ired_ind(MAX_IRNKPT) 
!      integer       seq2ired_ind(MAX_NKPT)
!      dimension     band_w(MAX_NBANDS,2), eigen(MAX_NBANDS,MAX_IRNKPT)
      integer,allocatable ::       coor(:,:), k_ired_ind(:)      
      integer,allocatable ::       fileseq_ired_ind(:) 
      integer,allocatable ::       seq2ired_ind(:)
      real*8,allocatable ::     band_w(:,:), eigen(:,:)
      dimension     recVec(3,3)

!                        123456
      data match_line1 /'    G1'/
!                        1234567890123
      data match_line2 /'  NO. OF MESH'/

      call getarg(2,filename)
      if(filename.eq.'      ') call getarg(1,filename)
      OPEN(11,FILE=filename,STATUS='OLD',ERR=8000)
 8003 READ(11,*,END=8001) IUNIT,FILENAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FILENAME,STATUS=STATUS,FORM=FORM,ERR=8002)
      GOTO 8003
 8000 WRITE(*,*) ' ERROR IN OPENING BAKGEN.DEF !!!!'
      STOP 'SPAG.DEF'
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FILENAME,'  STATUS: ',STATUS,
     &            '  FORM:',FORM
      STOP 'OPEN FAILED'
 8001 CONTINUE
      label ='            '
      label1='            '
      do i=1,12
         if(filename(i:i).ne.'.') then
            label(i:i)=filename(i:i)
         else
            goto 100
         endif
      enddo

 100  continue

!
!     read IRREDUBIBLE SET of K-POINT  ENERGIES
!   
      read(7,'(a)') line  ! 'IRREDUCIBLE-KPOINT-SET'
      read(7,*) emin, emax
      read(7,*) nbands, n_ir_kpt
      allocate ( coor(3,n_ir_kpt*48), k_ired_ind(n_ir_kpt*48) )     
      allocate ( fileseq_ired_ind(n_ir_kpt) )
      allocate ( seq2ired_ind(n_ir_kpt*48))
      allocate ( band_w(nbands,2), eigen(nbands,n_ir_kpt))
      do i=1,nbands
         read(7,*) iband, band_w(i,MINE), band_w(i,MAXE)
         do j=1,n_ir_kpt
            read(7,*) ik, eigen(i,j)
!            write(*,*) ik, eigen(i,j)
         enddo
      enddo
!      read(7,*) line ! 'END-IRREDUCIBLE-KPOINT-SET'
      close(7)
!     *****************************
!     OUTPUTKGEN file is closed !!!
!     *****************************

!     read reciprocal vectors
!           |11   21   31|           |x1  x2  x3|
!     VEC = |12   22   32|  => VEC = |y1  y2  y3|
!           |13   23   33|           |z1  z2  z3|
!     -------------------------------
!     VEC(i,j) = VEC(#vec,xyz) !!!!!!
!     -------------------------------
      not_match=.true.
      do while(not_match)
         read(8,'(a80)') line
         if (match_line1 .eq. line(1:6)) then
            not_match=.false.
!     recvec are written in columns in OUTPUTKGEN file !!!
            do j=1,3
               read(8,*) (recVec(i,j),i=1,3)
            enddo
         endif
      enddo

!     find the number of K-points
      not_match=.true.
      do while(not_match)
         read(8,'(a80)') line
         if (match_line2 .eq. line(1:13))then
            not_match=.false.
            read(line,'(43x,i7)',err=111) nkpt
            goto 112
 111        read(line,'(44x,i6)') nkpt     ! old format
 112        continue           
         endif
      enddo
!     read the division factors
      read(8,'(53x,3i5)') (n(i),i=1,3)
      read(8,*) line
!      
!     read the k-coordinates(integer) & ireducible-point-indeces
!
      ind=0
      do i=1,n(1)+1
         do j=1,n(2)+1
            do k=1,n(3)+1
               ind=ind+1
               read(8,*) itk,(coor(l,ind),l=1,3),k_ired_ind(ind)
            enddo
         enddo
      enddo

!     read "  weights of k-points:" line !!!
      read(8,*) line
      
      do i=1,n_ir_kpt
         read(8,*) itk, fileseq_ired_ind(i)
      enddo
      close(8)
!     ****************************
!     OUTPUTKGEN file is closed!!!
!     ****************************

!
!     reassign k_ired_ind() as is done in "$WIENROOT/SRC_kgen/zuord.f"
!      
      ind1=0
      ind2=0
      do i=1,n(1)+1
         do j=1,n(2)+1
            do k=1,n(3)+1
               ind1=ind1+1
               if (ind1 .eq. k_ired_ind(ind1)) then
                  ind2=ind2+1
!                  k_ired_ind(ind1)=ind2
                  seq2ired_ind(ind1)=ind2 !!!!
!                  print *,k_ired_ind(ind1),
!     $                 (coor(l,k_ired_ind(ind1)),l=1,3)
               else
                  seq2ired_ind(ind1) = k_ired_ind( k_ired_ind(ind1) )
               endif
            enddo
         enddo
      enddo

!
!     write ENERGIES of K-POINTS in XSF file
!   
      write(10,*) 'BEGIN_BLOCK_BANDGRID3D'
      write(10,*) 'band_energies'
      write(10,*) 'BANDGRID_3D_BANDS'
      write(10,*) nbands
      write(10,*) n(1)+1,n(2)+1,n(3)+1
      write(10,*) 0.0,0.0,0.0  !origin
!     -------------------------------
!     VEC(i,j) = VEC(#vec,xyz) !!!!!!
!     -------------------------------
      write(10,*) (recVec(1,j),j=1,3)    !reciprocal vectors
      write(10,*) (recVec(2,j),j=1,3)    !reciprocal vectors
      write(10,*) (recVec(3,j),j=1,3)    !reciprocal vectors
      do ib=1,nbands
         write(10,*) 'BAND: ',ib
         ic=0
         ind1=0
         do i=1,n(1)+1
            do j=1,n(2)+1
               do k=1,n(3)+1
                  ic=ic+1
                  ind1 = ind1 + 1
                  ind2 = k_ired_ind(ind1)
                  index(ic) = seq2ired_ind(ind2)
!                  index(ic) = ind2
!                  print *,ind1,index(ic),
!                  do l=1,n_ir_kpt
!                     if(fileseq_ired_ind(l) .eq. ind2) index(ic)=l
!                  enddo
                  if(ic.eq.6)then
!                     print *,(index(ii),ii=1,6)
                  write(10,'(6(E12.6,1x))') (eigen(ib,index(ii)),ii=1,6)
                     ic=0                     
                  endif
               enddo
            enddo
         enddo
         if (ic.gt.0) write(10,'(6(E12.6,1x))') 
     &                 (eigen(ib,index(ii)),ii=1,ic)
      enddo
      write(10,*) 'END_BANDGRID_3D'
      write(10,*) 'END_BLOCK_BANDGRID3D'
      goto 1000

 999  print *,'could not open file: ',filename
      stop 'OPEN FAILED'

 1000 continue
      END


!
!     THIS IS THE RIGHT ORDER, STUDY IT !!!
!
!      if(iswitch.eq.1)write(66,*)' internal and cartesian k-vectors:'
!      DO 10 I1=1,N(1)+1                                                 
!      DO 10 I2=1,N(2)+1                                                 
!      DO 10 I3=1,N(3)+1                                                 
!      I=I+1                                                             
!      IF(I.GT.NMSHP) THEN                                               
!        PRINT*,'I.GT.NMSHP,STOP',I,I1,I2,I3                             
!        STOP                                                            
!      END IF                                                            
!      IF(I.EQ.NUM(I))THEN                                               
!        NDIM=NDIM+1                                                     
!        IF(NDIM.GT.IDKP) THEN                                           
!          WRITE(66,1000)                                                
!1000      FORMAT(1H ,'NUMBER OF INEQUIVALENT POINTS ECCEEDS IDKP')      
!          STOP                                                          
!        END IF                                                          
!        NUM(I)=NDIM                                                     
!        RINDA=(DBLE(I1-1)+DBLE(ISHIFT(1))/2.D0)/DBLE(N(1))              
!        RINDB=(DBLE(I2-1)+DBLE(ISHIFT(2))/2.D0)/DBLE(N(2))              
!        RINDC=(DBLE(I3-1)+DBLE(ISHIFT(3))/2.D0)/DBLE(N(3))              
!        BK(1,NDIM)=GBAS(1,1)*RINDA+GBAS(2,1)*RINDB+GBAS(3,1)*RINDC      
!        BK(2,NDIM)=GBAS(1,2)*RINDA+GBAS(2,2)*RINDB+GBAS(3,2)*RINDC      
!        BK(3,NDIM)=GBAS(1,3)*RINDA+GBAS(2,3)*RINDB+GBAS(3,3)*RINDC      
!c       BK(1,NDIM)=GBAS(1,1)*RINDA+GBAS(1,2)*RINDB+GBAS(1,3)*RINDC      
!c       BK(2,NDIM)=GBAS(2,1)*RINDA+GBAS(2,2)*RINDB+GBAS(2,3)*RINDC      
!c       BK(3,NDIM)=GBAS(3,1)*RINDA+GBAS(3,2)*RINDB+GBAS(3,3)*RINDC      
!        bki(1,ndim)=rinda
!        bki(2,ndim)=rindb
!        bki(3,ndim)=rindc
!        if(iswitch.eq.1) 
!     *  write(66,100) rinda,rindb,rindc,bk(1,ndim),bk(2,ndim),bk(3,ndim)
! 100     format(3f10.5,10x,3f10.5)
!      ELSE                                                              
!        IF(NUM(I).GT.NMSHP) PRINT*,'ERROR'                              
!        NUM(I)=NUM(NUM(I))
!                                             
!      END IF                                                            
!10    CONTINUE                                                          
!      NKP=NDIM                                                          
!      RETURN                                                            
!      END                                                               
