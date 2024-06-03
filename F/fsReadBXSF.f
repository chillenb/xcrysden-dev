c
c     reads BXSF file (BEGIN_BLOCK_BANDGRID3D section) and write
c     ordinary XSF file with lattice vectors and BZ
c
      program fsReadBXSF
      implicit none
c
c     Usage: fsreadBXSF input output
c
      character line*80, input*256, output*256
      logical   finished
      REAL*8    invvec(3,3), dir(3,3), vec(3,3), o1, o2, o3
      integer   string_length, i, j, k, iband

      finished=.false.
      
      if (iargc().ne.2) then
         write(*,*) 'Usage:  fsreadBXSF  input output'
         stop
      endif
      call getarg(1,input)
      call getarg(2,output)
      open(unit=11, file=input,  status='old')
      open(unit=12, file=output, status='unknown')

      do while (.not.finished)
         read(11,'(a20)') line
         i = string_length(line)
         if ( line(1:11) .eq. 'BANDGRID_3D'
     $        .or. line(1:17) .eq. 'BEGIN_BANDGRID_3D' ) then
            read(11,*) iband
            read(11,*) i, j, k
            read(11,*) o1, o2, o3
            read(11,*) ((invvec(i,j),j=1,3),i=1,3)
            finished=.true.
         endif
      enddo

      call ZeroMat       (dir, 3, 3)
      call InvertSqMat123 (invvec, dir, 3)
      call MatTranspose  (dir, vec, 3, 3)

      write(12,*) 'DIM-GROUP'
      write(12,*) '3 1'
      call WVec ('PRIMVEC', vec)
      call WVec ('RECIP-PRIMVEC', invvec)      
      call WignerSeitz(invvec,  'BRILLOUIN-ZONE-PRIMCELL')
      
      write(12,*) 'PRIMCOORD'
      write(12,*) '1 1'
      write(12,*) '1   0.0  0.0  0.0'

      END
