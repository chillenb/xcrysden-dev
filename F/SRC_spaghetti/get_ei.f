      subroutine get_ei (aline,ener,ne)
c     ***************************************************
c
      IMPLICIT REAL*8 (A-H,O-Z)
cc      INCLUDE 'param.inc'
c
      dimension ener(*)
      character  aline*80
      character  ram_fi*80
c-----------------------------------------------------------------------
c
c.....COPY THE LINE TO AN INTERNAL FILE
      write(ram_fi,'(a80)') aline
c
c.....READ THE EIGENVALUES
      read(ram_fi,'(2x,5f13.7)')  (ener(j), j=ne+1,ne+5)
      if (ener(ne+1).eq. 0.000000)  then
         return
      elseif (ener(ne+2).eq. 0.000000)  then
c        1  eigenvalue has been read
         ne=ne+1
         return
      elseif (ener(ne+3).eq. 0.000000)  then
c        2  eigenvalue have been read
         ne=ne+2
         return
      elseif (ener(ne+4).eq. 0.000000)  then
c        3  eigenvalue have been read
         ne=ne+3
         return
      elseif (ener(ne+5).eq. 0.000000)  then
c        4  eigenvalue have been read
         ne=ne+4
         return
      else
c        5  eigenvalues have been read
         ne=ne+5
         return
      endif
      end
