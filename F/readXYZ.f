      SUBROUTINE ReadXYZ(NATR,NAT,X,Y,Z)
      include 'param.inc'
      INTEGER NAT(NAC)
      REAL*8 X(NAC),Y(NAC),Z(NAC)
      CHARACTER*80 DUMMY
      CHARACTER*4 ATSYM
      CHARACTER*10 ATOM(100), ATOM_UPPER(100)

      include 'atoms.inc'

      do while (.true.)
         read(33,*,end=200) natr
         read(33,*) dummy

         do i=1,natr
            read(33,*) atsym, x(i), y(i), z(i)
c     from "ATOMIC SYMBOL" to "ATOMIC NUMBER"
            len0=Index(atsym,' ')-1
            do j=1,100
               len=Index(atom(j),' ')-1
               if ( (atsym(1:len0).eq.atom(j)(1:len))
     $              .or. (atsym(1:len0).eq.atom_upper(j)(1:len)) ) then
                  nat(i)=j
                  goto 100
               endif
            enddo
 100        continue
         enddo
      enddo

 200  continue

      RETURN
      END
