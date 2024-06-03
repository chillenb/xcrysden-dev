c --------------------------------------------------------------------------
c Peter BLAHA, Inst.f.Techn.Elektrochemie, TU Vienna, A-1060 Vienna
c Phone: +43-1-58801-5187             FAX: +43-1-5868937
c Email: pblaha@email.tuwien.ac.at    WWW: http://www.tuwien.ac.at/theochem/
c --------------------------------------------------------------------------

c     ******************************************************************
c     Changes made by Tone Kokalj, Date: Mon Jul 13 15:36:23 CEST 1998
c
c     1. filename is now specified as command line argument
c     2. conversion from fractional coordinates to Cartesian was changed
c
c     where change was made -> I signed as t.k
c     ******************************************************************

      program str2xcr
C     
      implicit real*8 (a-h,o-z)
      parameter (nato=100000)
      character*4 typ
      dimension x(4,nato),klz(nato)
      dimension con(3,3),bra(3,3) !t.k
      logical rhombo !t.k
C     x.....Atomkoordinaten, Kugelradius
C     klz...Kernladungszahl      
c     *
c     * rhombo = .true. when system is rhombohedral;
c     *          because constants are specified as conventional (hexagonal),
c     *          but fractional coordinates are specified as 
c     *          primitive (rhombohedral)
c     *

      rhombo=.false. !t.k
      call o_file
C     oeffnet die noetigen files (fehlerkontrolle)!                
      call matrix(nat,typ,a,b,c,bra,con,rhombo) !t.k
C     schreibt bravais-matrix, uebergibt nat (number of atoms)
      call koord(nat,nato,na,x,klz)
      call writeStr(na,x,klz,typ,a,b,c,bra,con,rhombo) !t.k
      call option
      end
      

      subroutine o_file
      character*256 namstr,namlat
      character*10 test
      call GetArg(1,namstr) !t.k
      do 1000 i=256,1,-1
1000  if(namstr(i:i).ne.' ') goto 1010
1010  if(i.gt.110)goto 9992
      namlat=namstr      
      namstr(i+1:i+7)='.struct'
      namlat(i+1:i+4)='.xcr'
      open(unit=1,file=namstr,form='formatted',status='old',err=9993)
      open(unit=2,file=namlat,form='formatted',
     $     status='unknown',err=9994)
      return
9992  write(*,*)'Filename too long' 
      stop      
9993  write(*,*)'struct-File does not exist'      
      stop
9994  write(*,*)'File ',namlat(1:i+4),' already exists.'
      write(*,*)' Overwrite (y/n)?'
      read(*,99901)test
      if(test.eq.'y' .or. test.eq.'Y') then
       open(unit=2,file=namlat,form='formatted')
       return
      else
       stop
      endif    
99900 format(a80) 
99901 format(a1)     
      end
      

      subroutine matrix(nat,typ,a,b,c,bra,con,rhombo) !t.k
      implicit real*8 (a-h,o-z)      
      dimension bra(3,3),con(3,3) !t.k
      character*4 typ
      logical rhombo !t.k
      character*79 title
c     *
c     * con(3,3) is matrix that holds the conventional vectors
c     *
      do 10 i=1,3
         do 10 j=1,3
            con(j,i)=0.0d0      !t.k
            bra(j,i)=0.0d0
 10   continue
   
      read(1,99900) title
      write(2,*) 'INFO' !t.k
      write(2,*) title
      write(2,*) 'END_INFO'
      read(1,99901) typ,nat
      read(1,99902)
      read(1,99903) a,b,c,alpha,beta,gamma
      if(alpha.eq.0d0) alpha=90d0
      if(beta .eq.0d0) beta =90d0
      if(gamma.eq.0d0) gamma=90d0

      pi    = 4.d0*atan(1.d0)
      gamma1= gamma*pi/180.d0
      beta1 = beta*pi/180.d0
      alpha1= alpha*pi/180.d0
      cosg1 =(cos(gamma1)-cos(alpha1)*cos(beta1))/sin(alpha1)/sin(beta1)
      gamma0=acos(cosg1)
c      gamma=gamma*acos(0d0)/90d0

cIGROUP:	
c       1.....P centered
c	2.....A centered
c	3.....B centered (not body, but B-face)
c	4.....C centered
c	5.....F centered
c	6.....I centered (body centered)
c	7.....R centered (rhomohedral)
c	8.....H centered (hexagonal)   (space groups: 168-194)
c	9.....TRIGONAL_NOT_RHOMBOHEDRAL (space groups: 143-167)
c     umrechnung ins bogenmass      

c     this is ussualy used; if not it will be changed when needed
      con(1,1)=a                !t.k
      con(2,2)=b                !t.k
      con(3,3)=c                !t.k
      if(typ(1:3).eq.'CXY') then
         bra(1,1)=a/2d0
         bra(2,1)=-b/2d0
         bra(1,2)=a/2d0
         bra(2,2)=b/2d0
         bra(3,3)=c       
         igroup=4
      else if (typ(1:3).eq.'CYZ') then
         bra(1,1)=a
         bra(2,2)=-b/2d0
         bra(3,2)=c/2
         bra(2,3)=b/2d0
         bra(3,3)=c/2d0
         igroup=2
      else if(typ(1:3).eq.'CXZ') then
         con(1,1)=a*sin(gamma1)  !t.k
         con(2,1)=a*cos(gamma1)  !t.k
         bra(1,1)=a*sin(gamma1)/2d0
         bra(2,1)=a*cos(gamma1)/2d0
         bra(3,1)=-c/2d0
         bra(2,2)=b
         bra(1,3)=a*sin(gamma1)/2d0
         bra(2,3)=a*cos(gamma1)/2d0
         bra(3,3)=c/2d0
         igroup=3
      else if(typ(1:1).eq.'P') then
         bra(1,1)=a*sin(gamma0)*sin(beta1)                
         bra(2,1)=a*cos(gamma0)*sin(beta1)                 
         bra(3,1)=a*cos(beta1)                                   
         bra(1,2)=0.0d0                                                      
         bra(2,2)=b*sin(alpha1)
         bra(3,2)=b*cos(alpha1)
         bra(1,3)=0.0d0                                                      
         bra(2,3)=0.0d0                                                      
         bra(3,3)=c
         do i=1,3
            do j=1,3
               con(i,j)=bra(i,j) !conventional vectors equal to primitive
            enddo
         enddo
c     *** OLD (no triclini lattice) ***
c         con(1,1)=a*sin(gamma)  !t.k
c         con(2,1)=a*cos(gamma)  !t.k
c         bra(1,1)=a*sin(gamma)
c         bra(2,1)=a*cos(gamma)
c         bra(2,2)=b
c         bra(3,3)=c
         igroup=1
      else if(typ(1:1).eq.'F') then
         bra(1,1)=a/2d0
         bra(2,1)=b/2d0
         bra(1,2)=a/2d0
         bra(3,2)=c/2d0
         bra(2,3)=b/2d0
         bra(3,3)=c/2d0
         igroup=5
      else if(typ(1:1).eq.'B') then
c     Florent suggested this:
         bra(1,1)=-a/2d0
         bra(2,1)=b/2d0
         bra(3,1)=c/2d0
         bra(1,2)=a/2d0
         bra(2,2)=-b/2d0
         bra(3,2)=c/2d0
         bra(1,3)=a/2d0
         bra(2,3)=b/2d0
         bra(3,3)=-c/2d0
         igroup=6
c         bra(1,1)=a/2d0
c         bra(2,1)=-b/2d0
c         bra(3,1)=c/2d0
c         bra(1,2)=a/2d0
c         bra(2,2)=b/2d0
c         bra(3,2)=-c/2d0
c         bra(1,3)=-a/2d0
c         bra(2,3)=b/2d0
c         bra(3,3)=c/2d0
c         igroup=6
      else if(typ(1:1).eq.'R') then
         rhombo=.true.          !t.k
c     conventional vectors are hexagonal
         con(1,1)=a*sqrt(3d0)/2d0
         con(2,1)=-a/2d0
         bra(1,1)=a/sqrt(3d0)/2d0
         bra(2,1)=-a/2d0
         bra(3,1)=c/3d0
         bra(1,2)=a/sqrt(3d0)/2d0
         bra(2,2)=a/2d0
         bra(3,2)=c/3d0
         bra(1,3)=-a/sqrt(3d0)
         bra(3,3)=c/3d0
         igroup=7
      else if(typ(1:1).eq.'H') then
c     conventional and primitive vectors are the same in this case
         con(1,1)=a*sqrt(3d0)/2d0 !t.k
         con(2,1)=-a/2d0        !t.k
         bra(1,1)=a*sqrt(3d0)/2d0
         bra(2,1)=-a/2d0
         bra(2,2)=a
         bra(3,3)=c
         igroup=8
      else 
         goto 9990
      end if
      write(2,*) 'DIM-GROUP'
      write(2,*) 3,igroup
      write(2,*) 'PRIMVEC'
      write(2,99904) bra
      write(2,*) 'CONVVEC'
      write(2,99904) con
      write(2,*) 'SYMMOP'
      write(2,*) 1
      write(2,*) 1.d0,0.d0,0.d0
      write(2,*) 0.d0,1.d0,0.d0
      write(2,*) 0.d0,0.d0,1.d0
      write(2,*) 0.d0,0.d0,0.d0
      return
 9990 write(*,*)'Unknown lattice type'
      stop      
99900 format(a79)      
99901 format(a4,23x,i3)      
99902 format(13x,a4)
99903 format(6F10.7)
99904 format(2(3F15.9/),3f15.9)
      end      


      subroutine koord(nat,nato,na,x,klz)
      implicit real*8 (a-h,o-z)
      dimension x(4,*),klz(*)
      na=0
      do 10 i=1,nat
       na=na+1
       if(na.gt.nato) goto 9990
       read(1,99990)x(1,na),x(2,na),x(3,na)
       read(1,99991)mult
       do 11 im=2,mult
        na=na+1
        if(na.gt.nato) goto 9990
11      read(1,99990)x(1,na),x(2,na),x(3,na)        
       read(1,99992)rmt,z
       do 12 im=na-mult+1,na
        x(4,im)=rmt/3d0
12      klz(im)=z
       read(1,99993)
       read(1,99993)
10     read(1,99993)
      return      
9990  write(*,*)'Too many atoms. Increase parameter nato in source ',
     c'code.'
      stop
99990 format(5x,3x,4x,f10.7,3x,f10.7,3x,f10.7)
99991 format(15x,i2,17x,i2)
99992 format(10x,5x,5x,5x,10x,5x,f10.5,5x,f5.2)
99993 format(20x,3f10.8)
      end


      subroutine writeStr(na,x,klz,typ,a,b,c,bra,con,rhombo) !t.k
      implicit real*8 (a-h,o-z)
      real*8 ibra(3,3) !t.k
      character*4 typ
      dimension x(4,*),klz(*),bra(3,3),con(3,3),tmp(3,3)
      logical rhombo !t.k

      write(2,*) 'PRIMCOORD'
      write(2,99900) na, 1

c     t.k: changed completely!!!
      do i=1,na
         if ( .not. rhombo ) then
c     from conventional fracional coordinates to Cartesian
            xx=x(1,i)*con(1,1) + x(2,i)*con(1,2) + x(3,i)*con(1,3)
            yy=x(1,i)*con(2,1) + x(2,i)*con(2,2) + x(3,i)*con(2,3)
            zz=x(1,i)*con(3,1) + x(2,i)*con(3,2) + x(3,i)*con(3,3)

c     *********************************************************
c     some obla-di obla-da manipulation will be neccessary !!!
c     *********************************************************
c     AIM: primitive fractional coordinates must be 
c     between [-0.5,0.5);
c     we must get the corresponding Cartesin coordinates
c     *********************************************************

c     initialise ibra & tmp matrix; tmp is dummy anyway
            call Invert3x3(bra,ibra)
c     goto primitive fractional coordinates
            aa=xx*ibra(1,1) + yy*ibra(1,2) + zz*ibra(1,3)         
            bb=xx*ibra(2,1) + yy*ibra(2,2) + zz*ibra(2,3)
            cc=xx*ibra(3,1) + yy*ibra(3,2) + zz*ibra(3,3)
         else
c     rhobohedral fractional coordinates are already specified in 
c     primitive (rhombohedral) coordinates
            aa=x(1,i)
            bb=x(2,i)
            cc=x(3,i)
         endif
c     aa,bb,cc must be in range of (-1,1);
         aa=aa-dnint(aa)         !is this good enough!!!!
         bb=bb-dnint(bb)
         cc=cc-dnint(cc)
c     go back to Cartesian coordinate
         xx=aa*bra(1,1) + bb*bra(1,2) + cc*bra(1,3)         
         yy=aa*bra(2,1) + bb*bra(2,2) + cc*bra(2,3)
         zz=aa*bra(3,1) + bb*bra(3,2) + cc*bra(3,3)
         write(2,99901)klz(i),xx,yy,zz
      enddo
c         write(2,99901)x(1,i),x(2,i),x(3,i),x(4,i),klz(i)
      return
99900 format(2i5)
99901 format(i4,3f15.9)
      end
      

      subroutine option
      implicit real*8 (a-h,o-z)
      close(2)
      write(*,*)'Translation has been finished.'
      return
      end      


