      SUBROUTINE ConvexHull(
     *     index,        !index of plane/polygon
     *     iverplane,    !data of vertex collection for index's plane
     *     nvrt,         !number of vertices in index's plane     
     *     allvertex)    !all vertex coordinates

      IMPLICIT NONE
      include 'voronoi.inc'

      INTEGER index
      INTEGER iverplane(MAX_PL_VER,MAX_POLY)
      INTEGER nvrt
      REAL*8  allvertex(3,MAX_VERT)

      REAL*8  
     *     vertex(3,MAX_PL_VER), !vertex of polygon
     *     vec1(3),vec2(3),nml(3), !nml == normal of plane/polygon
     *     normal(4),           !normal of polygon
     *     x(3),y(3),z(3)       !arbitrarily esed
      INTEGER iorder(MAX_NORD)  !order of sorted vertices
      LOGICAL lfirst,notend
      REAL*8 cosyfi, sinyfi, coszfi, sinzfi, minfi, fi, yfi, zfi,
     $     ver1, xnml1, PolarAngle, sign
      INTEGER i, j, iTmp, iYmax, iYmin, nord

c     We have plane in 3D; it maybe parallel to Y axis or near parallel and
c     that is something that we not want, we will rotate the plane such that 
c     it will be parallel to XY plane, and than we will get the order of 
c     vertices. If we have the order of vertices, we are able do draw a convex
c     poligon.

c     take three points, and calculate normal of plane
      do i=1,3
         vec1(i) = allvertex(i,iverplane(3,index)) -
     *        allvertex(i,iverplane(1,index))
         vec2(i) = allvertex(i,iverplane(2,index)) -
     *        allvertex(i,iverplane(1,index))         
      enddo

      CALL VecProduct(vec1,vec2,nml)
      
c     --- --- DEBUG_BEGIN --- ---
c      PRINT *,'VEC1:: ',(vec1(i),i=1,3)
c      PRINT *,'VEC2:: ',(vec2(i),i=1,3)
c      PRINT *,'NML :: ',(nml(i), i=1,3)
c      PRINT *,'nmlsize = ',VecSize(nml)
c      do i=1,nvrt
c         PRINT *,iverplane(i,index),';   ',
c     *              (allvertex(j,iverplane(i,index)),j=1,3)
c      enddo
c     --- --- DEBUG_END --- ---

c     maybe nml is already Z-oriented
      if(abs(nml(1)).lt.MAXERR2 .and. abs(nml(2)).lt.MAXERR2) then
         do i=1,nvrt
            vertex(1,i) = allvertex(1,iverplane(i,index))
            vertex(2,i) = allvertex(2,iverplane(i,index))
            vertex(3,i) = allvertex(3,iverplane(i,index))         
         enddo
         goto 99
      endif

c     first rotate around Z axis
      zfi = dacos(nml(1) / dsqrt(nml(1)*nml(1) + nml(2)*nml(2)))
      if(nml(2).gt.0.0d0) zfi=-zfi
      coszfi=dcos(zfi)
      sinzfi=dsin(zfi)
      xnml1  = nml(1)*coszfi - nml(2)*sinzfi
      nml(2) = nml(1)*sinzfi + nml(2)*coszfi
      nml(1) = xnml1            !this must be so
c     now rotate around Y-axis
      yfi = dacos(nml(3) / dsqrt(nml(1)*nml(1) + nml(2)*nml(2) + 
     *     nml(3)*nml(3)))
      if(nml(1).gt.0.0d0) yfi=-yfi
      cosyfi=dcos(yfi)
      sinyfi=dsin(yfi)

c     --- --- DEBUG_BEGIN --- ---
      nml(1) =  xnml1*cosyfi + nml(3)*sinyfi
      nml(3) = -xnml1*sinyfi + nml(3)*sinyfi

c      PRINT *,'zfi     =',zfi,   ';  yfi     =',yfi
c      PRINT *,'cos(zfi)=',coszfi,';  cos(yfi)=',cosyfi
c      PRINT *,'sin(zfi)=',sinzfi,';  sin(yfi)=',sinyfi
c      PRINT *,'nml>', nml(1),nml(2),nml(3)
      if(abs(nml(1)).gt.1.d-5 .or. abs(nml(2)).gt.1.d-5)
     *     WRITE(*,*)'ALGORITHM ERROR: nml:',(nml(i),i=1,3)
c     --- --- DEBUG_END --- ---

      do i=1,nvrt
c     ---Z axis rotation---
         ver1 = allvertex(1,iverplane(i,index))*coszfi -
     *        allvertex(2,iverplane(i,index))*sinzfi

         vertex(2,i) = allvertex(1,iverplane(i,index))*sinzfi +
     *        allvertex(2,iverplane(i,index))*coszfi
         vertex(3,i) = allvertex(3,iverplane(i,index))

c     ---Y axis rotation---
         vertex(1,i) =  ver1*cosyfi + vertex(3,i)*sinyfi
         vertex(3,i) = -ver1*sinyfi + vertex(3,i)*cosyfi
      enddo

 99   continue

c     --- --- DEBUG_BEGIN --- ---
c      do i=1,nvrt
c         PRINT *,'rotated-vertex(',nvrt,') = ',(vertex(j,i),j=1,3)
c      enddo
c     --- --- DEBUG_END --- ---

c     this is Jarvis's march algorithm for Covex Hull
c     -----------------------------------------------
c     first we find two points; Ymin & Ymax points

c     I think that colinear points are omitted, because of our 
c     previous construction
      iYmin = 1
      iYmax = 1
      do i=2,nvrt
         if(vertex(2,i) .lt. vertex(2,iYmin)) 
     *        iYmin = i
         if(vertex(2,i) .gt. vertex(2,iYmax))
     *        iYmax = i
      enddo

c     --- --- DEBUG_BEGIN --- ---
c      PRINT  *,'iYmin=',iYmin,';  iYmax=',iYmax
c     --- --- DEBUG_END --- ---

      nord=1
      iorder(1)=iYmin
      lfirst=.true.
      notend=.true.
      sign=1.0d0

      do while (notend .and. nord.lt.MAX_NORD)
         do 100 i=1,nvrt
c     --- --- DEBUG_BEGIN --- ---
c            PRINT *,'JARVIS MARCH: i=',i,'; iYmin=',iYmin
c            PRINT *,'              nord=',nord,
c     *           ';    iorder(nord)=',iorder(nord)
c    --- --- DEBUG_END --- --- 

            if(i.eq.iorder(nord)) goto 100
            if(lfirst) then
c               print *,'DEBUG> i-vertex           :',i,'::',
c     $              (vertex(ipol,i),ipol=1,3)
c               print *,'DEBUG> iorder(nord)-vertex:',iorder(nord),'::',
c     $              (vertex(ipol,iorder(nord)),ipol=1,3)
               minfi=PolarAngle(i,iorder(nord),vertex,sign)
               iTmp=i
               lfirst=.false.
c            endif
            else
               fi=PolarAngle(i,iorder(nord),vertex,sign)
c               print *,'DEBUG> i,iorder(nord)::',i,iorder(nord),
c     $              '; fi=',fi, minfi
               if(fi.lt.minfi) then
                  minfi=fi
                  iTmp=i
               endif
            endif
 100     continue
         
c         PRINT *,'after 100; iTmp=',iTmp,' iYmin=',iYmin
         if(iTmp.eq.iYmin)then
            notend=.false.
            goto 101
         endif

         lfirst=.true.
         nord=nord+1
c         print *, 'iorder(',nord,')=',iTmp
         iorder(nord)=iTmp
         if(iTmp.eq.iYmax) sign=-1.0d0

 101     continue
      enddo

      if (nord.ge.MAX_NORD)then
         write(*,*) 'ALGORITHM ERROR: nord became too large'
      endif            

 102  continue

c     ==================================================================
c     the sign of normal is determined by requiring that the origin 
c     is inside the convex polygon
      do i=1,3
         x(i) = allvertex(1,iverplane(iorder(i),index))
         y(i) = allvertex(2,iverplane(iorder(i),index))
         z(i) = allvertex(3,iverplane(iorder(i),index))
      enddo

      normal(1) = (y(1)-y(2))*z(3) + (y(3)-y(1))*z(2) + (y(2)-y(3))*z(1)
      normal(2) = (x(2)-x(1))*z(3) + (x(1)-x(3))*z(2) + (x(3)-x(2))*z(1)
      normal(3) = (x(1)-x(2))*y(3) + (x(3)-x(1))*y(2) + (x(2)-x(3))*y(1)
      normal(4) = 
     *     x(1)*(y(2)*z(3)-y(3)*z(2)) + 
     *     y(1)*(x(3)*z(2)-x(2)*z(3)) +
     *     z(1)*(x(2)*y(3)-x(3)*y(2))
      call normalize_vec(3,normal)
      write(12,*) nvrt
      if (normal(4).lt.0d0) then
         do i=nvrt,1,-1 
            write(12,'(2(3e15.7,1x))')
     $           (allvertex(j,iverplane(iorder(i),index)),j=1,3),
     *           (normal(j),j=1,3)            
         enddo
      else
         normal(1) = -normal(1)
         normal(2) = -normal(2)
         normal(3) = -normal(3)
         do i=1,nvrt
            write(12,'(2(3e15.7,1x))')
     $           (allvertex(j,iverplane(iorder(i),index)),j=1,3),
     *           (normal(j),j=1,3)
         enddo
      endif
c      WRITE(13,*) normal(4)
c      CALL VecProduct(vec1,vec2,normal) 
c     ===================================================================
      
      RETURN
      END


c     =========================================================
      REAL*8 FUNCTION PolarAngle(index, iref, vertex, sign)
c     =========================================================
c     fn. calculate polar angle between points index & iref
c     PI2 == 2*PI
      IMPLICIT NONE
      include 'voronoi.inc'
      REAL*8 vertex(3,*)
c      ,MAXERR      
c      PARAMETER( 
c     *     MAXERR    = 1.0d-10)
      REAL*8 v(3)
      REAL*8 vsize, sign, PI2
      PARAMETER (PI2=6.28318530717958647688d0)
      INTEGER index, iref, i

c      print *, 'Polar-angle: index,iref:', index, iref

      do i=1,3
         v(i)=vertex(i,index)-vertex(i,iref)
c         PRINT *,' VERTEX-index:: ', vertex(i,index),
c     $        '; VERTEX-iref:: ', vertex(i,iref)
      enddo
      
c     vsize shouldn't be ZERO
      vsize=dsqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3)) 
            
c     --- --- DEBUG_BEGIN --- ---
c      PRINT *,'vsize=',vsize,';  sign*v(1)/vsize=', sign*v(1)/vsize
      if(abs(vsize).lt.MAXERR)then
         PRINT *,'DEBUG> distance between two points lower than MAXERR'
         PRINT *,'          index = ', index, ';     iref = ', iref
         PRINT *,'       distance = ', abs(vsize)
         STOP
      endif
c     --- --- DEBUG_END --- ---
      
      PolarAngle = DACOS(sign*v(1)/vsize)
c      print *, 'Polar-angle.1=',PolarAngle
!      if(sign*v(2).lt.0.0d0 .and. abs(v(2)/vsize).gt.HALF001_DEG)
      if(sign*v(2).lt.0.0d0 .and. abs(v(2)/vsize).gt.MAXERR)
     *     PolarAngle=PI2 - PolarAngle 

c      print *, 'Polar-angle.2=',PolarAngle
      RETURN
      END

