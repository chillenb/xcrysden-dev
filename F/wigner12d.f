      SUBROUTINE WignerSeitz1D(vec,word)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 vec(3,3)
      CHARACTER*(*) word
      REAL*8 vertex(3,3), nml(3)

      Call ZeroMat(vertex, 3, 3)
      vertex(1,1) = -0.5d0 * vec(1,1)
      vertex(2,1) =  0.5d0 * vec(1,1)
      nml(1)      =  0.0d0
      nml(2)      =  1.0d0
      nml(3)      =  0.0d0

      Write(12,'(1x,a)') word
      write(12,*) 2
      do i=1,2
         write(12,*) (vertex(i,j),j=1,3), (nml(j),j=1,3)
      enddo
      Write(12,'(1x,a,"_END")') word

      RETURN
      END


      SUBROUTINE WignerSeitz2D(srcvec,word)
      IMPLICIT REAL*8 (a-h,o-z)

      include 'voronoi.inc'
      REAL*8 srcvec(3,3), vec(3,3)
      CHARACTER*(*) word

c     this SUBROUTINE is not efficient but is very simple !!!
      REAL*8 
     *     vertex(3,MAX_VERT),!VERTICES OF VORONOI POLIHEDRA
     *     point(3,26),       !equation of plane: ax+by+cz=d; point = (a,b,c,d)
     *     vect1(3),          !arbitrary vector 1
     *     vect2(3),          !arbitrary vector 2
     *     nml(3),            ! nml=vect1 x vect2
     *     t(3)
      INTEGER 
     *     iverplane2(MAX_PL_VER,1), !copy of iverplane
     *     icc(3)
      LOGICAL 
     *     activ(MAX_PL_VER),
     *     newvertex, lnew, lline, equal12, equal13, equal23, equalvert

c     NPOI.........number of point
c     NVER.........number of vertices

C     ===========================================
C           |11   21|   31|           |x1  x2|
C     VEC = |12   22|   32|  => VEC = |y1  y2|
C     ===========================================

      
!     
!     if lattice vectors are not reduced, reduce them
!
      call ReduceBasis2d(srcvec,vec)

!     ***

      Write(12,'(1x,a)') word

c     make all neighbour points of origin (0,0); there are 8 such points
      NPOI=0
      do i = -1, 1
         do j = -1, 1
            if(i.ne.0 .or. j.ne.0)then !disregard point(0,0)
               di=dble(i)
               dj=dble(j)
               NPOI=NPOI+1
c     equation of line: ax+by=c, point = (a,b,c)
               point(1,NPOI)=di*vec(1,1) + dj*vec(2,1)
               point(2,NPOI)=di*vec(1,2) + dj*vec(2,2)
               point(3,NPOI)=0.5d0 * (
     *              point(1,NPOI)*point(1,NPOI) +
     *              point(2,NPOI)*point(2,NPOI))
            endif
         enddo
      enddo

c     --- --- DEBUG_BEGIN --- ---
c      PRINT *,'NUMBER OF POINTS: ',NPOI
c      do i=1,NPOI
c         write(*,*) (point(j,i),j=1,4)
c      enddo
c      PRINT *,''
c     --- --- DEBUG_END --- ---

C     find vertices as intersections of two lines;
C     lines are made at the middle of origin and points
      NVER=0
      do i=1,NPOI-1
         do j=i+1,NPOI
            det=det2x2d(
     $           point(1,i), point(1,j),
     $           point(2,i), point(2,j) )
c     print *,'DEBUG> det = ',det
            if(abs(det).gt.MAXERR)then
               detx=det2x2d(
     *              point(3,i), point(3,j), 
     *              point(2,i), point(2,j) )
 
               dety=det2x2d(
     *              point(1,i), point(1,j),
     *              point(3,i), point(3,j))

               px=detx/det
               py=dety/det               
c     PRINT *,'px,py,pz>',px,py,pz 
c     print *,'DEBUG> px,py,pz:: ',px,py,pz

C     if vertex is closer to origin than to any other point, 
C     than we have new vertex
               l=0
               newvertex=.true.
               do while (newvertex .and. l.lt.NPOI)
                  l=l+1
                  dist1=dsqrt(px*px + py*py)
                  dist2=dsqrt(
     *                 (px-point(1,l))*(px-point(1,l)) +
     *                 (py-point(2,l))*(py-point(2,l)) )
c     print *,
c     *                    'DEBUG> l:',l,' ==>dist1 = ',
c     *                    dist1,'; dist2 = ',dist2 
                  if(dist1.gt.(dist2 + MAXERR)) newvertex=.false.
               enddo
c     print *,'DEBUG> l=',l,'newvertex=',newvertex
               if(newvertex) then
                  NVER=NVER+1
                  vertex(1,NVER) = px
                  vertex(2,NVER) = py
                  vertex(3,NVER) = 0.0d0
                  activ(NVER)    = .true.
               endif
            endif
         enddo
      enddo

c     --- --- DEBUG_BEGIN --- ---
*      PRINT *,'ALL VERICES OF VORONOI: NVER=',NVER
*      do i=1,NVER
*         WRITE(*,*) (vertex(j,i),j=1,3),'; --> ',
*     *        (vindex(j,i),j=1,3)
*      enddo
*      WRITE(*,*) ''
c     --- --- DEBUG_END --- ---

c     some vertices that belong to one line may be identical or colinear, 
c     get rid of that vertices
      do ic=1,NVER-2
         if ( activ(ic) ) then
            do jc=ic+1,NVER-1
               if ( activ(jc) ) then
                  do kc=jc+1,NVER
                     if ( activ(kc) ) then
c     get rid of identical points
                        equal12 =
     *           abs(vertex(1,ic)-vertex(1,jc)).lt.MAXERR .and.
     *           abs(vertex(2,ic)-vertex(2,jc)).lt.MAXERR 
                        equal13 = 
     *           abs(vertex(1,ic)-vertex(1,kc)).lt.MAXERR .and.
     *           abs(vertex(2,ic)-vertex(2,kc)).lt.MAXERR

                        equal23 =
     *           abs(vertex(1,jc)-vertex(1,kc)).lt.MAXERR .and.
     *           abs(vertex(2,jc)-vertex(2,kc)).lt.MAXERR
                        
                        equalvert=.true.
                        if(equal12.and.equal13)then 
c     all three points are identical
c     disregard with lower indices: ic < jc < kc
                           activ(ic)=.false.
                           activ(jc)=.false.
                        elseif(equal12)then
                           activ(ic)=.false.
                        elseif(equal13)then
                           activ(ic)=.false.
                        elseif(equal23)then
                           activ(jc)=.false.
                        else
                           equalvert=.false.
                        endif
                              
c     --- --- DEBUG_BEGIN --- ---
c     PRINT *,'i=',i,';  equalvert=',equalvert
c     --- --- DEBUG_BEGIN --- ---
                        
                        if(.not.equalvert)then
                           do lc=1,3
                              vect1(lc) = vertex(lc,jc)
     *                             - vertex(lc,ic) 
                              vect2(lc) = vertex(lc,kc)
     *                             - vertex(lc,ic) 
                           enddo

                           CALL VecProduct(vect1,vect2,nml)
                           size = VecSize(nml)
c     --- --- DEBUG_BEGIN --- ---
c                              PRINT *,'SIZE>',
c     *                                (vertex(lc,iverplane(ic,i))
c     *                             ,vertex(lc,iverplane(jc,i))
c     *                             ,vertex(lc,iverplane(kc,i)),lc=1,3)
c                              PRINT *,'SIZE::',size,'MAXERR::',MAXERR
c     --- --- DEBUG_END --- ---
                           if(abs(size).lt.MAXERR)then
c     colinear points; disregard middle point
                              icc(1)=ic
                              icc(2)=jc
                              icc(3)=kc
                              t(1)=0.0d0
                              t(2)=1.0d0
                              t(3)=gett3(
     *                             vertex(1,ic),
     *                             vertex(1,jc),
     *                             vertex(1,kc))
                              do lc=2,3
                                 icc1=icc(lc)
                                 do lc1=lc-1,1,-1
                                    if(t(lc1).lt.t(icc1)) 
     *                                   goto 101
                                    icc(lc1+1)=icc(lc1)
                                 enddo
                                 lc1=0
 101                             icc(lc1+1)=icc1
                              enddo
                              
c     --- --- DEBUG_BEGIN --- ---
c                                    PRINT *,'i=',i,
c     *                                   ';SORTING ORDER OF VERT:'
c                                    PRINT *,'ic=',ic,';  ',
c     *    (vertex(iver,iverplane(icc(1),i)),iver=1,3)
c                                    PRINT *,'jc=',jc,';  ',
c     *    (vertex(iver,iverplane(icc(2),i)),iver=1,3)
c                                    PRINT *,'kc=',kc,';  ',
c     *    (vertex(iver,iverplane(icc(3),i)),iver=1,3)

c     --- --- DEBUG_END --- ---
                              activ(icc(2))=.false. !middle point disregarded
                           endif !if(abs(size).lt.MAXERR)then
                        endif   !if(.not.equalvert)then
                     endif
                  enddo         !kc
               endif
            enddo               !jc
         endif
      enddo                     !ic
      
      NVER2=0
      do ic=1,NVER
c     --- --- DEBUG_BEGIN --- ---
c                  PRINT *,'ic=',ic,';  activ(ic,i)=',activ(ic,i)
c                  PRINT *,(vertex(lc,iverplane(ic,i)),lc=1,3)
c     --- --- DEGUB_END --- ---
         if(activ(ic))then
            NVER2=NVER2+1
            iverplane2(NVER2,1)=ic
         endif
      enddo
         
      if(NVER2.gt.2) then         
c     --- --- DEBUG_BEGIN --- ---
c            WRITE(*,*)'PLANE N.:',i,'; NVPL2=',NVPL2,'; NVPL',NVPL,
c     *           ';  INDEX LIST/VERTEX LIST'
c         do ideb=1,NVPL2
c     WRITE(*,*) iverplane2(ideb,nplane),';   ',
c     *              (vertex(j,iverplane2(ideb,nplane)),j=1,3)
c         enddo
c            WRITE(*,*)'NUMBER OF PLANES: ',nplane
c     --- --- DEBUG_END --- ---
         
c     now we have nverplane(i) vertices in plane i; make Convex Hull
c     for that plane
         CALL ConvexHull(1,iverplane2,NVER2,vertex)
      endif

      Write(12,'(1x,a,"_END")') word

      RETURN
      END
