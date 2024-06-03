c     ------------------------------------------------------------------------
      subroutine DIMIFY_VEC(DVEC,IDIM)
c     ------------------------------------------------------------------------
      real*8 DVEC(3,3)

      if (idim.eq.1)then
         dvec(1,2)=0.0d0
         dvec(1,3)=0.0d0
         
         dvec(2,1)=0.0d0
         dvec(2,2)=1.0d0
         dvec(2,3)=0.0d0

         dvec(3,1)=0.0d0
         dvec(3,2)=0.0d0
         dvec(3,3)=1.0d0
      elseif (idim.eq.2)then
         dvec(1,3)=0.0d0
         dvec(2,3)=0.0d0

         dvec(3,1)=0.0d0
         dvec(3,2)=0.0d0

         dvec(3,3)=1.0d0
      endif
      return
      END

c     ----------------------------------------------------------------
      INTEGER function string_length(word)
c     trims the string from both sides and returns its trimmed length
c     ----------------------------------------------------------------
      CHARACTER word*(*)
      word          = adjustl(word)
      string_length = len_trim(word)
      RETURN
      END

c     ===============================
      SUBROUTINE WVEC(WORD,VEC)
c     subroutine write vector to unit12
      CHARACTER*(*) WORD
      REAL*8 VEC(3,3)

      WRITE(12,'(1x,a)') word
      DO J=1,3
         WRITE(12,'(3(1x,f15.10)/3(1x,f15.10)/3(1x,f15.10))')
     $        (VEC(J,I),I=1,3)
      ENDDO
      RETURN
      END

c     ===============================
      SUBROUTINE WIVEC(WORD,VEC)
c     WARNING:: reciprocal vectors are written in transpose matrix; 
c     take care of that
c     subroutine write vector to unit12
      CHARACTER*(*) WORD
      REAL*8 VEC(3,3)

      WRITE(12,'(1x,a)') word
      DO J=1,3
         WRITE(12,'(3(1x,f15.10)/3(1x,f15.10)/3(1x,f15.10))')
     $        (VEC(I,J),I=1,3)
      ENDDO
      RETURN
      END

c     =========================
      INTEGER FUNCTION C2I(ARG)
      CHARACTER*(*) ARG

c      print *,'arg>',arg
      I=INDEX(ARG,' ')-1
      IF(I.EQ.0) I=LEN(ARG)               
      C2I=0
      DO M=1,I
         C2I=C2I + (ICHAR(ARG(M:M)) - 48)*10**(I-M)
      ENDDO
c      print *,'iarg>',c2i
      RETURN 
      END

c     ============================================================
      REAL*8 FUNCTION DET3X3D(x11,x12,x13,x21,x22,x23,x31,x32,x33)
      REAL*8 x11,x12,x13,x21,x22,x23,x31,x32,x33
      
      DET3X3D=x11*x22*x33 + x12*x23*x31 + x13*x21*x32-
     *        x31*x22*x13 - x32*x23*x11 - x33*x21*x12
      RETURN
      END



c     ============================================================
      REAL*8 FUNCTION DET2X2D(x11,x12,x21,x22)
      REAL*8 x11,x12,x21,x22
      
      DET2X2D=x11*x22 - x21*x12
      RETURN
      END



c     =======================================
      SUBROUTINE VecProduct(vec1,vec2,resvec)
      REAL*8 vec1(3),vec2(3),resvec(3)
      
      resvec(1) = vec1(2)*vec2(3) - vec2(2)*vec1(3)
      resvec(2) = vec1(3)*vec2(1) - vec2(3)*vec1(1)
      resvec(3) = vec1(1)*vec2(2) - vec2(1)*vec1(2)
      RETURN
      END

c     ====================================
      REAL*8 FUNCTION ScalarProduct(v1,v2)
      REAL*8 v1(3), v2(3)

      ScalarProduct = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
      RETURN
      END

c     ============================
      REAL*8 FUNCTION VecSize(vec)
      REAL*8 vec(3)
      
      VecSize = dsqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
      RETURN
      END




c     ==================================
      SUBROUTINE MATMULT(A,NI,NJ,B,NK,C)
      REAL*8 A(NJ,NI),B(NK,NJ),C(NK,NI)
C     *** THIS SUBROUTINE MULTIPLIES TWO MATRICES ***
C     
C     C(KCOL,IROW)=SUM(J){ A(JCOL,IROW)*B(KCOL,JROW) }
C     
C     MATRIX C MUST BE ALREADY INITIALISED
C     ***********************************************

      DO K=1,NK
         DO I=1,NI
            DO J=1,NJ
               C(K,I)=C(K,I)+A(J,I)*B(K,J)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END



c     ================================
      REAL*8 FUNCTION VECSIZE3(DVEC,N)
      REAL*8 DVEC(3,3),DSQRT

      VECSIZE3=DSQRT(DVEC(N,1)**2+DVEC(N,2)**2+DVEC(N,3)**2)
      RETURN
      END



c     =======================================
      REAL*8 FUNCTION DIST(X1,Y1,Z1,X2,Y2,Z2)
      REAL*8 X1,Y1,Z1,X2,Y2,Z2,DSQRT

      DIST=DSQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
      RETURN
      END


      
c     ===========================
      LOGICAL FUNCTION EQUAL(A,B)
      REAL*8 A,B,MINTOL
      PARAMETER (MINTOL=1.d-6)
      
C     print *,'A=',A,'B=',B
      IF(ABS(A-B).LT.MINTOL) THEN
         EQUAL=.TRUE.
      ELSE
         EQUAL=.FALSE.
      ENDIF
      RETURN
      END



c     =========================
      SUBROUTINE MatTranspose(src,dst,m,n)
      REAL*8 src(m,n), dst(n,m)
      
      do i=1,n
         do j=1,m
            dst(i,j)=src(j,i)
         enddo
      enddo
      RETURN
      END


c     =========================
      SUBROUTINE MATCOPY(a,b,m,n)
      REAL*8 A(m,n), B(m,n)
      
      do i=1,n
         do j=1,m
            b(j,i)=a(j,i)
         enddo
      enddo
      RETURN
      END


c     ===========================
      SUBROUTINE ZeroMat(a,m,n)
      REAL*8 A(m,n)
      
      do i=1,n
         do j=1,m
            a(j,i) = 0.0d0
         enddo
      enddo
      RETURN
      END

c     ===========================
      SUBROUTINE IdentityMat(a,m,n)
      REAL*8 A(m,n)
      
      do i=1,n
         do j=1,m
            a(j,i) = 0.0d0
         enddo
         a(i,i) = 1.0d0
      enddo
      RETURN
      END

c     ===========================
      SUBROUTINE  Mat33ToMat22(m)
      real*8 m(3,3)
      do i=1,3
         m(3,i) = 0.0d0
         m(i,3) = 0.0d0
      enddo
      m(3,3) = 1.0d0
      return
      end


c     ===========================
      SUBROUTINE  Mat33ToMat11(m)
      real*8 m(3,3)
      call Mat33ToMat22(m)
      do i=1,3
         m(2,i) = 0.0d0
         m(i,2) = 0.0d0
      enddo
      m(2,2) = 1.0d0
      return
      end

c     ===========================
      SUBROUTINE INVERT2x2(vec)
      REAL*8 vec(3,3), tmp(3,3)
      REAL*8 v1(3), v2(3), v3(3)

      call Mat33ToMat22(vec)
      Call IdentityMat(tmp, 3, 3)
      do i=1,3
         v1(i) = vec(1,i)
         v2(i) = vec(2,i)
         v3(i) = 0.0d0
      enddo
      v3(3) = 1.0d0

      Call RecipVec(tmp, 1, v1, v2, v3)
      Call RecipVec(tmp, 2, v2, v3, v1)
      Call MatTranspose(tmp, vec, 3, 3)
      
      RETURN
      END

c     ===========================
      SUBROUTINE RecipVec(vec, i, v1, v2, v3)
      REAL*8 vec(3,3), v1(3), v2(3), v3(3)
      REAL*8 vp(3), sp, ScalarProduct

      Call VecProduct( v2, v3, vp )

      sp = ScalarProduct(v1, vp)
      do jj=1,3
         vec(i,jj) = vp(jj) / sp ! devision by zero can occur here
      enddo

      RETURN
      END


      SUBROUTINE InvertSqMat123(mat, invmat, ldim) 
      real*8 mat(3,3), invmat(3,3), tmpmat(3,3)
      
      call MatCopy(mat, invmat, 3, 3)      
      if (ldim .le. 1) then
         call Mat33ToMat11(invmat)
         if (ldim .eq. 0) then
            invmat(1,1) = 1.0d0
         else
            invmat(1,1) = 1.0d0/mat(1,1)
         endif
      else if(ldim .eq. 2) then
         call Mat33ToMat22(invmat)
         call Invert2x2(invmat)
      else if(ldim .eq. 3) then
         call Invert3x3(mat,invmat)
      else
         stop 'dimension of matrix greater then 3'
      endif
      return
      END

c     -------------------------------------------------
      integer function iCountFields(word)
c     counts the number fo fields in character
c     -------------------------------------------------
      character word*(*), auxword*255
      logical lspace, WhiteSpace

      iCountFields=0
      lspace=.false.
      ilen=len(word)

c     left trim
      do i=1,ilen
         if ( .not. WhiteSpace(word(i:i)) ) goto 1
      enddo
 1    continue
      ifirst=i

c     right trim
      do i=ilen,1,-1
         if ( .not. WhiteSpace(word(i:i)) ) goto 2
      enddo
 2    continue
      ilast=i

      if (ifirst .lt. ilast) iCountFields=1
c      print *,'TT:', ilen,ifirst,ilast

      do i=ifirst,ilast
         if ( WhiteSpace(word(i:i)) ) then
            if (.not.lspace) then
               lspace=.true.
               iCountFields=iCountFields+1
            endif
         else
            lspace=.false.
         endif
      enddo
      return
      END

c     -------------------------------------------------------------------------
      logical function WhiteSpace(char)
c     returns true if character is white-space or TAB
c     -------------------------------------------------------------------------
      character char
      WhiteSpace = .false.
      if ((char.eq.' ') .or. (char.eq.'\t')) WhiteSpace = .true.
      return
      END


c     ------------------------------------------------------------------------
      subroutine Invert3x3(mat, invmat)
c     Calculate the 3x3 inversion matrix. This routine is based on
c     routine invmat3 of FPMP/CP90 packages. Copyright (C) 2002 FPMD
c     group and CP90 group. The routine is distributed under the terms of
c     the GNU General Public License.
c     ------------------------------------------------------------------------
      implicit none
      real*8 mat(3,3), invmat(3,3), norm
      real*8 d11, d12, d13, d22, d23, d33, d21, d31, d32

      d11=mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)
      d12=mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3)
      d13=mat(2,1)*mat(3,2)-mat(3,1)*mat(2,2)
      d22=mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1)
      d23=mat(3,1)*mat(1,2)-mat(1,1)*mat(3,2)
      d33=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
      d21=mat(3,2)*mat(1,3)-mat(1,2)*mat(3,3)
      d31=mat(1,2)*mat(2,3)-mat(2,2)*mat(1,3)
      d32=mat(1,3)*mat(2,1)-mat(1,1)*mat(2,3)

      norm=mat(1,1)*d11+mat(1,2)*d12+mat(1,3)*d13
c     check for singular matrix
      if (ABS(norm).lt.1.d-20) stop 'singular matrix in Invert3x3'

      invmat(1,1)=d11/norm
      invmat(2,2)=d22/norm
      invmat(3,3)=d33/norm
      invmat(1,2)=d21/norm
      invmat(1,3)=d31/norm
      invmat(2,3)=d32/norm
      invmat(2,1)=d12/norm
      invmat(3,1)=d13/norm
      invmat(3,2)=d23/norm

      return
      end


      subroutine normalize_vec(n,vec)
      implicit none
      integer n, i
      real*8  vec(*), norm
      
      norm=0.0d0
      do i=1,n
         norm = norm + vec(i)*vec(i)
      enddo
      norm = dsqrt(norm)
c     if norm is zero don't normalize
      if (norm.lt.1.d-20) return
      do i=1,n
         vec(i) = vec(i) / norm
      enddo
      return
      end
