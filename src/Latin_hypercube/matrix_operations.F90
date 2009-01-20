! $Id$
      module matrix_operations

      implicit none

      public :: gaussj, matmult

      private ! Defualt scope

      contains

      subroutine gaussj( a_in, n, np, a )

      implicit none

!     integer m, mp

      integer nmax

      integer, intent(in) :: n, np
! Output: a = Matrix inverse

      double precision, intent(out) :: a(np,np)
! Input: a_in = Original matrix

      double precision, intent(in)  :: a_in(np,np)

      parameter (NMAX=50)
      integer :: i, icol, irow, j, k, l, ll
      integer :: indxc(NMAX), indxr(NMAX), ipiv(NMAX)
      double precision big, dum, pivinv

      ! Default Initialization
      irow = 0 
      icol = 0

      a = a_in
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
!          do 15 l=1,m
!            dum=b(irow,l)
!            b(irow,l)=b(icol,l)
!            b(icol,l)=dum
!15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) stop 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
!        do 17 l=1,m
!          b(icol,l)=b(icol,l)*pivinv
!17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
!            do 19 l=1,m
!              b(ll,l)=b(ll,l)-b(icol,l)*dum
!19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue

      return
      end subroutine gaussj
!----------------------------------------------------------------------
! Subroutine to calculate the product of two matrices: 
!           C(rowsa,colsb) = A(rowsa,colsa) * B(rowsb,colsb). 
! We must have colsa=rowsb; otherwise the product is undefined.
! All multiplication is done in double precision.

! Input:  A = left matrix
!         B = right matrix
!         rowsa, colsa = actual number of rows and columns for A
!         rowsap,colsap = physical (maximum declared) dimensions of A
!         rowsb, colsb = number of rows and columns for B
!         rowsbp,colsbp = physical (maximum declared) dimensions of B

! Output: C = product matrix with size rowsa x colsb

!     Note: Could we just use f90 intrinsic function matmul?  Is this a
!     legacy Fortran 77 subroutine?
!     -dschanen 14 June 2008
!----------------------------------------------------------------------
        subroutine matmult( A, rowsa, colsa, rowsap, colsap, & 
                            B, rowsb, colsb, rowsbp, colsbp, C )

        implicit none

! Input variables
        integer, intent(in) :: rowsa, colsa, rowsb, colsb
        integer, intent(in) :: rowsap, colsap, rowsbp, colsbp
        double precision, intent(in) ::  & 
                        A(rowsap,colsap), B(rowsbp,colsbp)

! Output variables
        double precision, intent(out) :: C(rowsa,colsb) 

! Local variables
        integer i, j, k
!        double precision sum    ! reserved word, might cause problems later on
        double precision matsum

! Check whether matrices A and B have consistent dimensions
        if ( colsa /= rowsb ) then
          print *, 'Error: Matrix dims inconsistent in matmult.'
          stop
        endif

! Perform matrix multiplication
       do i = 1, rowsa
         do j = 1, colsb
           matsum = 0.d0
             do k = 1, colsa
               matsum = matsum + A(i,k) * B(k,j)
             enddo
             C(i,j) = matsum
         enddo
       enddo

       return
       end subroutine matmult

!-----------------------------------------------------------------------

        subroutine band_mult( trans, ndim, mdim, nsup, nsub, yinc, xinc, & 
                              alpha, beta, lhs, xvec, yvec )
!       Description:
!       Wrapper subroutine for banded matrix by vector multiplication in
!       the level 2 BLAS library.

!       References:
!       <http://www.netlib.org/blas/>
!-----------------------------------------------------------------------

        implicit none

        ! External
        ! Level 2 BLAS to multiply a vector by a band diagonal matrix
        external :: sgbmv, dgbmv

        character(len=1), intent(in) ::  & 
        trans  ! Whether to use the transposition of the lhs matrix

        integer, intent(in) :: & 
        ndim, mdim,   & ! Dimensions of the matrix when not compact
        nsup, nsub,   & ! Super and Sub diagonals
        yinc, xinc   ! Increments of y and x vector

        real, intent(in) :: & 
        alpha,        & ! Coefficient of the matrix lhs
        beta         ! Coefficient of Y vector

        real, dimension(nsup+nsub+1, ndim), intent(in) :: & 
        lhs          ! The matrix 'A' in the blas subroutine

        real, dimension(ndim), intent(in) :: & 
        xvec ! The vector X

        real, dimension(ndim), intent(inout) :: & 
        yvec ! The vector Y


!-----------------------------------------------------------------------
!       *** BLAS 2 routine ***
!       SUBROUTINE DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!-----------------------------------------------------------------------
        ! Multiply so that Y := alpha*A*X + beta*Y

        if ( kind( lhs(1,1) ) == 4 ) then
          call sgbmv( trans, ndim, mdim, nsub, nsup,  & 
                      alpha, lhs, nsup+nsub+1,  & 
                      xvec, xinc, beta, yvec, yinc )

        else if ( kind( lhs(1,1) ) == 8 ) then
          call dgbmv( trans, ndim, mdim, nsub, nsup,  & 
                      alpha, lhs, nsup+nsub+1,  & 
                      xvec, xinc, beta, yvec, yinc )
        else
          stop "Cannot multiply this precision"
        end if

        return
        end subroutine band_mult
!-----------------------------------------------------------------------

      end module matrix_operations
