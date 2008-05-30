      module matrix_operations

      implicit none

      public :: gaussj, matmult

      private ! Defualt scope

      contains

      subroutine gaussj( a_in, n, np, a )

      implicit none

      integer m, mp, nmax

      integer, intent(in) :: n, np
! Output: a = Matrix inverse

      double precision, intent(out) :: a(np,np)
! Input: a_in = Original matrix

      double precision, intent(in)  :: a_in(np,np)

      parameter (NMAX=50)
      integer :: i, icol, irow, j, k, l, ll
      integer :: indxc(NMAX), indxr(NMAX), ipiv(NMAX)
      double precision big, dum, pivinv
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
!----------------------------------------------------------------------
        subroutine matmult( A, rowsa, colsa, rowsap, colsap,
     .                      B, rowsb, colsb, rowsbp, colsbp, C )

        implicit none

! Input variables
        integer, intent(in) :: rowsa, colsa, rowsb, colsb
        integer, intent(in) :: rowsap, colsap, rowsbp, colsbp
        double precision, intent(in) :: 
     .                  A(rowsap,colsap), B(rowsbp,colsbp)

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
      end module matrix_operations
