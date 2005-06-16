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
