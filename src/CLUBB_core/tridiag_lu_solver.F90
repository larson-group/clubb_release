module tridiag_lu_solvers

  ! Description:
  !   These routines solve lhs*soln=rhs using LU decomp. 
  !   
  !   LHS is stored in band diagonal form. 
  !     lhs = | lhs(0,1)  lhs(-1,1)      0          0          0          0      0
  !           | lhs(1,2)  lhs( 0,2)  lhs(-1,2)      0          0          0      0
  !           |     0     lhs( 1,3)  lhs( 0,3)  lhs(-1,3)      0          0      0
  !           |     0         0      lhs( 1,4)  lhs( 0,4)  lhs(-1,4)      0      0
  !           |     0         0          0      lhs( 1,5)  lhs( 0,5)  lhs(-1,5)  0 ...
  !           |    ...                                                                   
  !
  !    U is stored in band diagonal form 
  !     U = |   1    upper(1)  0       0       0       0       0
  !         |   0      1     upper(2)  0       0       0       0
  !         |   0      0       1     upper(3)  0       0       0
  !         |   0      0       0       1     upper(4)  0       0
  !         |   0      0       0       0       1     upper(5)  0  ...
  !         |  ...   
  !
  !   L is also stored in band diagonal form, but the lowest most band is equivalent to the 
  !   lowermost band of LHS, thus we don't need to store it
  !     L = | l_diag(1)       0           0            0          0        0   
  !         |  lower(2)      l_diag(2)    0            0          0        0    
  !         |    0           lower(3)   l_diag(3)      0          0        0          
  !         |    0            0         lower(4)    l_diag(4)     0        0     
  !         |    0            0           0         lower(5)    l_diag(5)   0   ...
  !         |  ...   
  ! 
  !
  !   To perform the LU decomposition, we go element by element. 
  !   First we start by noting that we want lhs=LU, so the first step of calculating
  !   L*U, by multiplying the first row of L by the columns of U, gives us
  !
  !     l_diag(1)*1        = lhs( 0,1)  =>  l_diag(1)   = lhs( 0,1)
  !     l_diag(1)*upper(1) = lhs(-1,1)  =>  upper(1)    = lhs(-1,1) / l_diag(1)
  !
  !   Now that we're passed the k=1 step, each following step uses all the bands,
  !   allowing us to write the general step
  !     
  !     lower(k)*1 = lhs(1,k)                         =>  lower(k)    = lhs( 1,k)
  !     lower(k)*upper(k-1)+l_diag(k)*1  = lhs( 0,k)  =>  l_diag(k) = lhs( 0,k) - lower(k)*upper(k-1)
  !     l_diag(k)*upper(k)               = lhs(-1,k)  =>  upper(k)  = lhs(-1,k) / l_diag(k)
  !
  !   This general step is done for k from 2 to ndim-1 (do k = 2, ndim-), and the last 
  !   step is tweaked similarly to the first one, where we disclude the upper band
  !   since it becomes no longer relevant. Note from this general step that the lower band
  !   is always equivalent to first subdiagonal band of lhs, thus we do not need to 
  !   calculate or store lower. Also note that we only ever need l_diag so that we can divide 
  !   by it, so instead we compute lower_diag_invrs to reduce divide operations.
  !
  !   After L and U are computed, normally we do forward substitution using L,
  !   then backward substitution using U to find the solution. This is replicated
  !   for every right hand side we want to solve for.
  !
  !
  ! References:
  !   none
  !------------------------------------------------------------------------


  use clubb_precision, only:  &
    core_rknd ! Variable(s)

  implicit none

  public :: tridiag_lu_solve

  private :: tridiag_lu_solve_single_rhs_multiple_lhs, tridiag_lu_solve_multiple_rhs_lhs

  interface tridiag_lu_solve
    module procedure tridiag_lu_solve_single_rhs_lhs
    module procedure tridiag_lu_solve_single_rhs_multiple_lhs
    module procedure tridiag_lu_solve_multiple_rhs_lhs
  end interface

  private ! Default scope

  contains

  !=============================================================================
  subroutine tridiag_lu_solve_single_rhs_lhs( ndim, lhs, rhs, &
                                              soln )
    ! Description:
    !   Written for single RHS and single LHS.
    !------------------------------------------------------------------------

    implicit none
 
    ! ----------------------- Input Variables -----------------------
    integer, intent(in) :: &
      ndim  ! Matrix size 

    real( kind = core_rknd ), intent(in), dimension(ndim) ::  &
      rhs       ! 

    ! ----------------------- Input/Output Variables -----------------------
    real( kind = core_rknd ), intent(inout), dimension(-1:1,ndim) :: &
      lhs   ! Matrices to solve, stored using band diagonal vectors
            ! -2 is the uppermost band, 2 is the lower most band, 0 is diagonal

    ! ----------------------- Output Variables -----------------------
    real( kind = core_rknd ), intent(out), dimension(ndim) ::  &
      soln     ! Solution vector

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(ndim) ::  &
      upper,         & ! First U band 
      lower_diag_invrs ! Inverse of the diagonal of L

    integer :: k    ! Loop variables

    ! ----------------------- Begin Code -----------------------
       
    !$acc data create( upper, lower_diag_invrs ) &
    !$acc      copyin( rhs, lhs ) &
    !$acc      copyout( soln )
    
    !$acc kernels
    
    lower_diag_invrs(1) = 1.0_core_rknd / lhs(0,1)
    upper(1)            = lower_diag_invrs(1) * lhs(-1,1) 

    do k = 2, ndim-1
      lower_diag_invrs(k) = 1.0_core_rknd / ( lhs(0,k) - lhs(1,k) * upper(k-1)  )
      upper(k)            = lower_diag_invrs(k) * lhs(-1,k) 
    end do

    lower_diag_invrs(ndim) = 1.0_core_rknd / ( lhs(0,ndim) - lhs(1,ndim) * upper(ndim-1)  )

    soln(1)   = lower_diag_invrs(1) * rhs(1) 

    do k = 2, ndim
      soln(k) = lower_diag_invrs(k) * ( rhs(k) - lhs(1,k) * soln(k-1) )
    end do

    do k = ndim-1, 1, -1
      soln(k) = soln(k) - upper(k) * soln(k+1)
    end do

    !$acc end kernels

    !$acc end data

  end subroutine tridiag_lu_solve_single_rhs_lhs

  !=============================================================================
  subroutine tridiag_lu_solve_single_rhs_multiple_lhs( ndim, ngrdcol, lhs, rhs, &
                                                       soln )
    ! Description:
    !   Written for single RHS and multiple LHS.
    !------------------------------------------------------------------------

    implicit none
 
    ! ----------------------- Input Variables -----------------------
    integer, intent(in) :: &
      ndim,   & ! Matrix size 
      ngrdcol   ! Number of grid columns

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,ndim) ::  &
      rhs       ! 

    ! ----------------------- Input/Output Variables -----------------------
    real( kind = core_rknd ), intent(inout), dimension(-1:1,ngrdcol,ndim) :: &
      lhs   ! Matrices to solve, stored using band diagonal vectors
            ! -2 is the uppermost band, 2 is the lower most band, 0 is diagonal

    ! ----------------------- Output Variables -----------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,ndim) ::  &
      soln     ! Solution vector

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  &
      upper,         & ! First U band 
      lower_diag_invrs ! Inverse of the diagonal of L

    integer :: i, k    ! Loop variables

    ! ----------------------- Begin Code -----------------------
       
    !$acc data create( upper, lower_diag_invrs ) &
    !$acc      copyin( rhs, lhs ) &
    !$acc      copyout( soln )
    
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      lower_diag_invrs(i,1) = 1.0_core_rknd / lhs(0,i,1)
      upper(i,1)            = lower_diag_invrs(i,1) * lhs(-1,i,1) 
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do k = 2, ndim-1
      do i = 1, ngrdcol
        lower_diag_invrs(i,k) = 1.0_core_rknd / ( lhs(0,i,k) - lhs(1,i,k) * upper(i,k-1)  )
        upper(i,k)            = lower_diag_invrs(i,k) * lhs(-1,i,k) 
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      lower_diag_invrs(i,ndim) = 1.0_core_rknd / ( lhs(0,i,ndim) - lhs(1,i,ndim) * upper(i,ndim-1)  )
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol 

      soln(i,1)   = lower_diag_invrs(i,1) * rhs(i,1) 

      do k = 2, ndim
        soln(i,k) = lower_diag_invrs(i,k) * ( rhs(i,k) - lhs(1,i,k) * soln(i,k-1) )
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol 
      do k = ndim-1, 1, -1
        soln(i,k) = soln(i,k) - upper(i,k) * soln(i,k+1)
      end do

    end do
    !$acc end parallel loop

    !$acc end data

  end subroutine tridiag_lu_solve_single_rhs_multiple_lhs

  
  !=============================================================================
  subroutine tridiag_lu_solve_multiple_rhs_lhs( ndim, nrhs, ngrdcol, lhs, rhs, &
                                                soln )
    ! Description:
    !   Written for multiple RHS and multiple LHS.
    !------------------------------------------------------------------------

    implicit none
 
    ! ----------------------- Input Variables -----------------------
    integer, intent(in) :: &
      ndim,   & ! Matrix size 
      nrhs,   & ! Number of right hand sides
      ngrdcol   ! Number of grid columns

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,ndim,nrhs) ::  &
      rhs       ! 

    ! ----------------------- Input/Output Variables -----------------------
    real( kind = core_rknd ), intent(inout), dimension(-1:1,ngrdcol,ndim) :: &
      lhs   ! Matrices to solve, stored using band diagonal vectors
            ! -2 is the uppermost band, 2 is the lower most band, 0 is diagonal

    ! ----------------------- Output Variables -----------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,ndim,nrhs) ::  &
      soln     ! Solution vector

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  &
      upper,            & ! First U band 
      lower_diag_invrs    ! Inverse of the diagonal of L

    integer :: i, k, j    ! Loop variables

    ! ----------------------- Begin Code -----------------------
       
    !$acc data create( upper, lower_diag_invrs ) &
    !$acc      copyin( rhs, lhs ) &
    !$acc      copyout( soln )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      lower_diag_invrs(i,1) = 1.0_core_rknd / lhs(0,i,1)
      upper(i,1)            = lower_diag_invrs(i,1) * lhs(-1,i,1) 
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      do k = 2, ndim-1
        lower_diag_invrs(i,k) = 1.0_core_rknd / ( lhs(0,i,k) - lhs(1,i,k) * upper(i,k-1)  )
        upper(i,k)            = lower_diag_invrs(i,k) * lhs(-1,i,k) 
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      lower_diag_invrs(i,ndim) = 1.0_core_rknd / ( lhs(0,i,ndim) - lhs(1,i,ndim) * upper(i,ndim-1)  )
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do j = 1, nrhs
      do i = 1, ngrdcol 

        soln(i,1,j)   = lower_diag_invrs(i,1) * rhs(i,1,j) 

        do k = 2, ndim
          soln(i,k,j) = lower_diag_invrs(i,k) * ( rhs(i,k,j) - lhs(1,i,k) * soln(i,k-1,j) )
        end do
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do j = 1, nrhs
      do i = 1, ngrdcol 
        do k = ndim-1, 1, -1
          soln(i,k,j) = soln(i,k,j) - upper(i,k) * soln(i,k+1,j)
        end do
      end do
    end do
    !$acc end parallel loop

    !$acc end data

  end subroutine tridiag_lu_solve_multiple_rhs_lhs

end module tridiag_lu_solvers
