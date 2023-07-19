module penta_lu_solvers

  ! Description:
  !   These routines solve lhs*soln=rhs using LU decomp. 
  !   
  !   LHS is stored in band diagonal form. 
  !     lhs = | lhs(0,1)  lhs(-1,1)  lhs(-2,1)     0           0          0          0
  !           | lhs(1,2)  lhs( 0,2)  lhs(-1,2)  lhs(-2,2)      0          0          0
  !           | lhs(2,3)  lhs( 1,3)  lhs( 0,3)  lhs(-1,3)  lhs(-2,3)      0          0
  !           |     0     lhs( 2,4)  lhs( 1,4)  lhs( 0,4)  lhs(-1,4)  lhs(-2,4)      0
  !           |     0         0      lhs( 2,5)  lhs( 1,5)  lhs( 0,5)  lhs(-1,5)  lhs(-2,5) ...
  !           |    ...                                                                   
  !
  !    U is stored in band diagonal form 
  !     U = |   1    upper_1(1)  upper_2(1)      0           0           0           0
  !         |   0        1       upper_1(2)  upper_2(2)      0           0           0
  !         |   0        0           1       upper_1(3)  upper_2(3)      0           0
  !         |   0        0           0           1       upper_1(4)  upper_2(4)      0
  !         |   0        0           0           0           1       upper_1(5)  upper_2(5)  ...
  !         |  ...    
  !
  !   L is also stored in band diagonal form, but the lowest most band is equivalent to the 
  !   lowermost band of LHS, thus we don't need to store it
  !     L = | l_diag(1)        0            0           0             0       0   
  !         | lower_1(2)    l_diag(2)       0           0             0       0    
  !         |   l_2(3)      lower_1(3)   l_diag(3)      0             0       0          
  !         |     0           l_2(4)     lower_1(4)  l_diag(4)        0       0     
  !         |     0             0        l_2(5)      lower_1(5)    l_diag(5)  0   ...
  !         |  ...   
  ! 
  !
  !   To perform the LU decomposition, we go element by element. 
  !   First we start by noting that we want lhs=LU, so the first step of calculating
  !   L*U, by multiplying the first row of L by the columns of U, gives us
  !
  !     l_diag(1)*1          = lhs( 0,1)  =>  l_diag(1)  = lhs( 0,1)
  !     l_diag(1)*upper_1(1) = lhs(-1,1)  =>  upper_1(1) = lhs(-1,1) / l_diag(1)
  !     l_diag(1)*upper_2(1) = lhs(-2,1)  =>  upper_2(1) = lhs(-2,1) / l_diag(1)
  !
  !   Multiplying the second row of L by U now we get
  !     
  !     lower_1(2)*1                               = lhs(1,2)   =>  lower_1(2)  = lhs(1,2)
  !     lower_1(2)*upper_1(1)+l_diag(2)*1          = lhs(0,2)   =>  l_diag(2)   = lhs(0,2) - lower_1(2)*upper_1(1)
  !     lower_1(2)*upper_2(1)+l_diag(2)*upper_1(2) = lhs(-1,2)  =>  upper_1(2)  = ( lhs(-1,2)-lower_1(2)*upper_2(1) )
  !                                                                               / l_diag(2)
  !     l_diag(2)*upper_2(2)                       = lhs(-2,2)  =>  upper_2(2)  = lhs(-2,2) / l_diag(2)
  !
  !   Now that we're passed the k=1 and k=2 steps, each following step uses all the bands,
  !   allowing us to write the general step
  !     
  !     l_2(k)*1                                    = lhs(2,k)   =>  l_2(k)     = lhs(2,k)
  !     l_2(k)*upper_1(k-2)+lower_1(k)*1            = lhs(1,k)   =>  lower_1(k) = lhs(1,k) - l_2(k)*upper_1(k-2)
  !     l_2(k)*upper_2(k-2)+lower_1(k)*upper_1(k-1) = lhs( 0,k)  =>  l_diag(k)  = lhs(0,k) - l_2(k)*upper_2(k-2)
  !                    +l_diag(k)*1                                               + lower_1(k)*upper_1(k-1) 
  !                    
  !     lower_1(k)*upper_2(k-1)
  !     + l_diag(k)*upper_1(k) = lhs(-1,k)  =>  upper_1(k) = ( lhs(-1,k) - lower_1(k)*upper_2(k-1) )
  !                                                          / l_diag(k)
  !     l_diag(k)*upper_2(k)   = lhs(-2,k)  =>  upper_2(k) = lhs(-2,k) / l_diag(k)
  !
  !
  !   This general step is done for k from 3 to ndim-2 (do k = 3, ndim-2), and the last two
  !   steps are tweaked similarly to the first two, where we disclude one then two bands
  !   since they become no longer relevant. Note from this general step that the l_2 band
  !   is always equivalent to second subdiagonal band of lhs, thus we do not need to 
  !   calculate or store l_2. Also note that we only ever need l_diag so that we can divide 
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

  public :: penta_lu_solve

  private :: penta_lu_solve_single_rhs_multiple_lhs, penta_lu_solve_multiple_rhs_lhs

  interface penta_lu_solve
    module procedure penta_lu_solve_single_rhs_multiple_lhs
    module procedure penta_lu_solve_multiple_rhs_lhs
  end interface

  private ! Default scope

  contains

  !=============================================================================
  subroutine penta_lu_solve_single_rhs_multiple_lhs( ndim, ngrdcol, lhs, rhs, &
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
    real( kind = core_rknd ), intent(inout), dimension(-2:2,ngrdcol,ndim) :: &
      lhs   ! Matrices to solve, stored using band diagonal vectors
            ! -2 is the uppermost band, 2 is the lower most band, 0 is diagonal

    ! ----------------------- Output Variables -----------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,ndim) ::  &
      soln     ! Solution vector

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  &
      upper_1,          & ! First U band 
      upper_2,          & ! Second U band 
      lower_diag_invrs, & ! Inverse of the diagonal of L
      lower_1,          & ! First L band
      lower_2             ! Second L band

    integer :: i, k, j    ! Loop variables

    ! ----------------------- Begin Code -----------------------
       
    !$acc data create( upper_1, upper_2, lower_1, lower_2, lower_diag_invrs ) &
    !$acc      copyin( rhs, lhs ) &
    !$acc      copyout( soln )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      lower_diag_invrs(i,1) = 1.0_core_rknd / lhs(0,i,1)
      upper_1(i,1)          = lower_diag_invrs(i,1) * lhs(-1,i,1) 
      upper_2(i,1)          = lower_diag_invrs(i,1) * lhs(-2,i,1) 

      lower_1(i,2)          = lhs(1,i,2)
      lower_diag_invrs(i,2) = 1.0_core_rknd / ( lhs(0,i,2) - lower_1(i,2) * upper_1(i,1) )
      upper_1(i,2)          = lower_diag_invrs(i,2) * ( lhs(-1,i,2) - lower_1(i,2) * upper_2(i,1) )
      upper_2(i,2)          = lower_diag_invrs(i,2) * lhs(-2,i,2)
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      do k = 3, ndim-2
        lower_2(i,k) = lhs(2,i,k)
        lower_1(i,k) = lhs(1,i,k) - lower_2(i,k) * upper_1(i,k-2)

        lower_diag_invrs(i,k) = 1.0_core_rknd / ( lhs(0,i,k) - lower_2(i,k) * upper_2(i,k-2) &
                                                             - lower_1(i,k) * upper_1(i,k-1) )

        upper_1(i,k) = lower_diag_invrs(i,k) * ( lhs(-1,i,k) - lower_1(i,k) * upper_2(i,k-1) )
        upper_2(i,k) = lower_diag_invrs(i,k) * lhs(-2,i,k)
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      lower_2(i,ndim-1) = lhs(2,i,ndim-1)
      lower_1(i,ndim-1) = lhs(1,i,ndim-1) - lower_2(i,ndim-1) * upper_1(i,ndim-3)

      lower_diag_invrs(i,ndim-1) = 1.0_core_rknd  &
                                   / ( lhs(0,i,ndim-1) - lower_2(i,ndim-1) * upper_2(i,ndim-3) &
                                                       - lower_1(i,ndim-1) * upper_1(i,ndim-2) )

      upper_1(i,ndim-1)  = lower_diag_invrs(i,ndim-1) * ( lhs(-1,i,ndim-1) - lower_1(i,ndim-1) &
                                                                             * upper_2(i,ndim-2) )

      lower_2(i,ndim) = lhs(2,i,ndim)
      lower_1(i,ndim) = lhs(1,i,ndim) - lower_2(i,ndim) * upper_1(i,ndim-2)

      lower_diag_invrs(i,ndim) = 1.0_core_rknd  &
                                 / ( lhs(0,i,ndim-1) - lower_2(i,ndim) * upper_2(i,ndim-2) &
                                                     - lower_1(i,ndim) * upper_1(i,ndim-1) )
    end do
    !$acc end parallel loop
    
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol 

      soln(i,1)   = lower_diag_invrs(i,1) * rhs(i,1) 

      soln(i,2)   = lower_diag_invrs(i,2) * ( rhs(i,2) - lower_1(i,2) * soln(i,1) )

      do k = 3, ndim
        soln(i,k) = lower_diag_invrs(i,k) * ( rhs(i,k) - lower_2(i,k) * soln(i,k-2) &
                                                       - lower_1(i,k) * soln(i,k-1) )
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol 
      soln(i,ndim-1) = soln(i,ndim-1) - upper_1(i,ndim-1) * soln(i,ndim)

      do k = ndim-2, 1, -1
        soln(i,k) = soln(i,k) - upper_1(i,k) * soln(i,k+1) - upper_2(i,k) * soln(i,k+2)
      end do

    end do
    !$acc end parallel loop

    !$acc end data

  end subroutine penta_lu_solve_single_rhs_multiple_lhs

  
  !=============================================================================
  subroutine penta_lu_solve_multiple_rhs_lhs( ndim, nrhs, ngrdcol, lhs, rhs, &
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
    real( kind = core_rknd ), intent(inout), dimension(-2:2,ngrdcol,ndim) :: &
      lhs   ! Matrices to solve, stored using band diagonal vectors
            ! -2 is the uppermost band, 2 is the lower most band, 0 is diagonal

    ! ----------------------- Output Variables -----------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,ndim,nrhs) ::  &
      soln     ! Solution vector

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  &
      upper_1,          & ! First U band 
      upper_2,          & ! Second U band 
      lower_diag_invrs, & ! Inverse of the diagonal of L
      lower_1,          & ! First L band
      lower_2             ! Second L band

    integer :: i, k, j    ! Loop variables

    ! ----------------------- Begin Code -----------------------
       
    !$acc data create( upper_1, upper_2, lower_1, lower_2, lower_diag_invrs ) &
    !$acc      copyin( rhs, lhs ) &
    !$acc      copyout( soln )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      lower_diag_invrs(i,1) = 1.0_core_rknd / lhs(0,i,1)
      upper_1(i,1)          = lower_diag_invrs(i,1) * lhs(-1,i,1) 
      upper_2(i,1)          = lower_diag_invrs(i,1) * lhs(-2,i,1) 

      lower_1(i,2)          = lhs(1,i,2)
      lower_diag_invrs(i,2) = 1.0_core_rknd / ( lhs(0,i,2) - lower_1(i,2) * upper_1(i,1) )
      upper_1(i,2)          = lower_diag_invrs(i,2) * ( lhs(-1,i,2) - lower_1(i,2) * upper_2(i,1) )
      upper_2(i,2)          = lower_diag_invrs(i,2) * lhs(-2,i,2)
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      do k = 3, ndim-2
        lower_2(i,k) = lhs(2,i,k)
        lower_1(i,k) = lhs(1,i,k) - lower_2(i,k) * upper_1(i,k-2)

        lower_diag_invrs(i,k) = 1.0_core_rknd / ( lhs(0,i,k) - lower_2(i,k) * upper_2(i,k-2) &
                                                             - lower_1(i,k) * upper_1(i,k-1) )

        upper_1(i,k) = lower_diag_invrs(i,k) * ( lhs(-1,i,k) - lower_1(i,k) * upper_2(i,k-1) )
        upper_2(i,k) = lower_diag_invrs(i,k) * lhs(-2,i,k)
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      lower_2(i,ndim-1) = lhs(2,i,ndim-1)
      lower_1(i,ndim-1) = lhs(1,i,ndim-1) - lower_2(i,ndim-1) * upper_1(i,ndim-3)

      lower_diag_invrs(i,ndim-1) = 1.0_core_rknd  &
                                   / ( lhs(0,i,ndim-1) - lower_2(i,ndim-1) * upper_2(i,ndim-3) &
                                                       - lower_1(i,ndim-1) * upper_1(i,ndim-2) )

      upper_1(i,ndim-1) = lower_diag_invrs(i,ndim-1) * ( lhs(-1,i,ndim-1) - lower_1(i,ndim-1) &
                                                                            * upper_2(i,ndim-2) )

      lower_2(i,ndim) = lhs(2,i,ndim)
      lower_1(i,ndim) = lhs(1,i,ndim) - lower_2(i,ndim) * upper_1(i,ndim-2)

      lower_diag_invrs(i,ndim) = 1.0_core_rknd  &
                                 / ( lhs(0,i,ndim-1) - lower_2(i,ndim) * upper_2(i,ndim-2) &
                                                     - lower_1(  i,ndim) * upper_1(i,ndim-1) )
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do j = 1, nrhs
      do i = 1, ngrdcol 

        soln(i,1,j)   = lower_diag_invrs(i,1) * rhs(i,1,j) 

        soln(i,2,j)   = lower_diag_invrs(i,2) * ( rhs(i,2,j) - lower_1(i,2) * soln(i,1,j) )

        do k = 3, ndim
          soln(i,k,j) = lower_diag_invrs(i,k) * ( rhs(i,k,j) - lower_2(i,k) * soln(i,k-2,j) &
                                                             - lower_1(i,k) * soln(i,k-1,j) )
        end do
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do j = 1, nrhs
      do i = 1, ngrdcol 
        soln(i,ndim-1,j) = soln(i,ndim-1,j) - upper_1(i,ndim-1) * soln(i,ndim,j)

        do k = ndim-2, 1, -1
          soln(i,k,j) = soln(i,k,j) - upper_1(i,k) * soln(i,k+1,j) - upper_2(i,k) * soln(i,k+2,j)
        end do

      end do
    end do
    !$acc end parallel loop

    !$acc end data

  end subroutine penta_lu_solve_multiple_rhs_lhs

end module penta_lu_solvers
