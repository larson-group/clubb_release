module penta_lu_solvers

  ! Description:
  !   These routines solve lhs*sol=rhs using LU decomp. 
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
  !     U = |   1    u_1(1)  u_2(1)    0       0       0       0
  !         |   0      1     u_1(2)  u_2(2)    0       0       0
  !         |   0      0       1     u_1(3)  u_2(3)    0       0
  !         |   0      0       0       1     u_1(4)  u_2(4)    0
  !         |   0      0       0       0       1     u_1(5)  u_2(5)  ...
  !         |  ...   
  !
  !   L is also stored in band diagonal form, but the lowest most band is equivalent to the 
  !   lowermost band of LHS, thus we don't need to store it
  !     L = | l_diag(1)       0         0           0         0       0   
  !         |  l_1(2)    l_diag(2)      0           0         0       0    
  !         |  l_2(3)     l_1(3)    l_diag(3)       0         0       0          
  !         |    0        l_2(4)     l_1(4)     l_diag(4)     0       0     
  !         |    0          0        l_2(5)      l_1(5)    l_diag(5)  0   ...
  !         |  ...   
  ! 
  !
  !   To perform the LU decomposition, we go element by element. 
  !   First we start by noting that we want lhs=LU, so the first step of calculating
  !   L*U, by multiplying the first row of L by the columns of U, gives us
  !
  !     l_diag(1)*1      = lhs( 0,1)  =>  l_diag(1) = lhs( 0,1)
  !     l_diag(1)*u_1(1) = lhs(-1,1)  =>  u_1(1)    = lhs(-1,1) / l_diag(1)
  !     l_diag(1)*u_2(1) = lhs(-2,1)  =>  u_2(1)    = lhs(-2,1) / l_diag(1)
  !
  !   Multiplying the second row of L by U now we get
  !     
  !     l_1(2)*1                       = lhs(1,2)   =>  l_1(2)    = lhs(1,2)
  !     l_1(2)*u_1(1)+l_diag(2)*1      = lhs(0,2)   =>  l_diag(2) = lhs(0,2) - l_1(2)*u_1(1)
  !     l_1(2)*u_2(1)+l_diag(2)*u_1(2) = lhs(-1,2)  =>  u_1(2)    = ( lhs(-1,2)-l_1(2)*u_2(1) )
  !                                                                 / l_diag(2)
  !     l_diag(2)*u_2(2)               = lhs(-2,2)  =>  u_2(2)    = lhs(-2,2) / l_diag(2)
  !
  !   Now that we're passed the k=1 and k=2 steps, each following step uses all the bands,
  !   allowing us to write the general step
  !     
  !     l_2(k)*1                        = lhs(2,k)   =>  l_2(k)    = lhs(2,k)
  !     l_2(k)*u_1(k-2)+l_1(k)*1        = lhs(1,k)   =>  l_1(k)    = lhs(1,k) - l_2(k)*u_1(k-2)
  !     l_2(k)*u_2(k-2)+l_1(k)*u_1(k-1) = lhs( 0,k)  =>  l_diag(k) = lhs(0,k) - l_2(k)*u_2(k-2)
  !                    +l_diag(k)*1                                           + l_1(k)*u_1(k-1) 
  !                    
  !     l_1(k)*u_2(k-1)+l_diag(k)*u_1(k) = lhs(-1,k)  =>  u_1(k) = ( lhs(-1,k) - l_1(k)*u_2(k-1) )
  !                                                                / l_diag(k)
  !     l_diag(k)*u_2(k)                 = lhs(-2,k)  =>  u_2(k) = lhs(-2,k) / l_diag(k)
  !
  !
  !   This general step is done for k from 3 to ndim-2 (do k = 3, ndim-2), and the last two
  !   steps are tweaked similarly to the first two, where we disclude one then two bands
  !   since they become no longer relevant. Note from this general step that the l_2 band
  !   is always equivalent to second subdiagonal band of lhs, thus we do not need to 
  !   calculate or store l_2. Also note that we only ever need l_diag so that we can divide 
  !   by it, so instead we compute l_diag_invrs to reduce divide operations.
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

  public :: penta_lu_solve, penta_lu_solve_batched, penta_lu_solve_batched_gpu

  interface penta_lu_solve
#ifdef _OPENACC && CLUBB_GPU
    module procedure penta_lu_solve_batched_gpu
    module procedure penta_lu_solve_single_rhs_gpu
#else
    module procedure penta_lu_solve_batched
    module procedure penta_lu_solve_single_rhs
#endif
  end interface

  private ! Default scope

  contains

  !=============================================================================
  subroutine penta_lu_solve_single_rhs( ndim, ngrdcol, lhs, rhs, &
                                        sol )
    ! Description:
    !   Written for a single RHS and multiple LHS, optimized for CPUs.
    !------------------------------------------------------------------------

    implicit none

    ! ----------------------- Input Variables -----------------------
    integer, intent(in) :: &
      ndim,   & ! Matrix size 
      ngrdcol   ! Number of grid columns

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,ndim) ::  &
      rhs       ! Right hand side 

    ! ----------------------- Input/Output Variables -----------------------
    real( kind = core_rknd ), intent(inout), dimension(-2:2,ngrdcol,ndim) :: &
      lhs       ! Matrices to solve, stored using band diagonal vectors
                ! -2 is the uppermost band, 2 is the lower most band, 0 is diagonal

    ! ----------------------- Output Variables -----------------------  
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,ndim) ::  &
      sol       ! Solution vector

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  &
      u_1, & ! First U band 
      u_2    ! Second U band 

    real( kind = core_rknd ) ::  &
      l_diag_invrs, & ! Inverse of the diagonal of L
      l_1             ! First L band

    integer :: i, k, j    ! Loop variables

    ! ----------------------- Begin Code -----------------------

    do i = 1, ngrdcol
      l_diag_invrs = 1.0_core_rknd / lhs(0,i,1)
      sol(i,1)     = l_diag_invrs * rhs(i,1) 
      u_1(i,1)     = l_diag_invrs * lhs(-1,i,1) 
      u_2(i,1)     = l_diag_invrs * lhs(-2,i,1) 
    end do

    do i = 1, ngrdcol
      l_1          = lhs(1,i,2)
      l_diag_invrs = 1.0_core_rknd  / ( lhs(0,i,2) - l_1 * u_1(i,1) )
      sol(i,2)     = l_diag_invrs * ( rhs(i,2) - l_1 * sol(i,1) ) 
      u_1(i,2)     = l_diag_invrs * ( lhs(-1,i,2) - l_1 * u_2(i,1) ) 
      u_2(i,2)     = l_diag_invrs * lhs(-2,i,2) 
    end do
    
    do k = 3, ndim-2
      do i = 1, ngrdcol
        l_1      = lhs(1,i,k) - lhs(2,i,k) * u_1(i,k-2)
        l_diag_invrs = 1.0_core_rknd  / ( lhs(0,i,k) - lhs(2,i,k) * u_2(i,k-2) &
                                                     - l_1 *u_1(i,k-1) )
        sol(i,k) = l_diag_invrs * ( rhs(i,k) - lhs(2,i,k) * sol(i,k-2) - l_1 * sol(i,k-1) )
        u_1(i,k) = l_diag_invrs * ( lhs(-1,i,k) - l_1 * u_2(i,k-1) ) 
        u_2(i,k) = l_diag_invrs * lhs(-2,i,k) 
      end do
    end do
    
    do i = 1, ngrdcol
      l_1 = lhs(1,i,ndim-1) - lhs(2,i,ndim-1) * u_1(i,ndim-3)
      l_diag_invrs  = 1.0_core_rknd  / ( lhs(0,i,ndim-1) - lhs(2,i,ndim-1) * u_2(i,ndim-3) &
                                                         - l_1 * u_1(i,ndim-2) )
      sol(i,ndim-1) = l_diag_invrs * ( rhs(i,ndim-1) - lhs(2,i,ndim-1) * sol(i,ndim-3) &
                                                     - l_1 * sol(i,ndim-2) ) 
      u_1(i,ndim-1) = l_diag_invrs * ( lhs(-1,i,ndim-1) - l_1 * u_2(i,ndim-2) ) 
    end do

    do i = 1, ngrdcol
      l_1 = lhs(1,i,ndim) - lhs(2,i,ndim) * u_1(i,ndim-2)
      l_diag_invrs = 1.0_core_rknd  / ( lhs(0,i,ndim-1) - lhs(2,i,ndim) * u_2(i,ndim-2) &
                                                        - l_1 * u_1(i,ndim-1) )
      sol(i,ndim)  = l_diag_invrs * ( rhs(i,ndim) - lhs(2,i,ndim) * sol(i,ndim-2) &
                                                  - l_1 * sol(i,ndim-1) ) 
    end do

    do i = 1, ngrdcol
      sol(i,ndim-1) = sol(i,ndim-1) - u_1(i,ndim-1) * sol(i,ndim)
    end do

    do k = ndim-2, 1, -1
      do i = 1, ngrdcol
        sol(i,k) = sol(i,k) - u_1(i,k) * sol(i,k+1) - u_2(i,k) * sol(i,k+2)
      end do
    end do

  end subroutine penta_lu_solve_single_rhs

  !=============================================================================
  subroutine penta_lu_solve_batched( ndim, nrhs, ngrdcol, lhs, rhs, &
                                     sol )
    ! Description:
    !   Written for multiple RHS and multiple LHS, optimized to run 
    !   on CPUs,
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
      sol     ! Solution vector

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  &
      u_1, & ! First U band 
      u_2    ! Second U band 

    real( kind = core_rknd ) ::  &
      l_diag_invrs, & ! Inverse of the diagonal of L
      l_1             ! First L band

    integer :: i, k, j    ! Loop variables

    ! ----------------------- Begin Code -----------------------
       
    do i = 1, ngrdcol
      l_diag_invrs = 1.0_core_rknd / lhs(0,i,1)
      sol(i,1,:)   = l_diag_invrs * rhs(i,1,:) 
      u_1(i,1)     = l_diag_invrs * lhs(-1,i,1) 
      u_2(i,1)     = l_diag_invrs * lhs(-2,i,1) 
    end do

    do i = 1, ngrdcol
      l_1          = lhs(1,i,2)
      l_diag_invrs = 1.0_core_rknd  / ( lhs(0,i,2) - l_1 * u_1(i,1) )
      sol(i,2,:)   = l_diag_invrs * ( rhs(i,2,:) - l_1 * sol(i,1,:) )
      u_1(i,2)     = l_diag_invrs * ( lhs(-1,i,2) - l_1 * u_2(i,1) )
      u_2(i,2)     = l_diag_invrs * lhs(-2,i,2)
    end do
    
    do k = 3, ndim-2
      do i = 1, ngrdcol 
        l_1          = lhs(1,i,k) - lhs(2,i,k) * u_1(i,k-2)
        l_diag_invrs = 1.0_core_rknd  / ( lhs(0,i,k) - lhs(2,i,k) * u_2(i,k-2) - l_1 *u_1(i,k-1) )
        sol(i,k,:)   = l_diag_invrs * ( rhs(i,k,:) - lhs(2,i,k) * sol(i,k-2,:) &
                                                   - l_1 * sol(i,k-1,:) )
        u_1(i,k)     = l_diag_invrs * ( lhs(-1,i,k) - l_1 * u_2(i,k-1) )
        u_2(i,k)     = l_diag_invrs * lhs(-2,i,k)
      end do
    end do
    
    do i = 1, ngrdcol
      l_1 = lhs(1,i,ndim-1) - lhs(2,i,ndim-1) * u_1(i,ndim-3)
      l_diag_invrs    = 1.0_core_rknd  / ( lhs(0,i,ndim-1) - lhs(2,i,ndim-1) * u_2(i,ndim-3) &
                                                           - l_1 * u_1(i,ndim-2) )
      sol(i,ndim-1,:) = l_diag_invrs * ( rhs(i,ndim-1,:) - lhs(2,i,ndim-1) * sol(i,ndim-3,:) &
                                                         - l_1 * sol(i,ndim-2,:) ) 
      u_1(i,ndim-1)   = l_diag_invrs * ( lhs(-1,i,ndim-1) - l_1 * u_2(i,ndim-2) )
    end do

    do i = 1, ngrdcol 
      l_1 = lhs(1,i,ndim) - lhs(2,i,ndim) * u_1(i,ndim-2)
      l_diag_invrs   = 1.0_core_rknd  / ( lhs(0,i,ndim-1) - lhs(2,i,ndim) * u_2(i,ndim-2) &
                                                          - l_1 * u_1(i,ndim-1) )
      sol(i,ndim,:)  = l_diag_invrs * ( rhs(i,ndim,:) - lhs(2,i,ndim) * sol(i,ndim-2,:) &
                                                      - l_1 * sol(i,ndim-1,:) ) 
    end do

    do i = 1, ngrdcol
      sol(i,ndim-1,:) = sol(i,ndim-1,:) - u_1(i,ndim-1) * sol(i,ndim,:)
    end do

    do j = 1, nrhs 
      do k = ndim-2, 1, -1
        do i = 1, ngrdcol
          sol(i,k,j) = sol(i,k,j) - u_1(i,k) * sol(i,k+1,j) - u_2(i,k) * sol(i,k+2,j)
        end do
      end do
    end do

  end subroutine penta_lu_solve_batched



  !=============================================================================
  !                               GPU VERSIONS
  !=============================================================================


  subroutine penta_lu_solve_batched_gpu( ndim, nrhs, ngrdcol, lhs, rhs, &
                                         sol )
    ! Description:
    !   Written for multiple RHS and multiple LHS, optimized 
    !   to run on GPUs.
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
      sol     ! Solution vector

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  &
      u_1, & ! First U band 
      u_2, & ! Second U band 
      l_diag_invrs, & ! Inverse of the diagonal of L
      l_1             ! First L band

    integer :: i, k, j    ! Loop variables

    ! ----------------------- Begin Code -----------------------
       
    !$acc data create( u_1, u_2, l_1, l_diag_invrs ) &
    !$acc      copyin( rhs, lhs ) &
    !$acc      copyout( sol )

    !$acc parallel loop default(present) async(1)
    do i = 1, ngrdcol
      l_diag_invrs(i,1) = 1.0_core_rknd / lhs(0,i,1)
      u_1(i,1)          = l_diag_invrs(i,1) * lhs(-1,i,1) 
      u_2(i,1)          = l_diag_invrs(i,1) * lhs(-2,i,1) 

      l_1(i,2)          = lhs(1,i,2)
      l_diag_invrs(i,2) = 1.0_core_rknd / ( lhs(0,i,2) - l_1(i,2) * u_1(i,1) )
      u_1(i,2)          = l_diag_invrs(i,2) * ( lhs(-1,i,2) - l_1(i,2) * u_2(i,1) )
      u_2(i,2)          = l_diag_invrs(i,2) * lhs(-2,i,2)
    end do

    !$acc parallel loop default(present) async(1)
    do i = 1, ngrdcol 
      do k = 3, ndim-2
        l_1(i,k)          = lhs(1,i,k) - lhs(2,i,k) * u_1(i,k-2)

        l_diag_invrs(i,k) = 1.0_core_rknd / ( lhs(0,i,k) - lhs(2,i,k) * u_2(i,k-2) &
                                                         - l_1(  i,k) * u_1(i,k-1) )

        u_1(i,k)          = l_diag_invrs(i,k) * ( lhs(-1,i,k) - l_1(i,k) * u_2(i,k-1) )
        u_2(i,k)          = l_diag_invrs(i,k) * lhs(-2,i,k)
      end do
    end do

    !$acc parallel loop default(present) async(1)
    do i = 1, ngrdcol
      l_1(i,ndim-1) = lhs(1,i,ndim-1) - lhs(2,i,ndim-1) * u_1(i,ndim-3)

      l_diag_invrs(i,ndim-1) = 1.0_core_rknd  &
                               / ( lhs(0,i,ndim-1) - lhs(2,i,ndim-1) * u_2(i,ndim-3) &
                                                   - l_1(  i,ndim-1) * u_1(i,ndim-2) )

      u_1(i,ndim-1)   = l_diag_invrs(i,ndim-1) * ( lhs(-1,i,ndim-1) - l_1(i,ndim-1) * u_2(i,ndim-2) )

      l_1(i,ndim) = lhs(1,i,ndim) - lhs(2,i,ndim) * u_1(i,ndim-2)

      l_diag_invrs(i,ndim)   = 1.0_core_rknd  &
                               / ( lhs(0,i,ndim-1) - lhs(2,i,ndim) * u_2(i,ndim-2) &
                                                   - l_1(  i,ndim) * u_1(i,ndim-1) )
    end do

    
    !$acc parallel loop collapse(2) default(present) async(1)
    do j = 1, nrhs
      do i = 1, ngrdcol 

        sol(i,1,j)   = l_diag_invrs(i,1) * rhs(i,1,j) 

        sol(i,2,j)   = l_diag_invrs(i,2) * ( rhs(i,2,j) - l_1(i,2) * sol(i,1,j) )

        do k = 3, ndim
          sol(i,k,j)   = l_diag_invrs(i,k) * ( rhs(i,k,j) - lhs(2,i,k) * sol(i,k-2,j) &
                                                          - l_1(  i,k) * sol(i,k-1,j) )
        end do

        sol(i,ndim-1,j) = sol(i,ndim-1,j) - u_1(i,ndim-1) * sol(i,ndim,j)

        do k = ndim-2, 1, -1
          sol(i,k,j) = sol(i,k,j) - u_1(i,k) * sol(i,k+1,j) - u_2(i,k) * sol(i,k+2,j)
        end do

      end do
    end do

    !$acc wait(1)

    !$acc end data

  end subroutine penta_lu_solve_batched_gpu


  !=============================================================================
  subroutine penta_lu_solve_single_rhs_gpu( ndim, ngrdcol, lhs, rhs, &
                                            sol )
    ! Description:
    !   Written for a single RHS and multiple LHS, optimized 
    !   to run on GPUs.
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
      sol     ! Solution vector

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  &
      u_1, & ! First U band 
      u_2, & ! Second U band 
      l_diag_invrs, & ! Inverse of the diagonal of L
      l_1             ! First L band

    integer :: i, k, j    ! Loop variables

    ! ----------------------- Begin Code -----------------------
       
    !$acc data create( u_1, u_2, l_1, l_diag_invrs ) &
    !$acc      copyin( rhs, lhs ) &
    !$acc      copyout( sol )

    !$acc parallel loop default(present) async(1)
    do i = 1, ngrdcol
      l_diag_invrs(i,1) = 1.0_core_rknd / lhs(0,i,1)
      u_1(i,1)          = l_diag_invrs(i,1) * lhs(-1,i,1) 
      u_2(i,1)          = l_diag_invrs(i,1) * lhs(-2,i,1) 

      l_1(i,2)          = lhs(1,i,2)
      l_diag_invrs(i,2) = 1.0_core_rknd / ( lhs(0,i,2) - l_1(i,2) * u_1(i,1) )
      u_1(i,2)          = l_diag_invrs(i,2) * ( lhs(-1,i,2) - l_1(i,2) * u_2(i,1) )
      u_2(i,2)          = l_diag_invrs(i,2) * lhs(-2,i,2)

      do k = 3, ndim-2

        l_1(i,k)          = lhs(1,i,k) - lhs(2,i,k) * u_1(i,k-2)

        l_diag_invrs(i,k) = 1.0_core_rknd / ( lhs(0,i,k) - lhs(2,i,k) * u_2(i,k-2) &
                                                         - l_1(  i,k) * u_1(i,k-1) )

        u_1(i,k)          = l_diag_invrs(i,k) * ( lhs(-1,i,k) - l_1(i,k) * u_2(i,k-1) )
        u_2(i,k)          = l_diag_invrs(i,k) * lhs(-2,i,k)
      end do

      l_1(i,ndim-1) = lhs(1,i,ndim-1) - lhs(2,i,ndim-1) * u_1(i,ndim-3)

      l_diag_invrs(i,ndim-1) = 1.0_core_rknd  &
                               / ( lhs(0,i,ndim-1) - lhs(2,i,ndim-1) * u_2(i,ndim-3) &
                                                   - l_1(  i,ndim-1) * u_1(i,ndim-2) )

      u_1(i,ndim-1)   = l_diag_invrs(i,ndim-1) * ( lhs(-1,i,ndim-1) - l_1(i,ndim-1) * u_2(i,ndim-2) )

      l_1(i,ndim) = lhs(1,i,ndim) - lhs(2,i,ndim) * u_1(i,ndim-2)

      l_diag_invrs(i,ndim)   = 1.0_core_rknd  &
                               / ( lhs(0,i,ndim-1) - lhs(2,i,ndim) * u_2(i,ndim-2) &
                                                   - l_1(  i,ndim) * u_1(i,ndim-1) )
    end do

    
    !$acc parallel loop default(present) async(1)
    do i = 1, ngrdcol 

      sol(i,1)   = l_diag_invrs(i,1) * rhs(i,1) 

      sol(i,2)   = l_diag_invrs(i,2) * ( rhs(i,2) - l_1(i,2) * sol(i,1) )

      do k = 3, ndim
        sol(i,k)   = l_diag_invrs(i,k) * ( rhs(i,k) - lhs(2,i,k) * sol(i,k-2) &
                                                    - l_1(  i,k) * sol(i,k-1) )
      end do

      sol(i,ndim-1) = sol(i,ndim-1) - u_1(i,ndim-1) * sol(i,ndim)

      do k = ndim-2, 1, -1
        sol(i,k) = sol(i,k) - u_1(i,k) * sol(i,k+1) - u_2(i,k) * sol(i,k+2)
      end do

    end do

    !$acc wait(1)

    !$acc end data

  end subroutine penta_lu_solve_single_rhs_gpu

end module penta_lu_solvers
