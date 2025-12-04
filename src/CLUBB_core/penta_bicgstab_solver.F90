module bicgstab_solvers

  ! Description:
  !   Solve lhs*soln=rhs using the biconjugate gradient stabilized method. 
  !   This is an iterative solver, and mainly experimental. The convergence
  !   properties and error should be better analyzed before real use.
  !
  !   The benefit of this method is that uses almost entirely matrix-vector 
  !   calculations, which are highly parallelizable, and thus "well-suited" 
  !   for use on GPUs. However, benchmarks show this to be slower than the 
  !   direct LU solver in penta_lu_solver.F90, even on GPUs.
  !
  ! Notes: 
  !
  !   CURRENT VERSION NOT GPUIZED
  !   
  !   LHS needs to be stored in band diagonal form:
  !     lhs = | lhs(0,1)  lhs(-1,1)  lhs(-2,1)     0           0          0          0
  !           | lhs(1,2)  lhs( 0,2)  lhs(-1,2)  lhs(-2,2)      0          0          0
  !           | lhs(2,3)  lhs( 1,3)  lhs( 0,3)  lhs(-1,3)  lhs(-2,3)      0          0
  !           |     0     lhs( 2,4)  lhs( 1,4)  lhs( 0,4)  lhs(-1,4)  lhs(-2,4)      0
  !           |     0         0      lhs( 2,5)  lhs( 1,5)  lhs( 0,5)  lhs(-1,5)  lhs(-2,5) ...
  !           |    ...                                                                   
  !
  ! References:
  !   https://epubs.siam.org/doi/10.1137/0913035
  !   
  !   Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG for the 
  !              Solution of Nonsymmetric Linear Systems
  !   H. A. van der Vorst
  !
  !   Using the algorithm for the preconditioned Bi-CGSTAB on p. 638
  !------------------------------------------------------------------------

  use clubb_precision, only:  &
    core_rknd ! Variable(s)

  use constants_clubb, only:  &
    zero      ! Parameters(s)

  implicit none
  
  integer, parameter :: &
    n_pc_bands = 5            ! Number of bands to use in the preconditioner
                            ! The best number for this is unclear, it seems to 
                            ! depend on the matrix size, and almost certianly
                            ! depends on the number of lhs bands too. 
                            ! 5 seems like a good number, but this needs more investigation

  real( kind = core_rknd ), parameter :: &
    bicgstab_tol = max( 1.e-10, epsilon(bicgstab_tol) )   ! Reference tolerance used for convergence criteria.
                                                          ! This is somewhat arbitrary, feel free to experiment with it.
                                                          ! In theory, lower tolerance = less error, but too small
                                                          ! of a tolerance result in unnecessary iterations. 
                                                          ! max-epsilon statement is used to increase tolerance 
                                                          ! when using single precision

  public  :: penta_bicgstab_solve

  private ! Default scope

  contains

  !=============================================================================
  subroutine penta_bicgstab_solve( ndim, ngrdcol, lhs, rhs, &
                                   soln, &
                                   sol_init, &
                                   iters )
    ! Description:
    !   Written for one RHS and one LHS per grid column.
    !
    !   ndim      - matrix dimension
    !   ngrdcol   - number of systems to solve (one for each grid column)
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
      soln       ! Solution vector

    ! ----------------------- Optional Input Variables -----------------------
    real( kind = core_rknd ), intent(in), optional, dimension(ngrdcol,ndim) ::  &
      sol_init  ! Initial guess for solution

    integer, intent(out), optional ::  &
      iters  ! Number of iterations used

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  &
      residual, & ! Residual
      s_res,    & ! Residual after Bi-CG step
      res_hat,  & ! Saved initial residual
      p_search, & ! Search direction for next step
      Ay,       & ! A*s_res, A being our matrix stored in bands
      Az,       & ! A*p_search, A being our matrix stored in bands
      z_vec, &
      y_vec    

    real( kind = core_rknd ), dimension(ngrdcol) ::  &
      omega,            & ! Scale factor used in updating p_search and soln
      alpha,            & ! Scale factor used in updating s_res
      beta,             & ! Scale factor used in updating p_search
      res_dot_res_hat,  & ! Dot product of residual and res_hat
      res_tol,          & ! Tolerance used for stopping condition
      error               ! Error for every grid columns

    real( kind = core_rknd ) ::  &
      norm_rhs,             & ! Norm of each rhs
      res_dot_res_hat_old     ! Old value of dot product of residual and res_hat, temporary value

    integer :: iter, j, k, i    ! Loop Variables

    ! ----------------------- Begin Code -----------------------

    ! Calculate res_tol (error tolerance) using rhs norm.
    ! Some of our matrices solve for coupled fields of very different magnitudes,
    ! so it may make more sense to have a per-vertical level error tolerance
    ! which will be more expensive, but not a big deal for a GPU, which this
    ! solver is intended for. Another option could be to track the maximum % change
    ! in the residual, and consider that as the error instead. 
    do i = 1, ngrdcol
      norm_rhs = sqrt(sum(rhs(i,:)*rhs(i,:)))

      if ( norm_rhs == 0.0) then
        res_tol(i) = bicgstab_tol
      else
        res_tol(i) = bicgstab_tol*norm_rhs
      end if  
      
    end do

    ! Set initial solution if preset, otherwise just use 0
    if ( present(sol_init) ) then
      do k = 1, ndim
        do i = 1, ngrdcol
          soln(i,k) = sol_init(i,k)
        end do
      end do
    else
      soln = zero
    end if

    ! Calculate initial residual: residual = rhs - A * soln
    call penta_matmul_batched( ngrdcol,ndim, &
                               lhs, soln, &
                               residual )

    ! Finish residual calculation and set res_hat (which doesn't change) 
    ! and use this first residual vector as the first search vector
    do k = 1, ndim
      do i = 1, ngrdcol
        residual(i,k) = rhs(i,k) - residual(i,k)
        res_hat(i,k)  = residual(i,k)
        p_search(i,k) = residual(i,k)
      end do
    end do

    ! Initial error and vectors
    do i = 1, ngrdcol
      res_dot_res_hat(i)  = dot_product( residual(i,:), residual(i,:) )
      error(i)            = sqrt( res_dot_res_hat(i) / ndim )
    end do

    iter = 0

    ! Loop until error small, cap at ndim since this is supposed to 
    ! be an exact solver after ndim iterations
    do while( any(error > res_tol)  .and. iter <= ndim )

      iter = iter + 1

      ! Solve for y_vec: Ky=p_search
      call precond_psuedo_solve( ngrdcol,ndim, &
                                 lhs,p_search, &
                                 y_vec )

      ! Ay = A * y_vec
      call penta_matmul_batched( ngrdcol,ndim, &
                                 lhs,y_vec, &
                                 Ay )

      ! alpha = rho / ( rho_hat, v )
      do i = 1, ngrdcol
        alpha(i) = res_dot_res_hat(i) / dot_product(res_hat(i,:), Ay(i,:))
      end do

      ! s = residual - alpha * v
      do k = 1, ndim
        do i = 1, ngrdcol
          s_res(i,k) = residual(i,k) - alpha(i) * Ay(i,k)
        end do
      end do


      ! Solve for z_vec: Kz = s_res
      call precond_psuedo_solve( ngrdcol,ndim, &
                                 lhs,s_res, &
                                 z_vec )

      ! Ay = A * y_vec
      call penta_matmul_batched( ngrdcol,ndim, &
                                 lhs,z_vec, &
                                 Az )
          
      ! omega = (t,s)/(t,t) or <K1^-1 t, K1^-1 s>/<K1^-1 t, K1^-1 t> in the offocial algorithm
      do i = 1, ngrdcol
        omega(i) =  dot_product( Az(i,:), s_res(i,:) ) / dot_product( Az(i,:), Az(i,:) )
      end do

      ! update solution
      do k = 1, ndim
        do i = 1, ngrdcol
          if ( error(i) > res_tol(i) ) then
            soln(i,k) = soln(i,k) + alpha(i) * y_vec(i,k) &
                                  + omega(i) * z_vec(i,k)
          end if
        end do
      end do

      ! recompute residual
      do k = 1, ndim
        do i = 1, ngrdcol
          residual(i,k) = s_res(i,k) - omega(i) * Az(i,k)
        end do
      end do

      ! recompute error
      do i = 1, ngrdcol
        error(i) = sqrt(dot_product( residual(i,:), residual(i,:) ) / ndim)
      end do

      ! update beta, recalculate
      do i = 1, ngrdcol
        res_dot_res_hat_old = res_dot_res_hat(i)
        res_dot_res_hat(i)  = dot_product(res_hat(i,:), residual(i,:))

        beta(i) = ( res_dot_res_hat(i) / res_dot_res_hat_old ) * ( alpha(i) / omega(i) )
      end do

      ! update search vector
      do k = 1, ndim
        do i = 1, ngrdcol
          p_search(i,k) = residual(i,k) &
                          + beta(i) * ( p_search(i,k) - omega(i) * Ay(i,k) )
        end do
      end do

    end do

    if ( present(iters) ) iters = iter

  end subroutine penta_bicgstab_solve

  subroutine penta_matmul_batched( ngrdcol,ndim, &
                                   lhs, x_vec, &
                                   y_vec )
    ! Description:
    !   This performs batches of matrix multiplications, 
    !   calculating y_vec = lhs * x_vec, for ngrdcol different systems,
    !   and assuming lhs is penta-diagonal.
    !------------------------------------

    implicit none

    integer, intent(in) :: ngrdcol,ndim

    real( kind = core_rknd ), intent(in), dimension(-2:2,ngrdcol,ndim) :: &
      lhs

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,ndim) ::  &
      x_vec

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,ndim) ::  &
      y_vec

    integer :: i, j, k


    do i = 1, ngrdcol
      y_vec(i,1) =   lhs( 0,i,1) * x_vec(i,1) &
                   + lhs(-1,i,1) * x_vec(i,2) &
                   + lhs(-2,i,1) * x_vec(i,3)

      y_vec(i,2) =   lhs( 1,i,2) * x_vec(i,1) &
                   + lhs( 0,i,2) * x_vec(i,2) &
                   + lhs(-1,i,2) * x_vec(i,3) &
                   + lhs(-2,i,2) * x_vec(i,4)
    end do

    do k = 3, ndim-2
      do i = 1, ngrdcol
        y_vec(i,k) =    lhs( 2,i,k) * x_vec(i,k-2) &
                      + lhs( 1,i,k) * x_vec(i,k-1) &
                      + lhs( 0,i,k) * x_vec(i,k  ) &
                      + lhs(-1,i,k) * x_vec(i,k+1) &
                      + lhs(-2,i,k) * x_vec(i,k+2)
      end do
    end do

    do i = 1, ngrdcol
      y_vec(i,ndim-1) =   lhs( 2,i,ndim-1) * x_vec(i,ndim-3) &
                        + lhs( 1,i,ndim-1) * x_vec(i,ndim-2) &
                        + lhs( 0,i,ndim-1) * x_vec(i,ndim-1) &
                        + lhs(-1,i,ndim-1) * x_vec(i,ndim  )

      y_vec(i,ndim)   =   lhs(2,i,ndim) * x_vec(i,ndim-2) &
                        + lhs(1,i,ndim) * x_vec(i,ndim-1) &
                        + lhs(0,i,ndim) * x_vec(i,ndim  )
    end do


  end subroutine penta_matmul_batched

  subroutine precond_psuedo_solve( ngrdcol,ndim, &
                                   lhs, x_vec, &
                                   y_vec )

    ! Approximately solve for y_vec in:  lhs * y_vec = x_vec
    ! The "solve" is done as follows:
    ! 1) Denote lhs = L+D+U, 
    !     where L, D, U are lower, diagonal, and upper parts of lhs
    ! 2) Find diagonal matrix E such that K = (L+E)E^-1(U+E) ~= A
    !     these lower and upper parts are easily invertible 
    ! 3) Compute y_vec = lhs^-1 x_vec ~= K^-1 x_vec = (U+E)^-1 E (L+E)^-1 x_vec

    implicit none

    integer, intent(in) :: ngrdcol,ndim

    real( kind = core_rknd ), intent(in), dimension(-2:2,ngrdcol,ndim) :: &
      lhs

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,ndim) ::  &
      x_vec

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,ndim) ::  &
      y_vec


    ! ----------------------- Local Variables -----------------------

    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  &
      E_inv

    real( kind = core_rknd ), dimension(ngrdcol,n_pc_bands,ndim) ::  &
      lower, upper

    real( kind = core_rknd ), dimension(ngrdcol,ndim) ::  y_tmp


    integer :: j, k, i, b    ! Loop Variables

    ! ----------------------- Begin Code -----------------------

    lower = zero
    upper = zero

    !====================================================================================
    !   Find a diagonal matrix E (stored as vector) such that K = (L+E)E^-1(U+E) ~= lhs
    !===================================================================================

    ! we only calculate E_inv because we mostly divide by E terms
    if ( .true. ) then

      ! Better, but k-dependence
      do i = 1, ngrdcol

        E_inv(i,1) = 1.0_core_rknd / lhs(0,i,1)
        E_inv(i,2) = 1.0_core_rknd & 
                     / ( lhs(0,i,2) - ( lhs(1,i,2) * lhs(-1,i,1) ) * E_inv(i,1) &
                                    - ( lhs(1,i,2) * lhs(-2,i,1) ) * E_inv(i,1) )

        do k = 3, ndim
          E_inv(i,k) = 1.0_core_rknd & 
                       / ( lhs(0,i,k) - ( lhs(2,i,k) * lhs(-1,i,k-2) ) * E_inv(i,k-2) &
                                      - ( lhs(2,i,k) * lhs(-2,i,k-2) ) * E_inv(i,k-2) &
                                      - ( lhs(1,i,k) * lhs(-1,i,k-1) ) * E_inv(i,k-1) &
                                      - ( lhs(1,i,k) * lhs(-2,i,k-1) ) * E_inv(i,k-1) )
        end do

      end do

    else

      ! No k-dependence, but results in slower convergence because it's a worse estimate
      ! still interesting, and should be tested more
      do i = 1, ngrdcol
        E_inv(i,1) = 1.0_core_rknd / lhs(0,i,1)
        E_inv(i,2) = 1.0_core_rknd & 
                     / ( lhs(0,i,2) - ( lhs(1,i,2) * lhs(-1,i,1) ) / lhs(0,i,1) &
                                    - ( lhs(1,i,2) * lhs(-2,i,1) ) / lhs(0,i,1) )
        do k = 3, ndim
          E_inv(i,k) = 1.0_core_rknd & 
                       / ( lhs(0,i,k) - ( lhs(2,i,k) * lhs(-1,i,k-2) ) / lhs(0,i,k-2) &
                                      - ( lhs(2,i,k) * lhs(-2,i,k-2) ) / lhs(0,i,k-2) &
                                      - ( lhs(1,i,k) * lhs(-1,i,k-1) ) / lhs(0,i,k-1) &
                                      - ( lhs(1,i,k) * lhs(-2,i,k-1) ) / lhs(0,i,k-1) )
        end do
      end do

    end if

    ! It's interesting to see how close E_inv is to the lhs diagonal, it's usually (visually) very close, but 
    ! using it instead of E_inv requires so many more iterations (use it would be jacobi preconditioner)
    ! do k = 1, ndim
    !   print *, "E -(",k,")- D = ", 1.0_core_rknd / E_inv(1,k), " -- ", lhs(0,1,k)
    ! end do


    !======================================================================
    !   Calculate lower bands of (L+E)^-1 and upper bands of (U+E)^-1
    !======================================================================
    
    ! First lower band of (L+E)^-1
    do k = 2, ndim
      do i = 1, ngrdcol
          lower(i,1,k) = - E_inv(i,k) * E_inv(i,k-1) * lhs(1,i,k)
      end do
    end do 

    ! First upper band of (U+E)^-1
    do k = 1, ndim-1
      do i = 1, ngrdcol
          upper(i,1,k) = - E_inv(i,k) * E_inv(i,k+1) * lhs(-1,i,k)
      end do
    end do  

    ! Second lower band of (L+E)^-1
    do k = 3, ndim
      do i = 1, ngrdcol
        lower(i,2,k) = - E_inv(i,k) * ( lhs(2,i,k) * E_inv(i,k-2) &
                                      + lhs(1,i,k) * lower(i,1,k-1) )
      end do
    end do

    ! Second upper band of (U+E)^-1
    do k = 1, ndim-2
      do i = 1, ngrdcol
        upper(i,2,k) = - E_inv(i,k) * ( lhs(-2,i,k) * E_inv(i,k+2) &
                                      + lhs(-1,i,k) * upper(i,1,k+1) )
      end do
    end do

    ! Since lhs is 5 banded (2 lower, 2 upper), we can generalize this step
    ! for bands > 2
    do b = 3, n_pc_bands

      ! bth lower band of (L+E)^-1
      do k = 3, ndim
        do i = 1, ngrdcol
          lower(i,b,k) = - E_inv(i,k) * ( lhs(2,i,k) * lower(i,b-2,k-2) &
                                        + lhs(1,i,k) * lower(i,b-1,k-1) )
        end do
      end do

      ! bth upper band of (U+E)^-1
      do k = 1, ndim-2
        do i = 1, ngrdcol
          upper(i,b,k) = - E_inv(i,k) * ( lhs(-2,i,k) * upper(i,b-2,k+2) &
                                        + lhs(-1,i,k) * upper(i,b-1,k+1) )
        end do
      end do

    end do 

    !======================================================================
    !     Compute y_vec = A^-1 x_vec ~= K^-1 x_vec = (U+E)^-1 E (L+E)^-1 x_vec
    !======================================================================

    ! calculate y_tmp = (L+E)^-1 x_vec
    do k = 1, ndim
      do i = 1, ngrdcol

        y_tmp(i,k) = x_vec(i,k) * E_inv(i,k)

        do b = 1, min( k-1, n_pc_bands-1 )
          y_tmp(i,k) = y_tmp(i,k) + x_vec(i,k-b) * lower(i,b,k)
        end do

      end do
    end do

    ! calculate y_tmp = E y_tmp
    do k = 1, ndim
      do i = 1, ngrdcol
        y_tmp(i,k) = y_tmp(i,k) / E_inv(i,k)
      end do
    end do

    ! calculate y_vec = (U+E)^-1 y_tmp  
    ! this step is skipped in the offcial version of the algorithm, which uses K = K1 K2 and
    ! omega = <K1^-1 t, K1^-1 s> / <K1^-1 t, K1^-1 t>, which might be worth investigating, 
    ! we might want to make the above step (y_tmp = E y_tmp) into (y_tmp = E^(1/2) y_tmp) instead
    do k = 1, ndim
      do i = 1, ngrdcol

        y_vec(i,k) = y_tmp(i,k) * E_inv(i,k)

        do b = 1, min( n_pc_bands, ndim-k)
          y_vec(i,k) = y_vec(i,k) + y_tmp(i,k+b) * upper(i,b,k)
        end do

      end do
    end do

  end subroutine precond_psuedo_solve

end module bicgstab_solvers

