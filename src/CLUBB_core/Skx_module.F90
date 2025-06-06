!-------------------------------------------------------------------------
!$Id$
!===============================================================================
module Skx_module

  implicit none

  private ! Default Scope

  public :: Skx_func, &
            LG_2005_ansatz, &
            xp3_LG_2005_ansatz

  contains

  !-----------------------------------------------------------------------------
  subroutine Skx_func( nz, ngrdcol, xp2, xp3, &
                       x_tol, clubb_params, &
                       Skx )

    ! Description:
    ! Calculate the skewness of x

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd         ! Variable(s)

    use parameter_indices, only: &
      nparams,                 & ! Variable(s)
      iSkw_denom_coef,         &
      iSkw_max_mag

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    ! External
    intrinsic :: min, max

    ! Parameter Constants
    ! Whether to apply clipping to the final result
    logical, parameter ::  &
      l_clipping_kluge = .false.

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      xp2,   & ! <x'^2>               [(x units)^2]
      xp3      ! <x'^3>               [(x units)^3]

    real( kind = core_rknd ), intent(in) :: &
      x_tol     ! x tolerance value                       [(x units)]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      Skx      ! Skewness of x        [-]

    ! Local Variable
    real( kind = core_rknd ), dimension(ngrdcol) :: &
      Skx_denom_tol
      
    integer :: i, k

    ! ---- Begin Code ----

    !$acc data create( Skx_denom_tol )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      Skx_denom_tol(i) = clubb_params(i,iSkw_denom_coef) * x_tol**2
    end do
    !$acc end parallel loop

    !Skx = xp3 / ( max( xp2, x_tol**two ) )**three_halves
    ! Calculation of skewness to help reduce the sensitivity of this value to
    ! small values of xp2.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        Skx(i,k) = xp3(i,k) * sqrt( xp2(i,k) + Skx_denom_tol(i) )**(-3)
      end do
    end do
    !$acc end parallel loop

    ! This is no longer needed since clipping is already
    ! imposed on wp2 and wp3 elsewhere in the code

    ! I turned clipping on in this local copy since thlp3 and rtp3 are not clipped
    if ( l_clipping_kluge ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          Skx(i,k) = min( max( Skx(i,k), -clubb_params(i,iSkw_max_mag) ), &
                          clubb_params(i,iSkw_max_mag) )
        end do
      end do
      !$acc end parallel loop
    end if

    !$acc end data

    return

  end subroutine Skx_func

  !-----------------------------------------------------------------------------
  subroutine LG_2005_ansatz( nz, ngrdcol, Skw, wpxp, wp2, &
                             xp2, beta, sigma_sqd_w, x_tol, &
                             Skx )

    ! Description:
    ! Calculate the skewness of x using the diagnostic ansatz of Larson and
    ! Golaz (2005).

    ! References:
    ! Vincent E. Larson and Jean-Christophe Golaz, 2005:  Using Probability
    ! Density Functions to Derive Consistent Closure Relationships among
    ! Higher-Order Moments.  Mon. Wea. Rev., 133, 1023–1042.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,          & ! Variable(s)
        w_tol_sqd

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    !-------------------------- Input Variables --------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      Skw,         & ! Skewness of w                  [-]
      wpxp,        & ! Turbulent flux of x            [m/s (x units)]
      wp2,         & ! Variance of w                  [m^2/s^2]
      xp2,         & ! Variance of x                  [(x units)^2]
      sigma_sqd_w    ! Normalized variance of w       [-]
      
    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      beta           ! Tunable parameter              [-]
      
    real( kind = core_rknd ), intent(in) :: &
      x_tol          ! Minimum tolerance of x         [(x units)]

    !-------------------------- Output Variable --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      Skx            ! Skewness of x                  [-]

    !-------------------------- Local Variables --------------------------
    real( kind = core_rknd ) :: &
      nrmlzd_corr_wx, & ! Normalized correlation of w and x       [-]
      nrmlzd_Skw        ! Normalized skewness of w                [-]
      
    integer :: i, k

    !--------------------------Begin Code --------------------------

    ! weberjk, 8-July 2015. Commented this out for now. cgils was failing during some tests.

    ! Larson and Golaz (2005) eq. 16
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        nrmlzd_corr_wx = &
                wpxp(i,k) / sqrt( max( wp2(i,k), w_tol_sqd ) &
                             * max( xp2(i,k), x_tol**2 ) * ( one - sigma_sqd_w(i,k) ) )

        ! Larson and Golaz (2005) eq. 11
        nrmlzd_Skw = Skw(i,k) / ( ( one - sigma_sqd_w(i,k)) * sqrt( one - sigma_sqd_w(i,k) ) )

        ! Larson and Golaz (2005) eq. 33
        Skx(i,k) = nrmlzd_Skw * nrmlzd_corr_wx &
              * ( beta(i) + ( one - beta(i) ) * nrmlzd_corr_wx**2 )
      end do
    end do
    !$acc end parallel loop

    return

  end subroutine LG_2005_ansatz

  !-----------------------------------------------------------------------------
  subroutine xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, wpxp_zt, wp2_zt, &
                                 xp2_zt, sigma_sqd_w_zt, &
                                 clubb_params, x_tol, &
                                 xp3 )
    ! Description:
    ! Calculate <x'^3> after calculating the skewness of x using the ansatz of
    ! Larson and Golaz (2005).

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use parameter_indices, only: &
        nparams,                 & ! Variable(s)
        iSkw_denom_coef,         &
        ibeta

    implicit none
    
    !-------------------------- Input Variables--------------------------
    integer, intent(in) :: &
      nzt, &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      Skw_zt,         & ! Skewness of w on thermodynamic levels   [-]
      wpxp_zt,        & ! Flux of x  (interp. to t-levs.)         [m/s(x units)]
      wp2_zt,         & ! Variance of w (interp. to t-levs.)      [m^2/s^2]
      xp2_zt,         & ! Variance of x (interp. to t-levs.)      [(x units)^2]
      sigma_sqd_w_zt    ! Normalized variance of w (interp. to t-levs.)   [-]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    real( kind = core_rknd ), intent(in) :: &
      x_tol             ! Minimum tolerance of x                  [(x units)]

    !-------------------------- Return Variable --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      xp3    ! <x'^3> (thermodynamic levels)    [(x units)^3]

    !-------------------------- Local Variable --------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      Skx_zt       ! Skewness of x on thermodynamic levels    [-]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      Skx_denom_tol

    integer :: i, k

    !-------------------------- Begin Code --------------------------

    !$acc data create( Skx_zt, Skx_denom_tol )
    
    ! Calculate skewness of x using the ansatz of LG05.
    call LG_2005_ansatz( nzt, ngrdcol, Skw_zt, wpxp_zt, wp2_zt, &
                         xp2_zt, clubb_params(:,ibeta), sigma_sqd_w_zt, x_tol, &
                         Skx_zt )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      Skx_denom_tol(i) = clubb_params(i,iSkw_denom_coef) * x_tol**2
    end do
    !$acc end parallel loop

    ! Calculate <x'^3> using the reverse of the special sensitivity reduction
    ! formula in function Skx_func above.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        xp3(i,k) = Skx_zt(i,k) * ( xp2_zt(i,k) + Skx_denom_tol(i) ) &
                               * sqrt( xp2_zt(i,k) + Skx_denom_tol(i) )
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end subroutine xp3_LG_2005_ansatz

  !-----------------------------------------------------------------------------

end module Skx_module
