!-------------------------------------------------------------------------
! $Id$
!===============================================================================
module sigma_sqd_w_module

  implicit none

  public :: compute_sigma_sqd_w

  private ! Default scope

  contains

  subroutine compute_sigma_sqd_w( nzm, nzt, ngrdcol, gr, &
                                  wp3, wp2, thlp2, rtp2, &
                                  up2, vp2, wpthlp, wprtp, upwp, vpwp, &
                                  clubb_params, &
                                  l_predict_upwp_vpwp, &
                                  sigma_sqd_w )

    ! Description:
    ! Compute the variable sigma_sqd_w (PDF width parameter).
    !
    ! The value of sigma_sqd_w is restricted in the ADG1 PDF in order to keep
    ! the marginal PDFs of all responder variables (variables that do not set
    ! the mixture fraction) valid.  The limits on sigma_sqd_w in order to keep
    ! the PDF of a responder variable, x, valid are:
    !
    ! 0 <= sigma_sqd_w <= 1 - <w'x'>^2 / ( <w'^2> * <x'^2> ).
    !
    ! The overall limits on sigma_sqd_w must be applied based on the most
    ! restrictive case so that all Double Gaussian PDF responder variables, x,
    ! have realizable PDFs.  The overall limits on sigma_sqd_w are:
    !
    ! 0 <= sigma_sqd_w <= 1 - max( <w'x'>^2 / ( <w'^2> * <x'^2> ), for all x).
    !
    ! The equation used for sigma_sqd_w is:
    !
    ! sigma_sqd_w = gamma_Skw_fnc
    !               * ( 1 - max( <w'x'>^2 / ( <w'^2> * <x'^2> ), for all x) );
    !
    ! where 0 <= gamma_Skw_fnc <= 1.

    ! References:
    !   Eqn 22 in ``Equations for CLUBB''
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,         & ! Constant(s)
        w_tol,       &
        rt_tol,      &
        thl_tol,     &
        one_hundred, &
        w_tol_sqd,   &
        zero_threshold

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use grid_class, only: &
        grid,       & ! Type
        zt2zm_api,  & ! Procedure(s)
        zm2zt2zm      ! Procedure(s)

    use parameter_indices, only: &
        nparams

    use Skx_module, only: &
        Skx_func, &
        compute_gamma_Skw

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol

    type(grid), intent(in) :: &
      gr

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      wp3             ! Third moment of vertical velocity            [m^3/s^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      wp2,           & ! Variance of vertical velocity               [m^2/s^2]
      thlp2,         & ! Variance of liquid water potential temp.    [K^2]
      rtp2,          & ! Variance of total water mixing ratio        [kg^2/kg^2]
      up2,           & ! Variance of west-east horizontal velocity   [m^2/s^2]
      vp2,           & ! Variance of south-north horizontal velocity [m^2/s^2]
      wpthlp,        & ! Flux of liquid water potential temp.        [m/s K]
      wprtp,         & ! Flux of total water mixing ratio            [m/s kg/kg]
      upwp,          & ! Flux of west-east horizontal velocity       [m^2/s^2]
      vpwp             ! Flux of south-north horizontal velocity     [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params     ! Array of CLUBB tunable parameters           [units vary]

    logical, intent(in) :: &
      l_predict_upwp_vpwp ! Flag to predict <u'w'> and <v'w'> along with <u> and <v> alongside the
                          ! advancement of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>, and <w'sclr'> in
                          ! subroutine advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                          ! approximated by eddy diffusivity when <u> and <v> are advanced in
                          ! subroutine advance_windm_edsclrm.

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
      sigma_sqd_w ! PDF width parameter      [-]

    ! Local Variable
    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      wp3_zm,         & ! Third moment of w on momentum levels            [m^3/s^3]
      Skw_zm,         & ! Skewness of w on momentum levels                [-]
      gamma_Skw_fnc,  & ! Gamma as a function of skewness                 [-]
      sigma_sqd_w_tmp   ! Unsmoothed PDF width parameter                  [-]

    real( kind = core_rknd ) :: &
      max_corr_w_x_sqd   ! Max. val. of wpxp^2/(wp2*xp2) for all vars. x  [-]

    integer :: i, k

    ! ---- Begin Code ----

    !$acc enter data create( wp3_zm, Skw_zm, gamma_Skw_fnc, sigma_sqd_w_tmp )

    if ( nzm > 1 ) then
      wp3_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, wp3(:,:) )
    else
      ! Skip interpolation if nzm == 1, this only occurs for testing purposes
      wp3_zm(:,:) = wp3(:,:)
    end if

    call Skx_func( nzm, ngrdcol, wp2, wp3_zm, &
                   w_tol, clubb_params, &
                   Skw_zm )

    call compute_gamma_Skw( nzm, ngrdcol, Skw_zm, clubb_params, & ! In
                            gamma_Skw_fnc )                       ! Out

    !----------------------------------------------------------------
    ! Compute sigma_sqd_w with new formula from Vince
    !----------------------------------------------------------------

    ! Find the maximum value of <w'x'>^2 / ( <w'^2> * <x'^2> ) for all
    ! variables x that are Double Gaussian PDF responder variables.  This
    ! includes rt and theta-l.  When l_predict_upwp_vpwp is enabled, u and v are
    ! also calculated as part of the PDF, and they are included as well.
    ! Additionally, when sclr_dim > 0, passive scalars (sclr) are also included.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol

        max_corr_w_x_sqd = max( ( wpthlp(i,k) / ( sqrt( wp2(i,k) * thlp2(i,k) ) &
                                  + one_hundred * w_tol * thl_tol ) )**2, &
                                  ( wprtp(i,k) / ( sqrt( wp2(i,k) * rtp2(i,k) )  &
                                  + one_hundred * w_tol * rt_tol ) )**2 )

        if ( l_predict_upwp_vpwp ) then
          max_corr_w_x_sqd = max( max_corr_w_x_sqd, &
                                  ( upwp(i,k) / ( sqrt( up2(i,k) * wp2(i,k) ) &
                                  + one_hundred * w_tol_sqd ) )**2, &
                                  ( vpwp(i,k) / ( sqrt( vp2(i,k) * wp2(i,k) ) &
                                  + one_hundred * w_tol_sqd ) )**2 )
        endif ! l_predict_upwp_vpwp

        sigma_sqd_w_tmp(i,k) = gamma_Skw_fnc(i,k) * ( one - min( max_corr_w_x_sqd, one ) )

      end do
    end do
    !$acc end parallel loop

    if ( nzm > 1 ) then
      ! Smooth in the vertical using interpolation.
      sigma_sqd_w(:,:) = zm2zt2zm( nzm, nzt, ngrdcol, gr, sigma_sqd_w_tmp(:,:), &
                                  zero_threshold )
    else
      ! Skip interpolation if nzm == 1, this only occurs for testing purposes
      sigma_sqd_w(:,:) = sigma_sqd_w_tmp(:,:)
    end if

    !$acc exit data delete( wp3_zm, Skw_zm, gamma_Skw_fnc, sigma_sqd_w_tmp )

    return

  end subroutine compute_sigma_sqd_w

end module sigma_sqd_w_module
