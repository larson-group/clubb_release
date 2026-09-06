!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module lh_microphys_var_covar_module

  implicit none

  public :: lh_microphys_var_covar_driver_api

  private ! Default scope

  contains

  !-----------------------------------------------------------------------
  subroutine lh_microphys_var_covar_driver_api( &
               nzt, num_samples, ngrdcol, dt, lh_sample_point_weights, & ! In
               pdf_params, lh_rt_all, lh_thl_all, lh_w_all, & ! In
               lh_rcm_mc_all, lh_rvm_mc_all, lh_thlm_mc_all, & ! In
               l_lh_instant_var_covar_src, & ! In
               lh_rtp2_mc_zt, lh_thlp2_mc_zt, lh_wprtp_mc_zt, & ! Out
               lh_wpthlp_mc_zt, lh_rtpthlp_mc_zt ) ! Out

  ! Description:
  !   Computes the effect of microphysics on gridbox variances and covariances
  !   for all model columns, including single-column runs with ngrdcol=1.
  !
  ! More description:
  !   The equations for the (co)variance microphysical tendencies, when
  !   integrated forward in time explicitly, are:
  !
  !   rtp2_mc    = 2*covar(rt,rt_mc) + dt*var(rt_mc)
  !   thlp2_mc   = 2*covar(thl,thl_mc) + dt*var(thl_mc)
  !   wprtp_mc   = covar(w,rt_mc)
  !   wpthlp_mc  = covar(w,thl_mc)
  !   rtpthlp_mc = covar(thl,rt_mc) + covar(rt,thl_mc) + dt*covar(rt_mc,thl_mc)
  !
  !   This code can optionally take the limit of these equations at an
  !   infinitesimally small time step, such that the terms involving
  !   dt drop out. This configuration agrees with the KK upscaled analytic
  !   solution. (See clubb:ticket:753 for more discussion on this.)
  !
  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Constant

    use constants_clubb, only: &
      one, &
      two

    use math_utilities, only: &
      compute_sample_mean, &
      compute_sample_variance,    &
      compute_sample_covariance

    use pdf_parameter_module, only: &
      pdf_parameter

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nzt,           &                  ! Number of vertical levels
      num_samples, &       ! Number of SILHS sample points
      ngrdcol              ! Number of model columns

    real( kind = core_rknd ), intent(in) :: &
      dt                               ! Model time step                             [s]

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      lh_sample_point_weights, &  ! Weight of SILHS sample points
      lh_rt_all, &                     ! SILHS samples of total water                [kg/kg]
      lh_thl_all, &                    ! SILHS samples of potential temperature      [K]
      lh_w_all, &                      ! SILHS samples of vertical velocity          [m/s]
      lh_rcm_mc_all, &                 ! SILHS microphys. tendency of rcm            [kg/kg/s]
      lh_rvm_mc_all, &                 ! SILHS microphys. tendency of rvm            [kg/kg/s]
      lh_thlm_mc_all                   ! SILHS microphys. tendency of thlm           [K/s]

    type(pdf_parameter), intent(in) :: &
      pdf_params           ! The PDF parameters_silhs

    logical, intent(in) :: &
      l_lh_instant_var_covar_src       ! Produce instantaneous var/covar tendencies  [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      lh_rtp2_mc_zt,   &               ! SILHS microphys. est. tendency of <rt'^2>   [(kg/kg)^2/s]
      lh_thlp2_mc_zt,  &               ! SILHS microphys. est. tendency of <thl'^2>  [K^2/s]
      lh_wprtp_mc_zt,  &               ! SILHS microphys. est. tendency of <w'rt'>   [m*(kg/kg)/s^2]
      lh_wpthlp_mc_zt, &               ! SILHS microphys. est. tendency of <w'thl'>  [m*K/s^2]
      lh_rtpthlp_mc_zt                 ! SILHS microphys. est. tendency of <rt'thl'> [K*(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt) :: &
      lh_rt_mc_all

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      mean_rt,         &
      mean_rt_mc,      &
      covar_rt_rt_mc,   &
      mean_thl,        &
      mean_thl_mc,     &
      covar_thl_thl_mc, &
      mean_w,          &
      covar_w_rt_mc,    &
      covar_w_thl_mc,   &
      covar_thl_rt_mc,  &
      covar_rt_thl_mc, &
      var_rt_mc, &
      var_thl_mc, &
      covar_rt_mc_thl_mc

    !----- Begin Code -----

    lh_rt_mc_all = lh_rcm_mc_all + lh_rvm_mc_all

    ! Calculate means, variances, and covariances needed for the tendency terms
    mean_rt = pdf_params%mixt_frac * pdf_params%rt_1 &
      + ( one - pdf_params%mixt_frac ) * pdf_params%rt_2
    mean_thl = pdf_params%mixt_frac * pdf_params%thl_1 &
      + ( one - pdf_params%mixt_frac ) * pdf_params%thl_2
    mean_w = pdf_params%mixt_frac * pdf_params%w_1 &
      + ( one - pdf_params%mixt_frac ) * pdf_params%w_2

    mean_rt_mc = compute_sample_mean( nzt, num_samples, ngrdcol, &
                                       lh_sample_point_weights, lh_rt_mc_all )
    covar_rt_rt_mc = compute_sample_covariance( nzt, num_samples, ngrdcol, &
                                                 lh_sample_point_weights, lh_rt_all, mean_rt, &
                                                 lh_rt_mc_all, mean_rt_mc )
    mean_thl_mc = compute_sample_mean( nzt, num_samples, ngrdcol, &
                                        lh_sample_point_weights, lh_thlm_mc_all )
    covar_thl_thl_mc = compute_sample_covariance( nzt, num_samples, ngrdcol, &
                                                   lh_sample_point_weights, lh_thl_all, mean_thl, &
                                                   lh_thlm_mc_all, mean_thl_mc )
    covar_w_rt_mc = compute_sample_covariance( nzt, num_samples, ngrdcol, &
                                                lh_sample_point_weights, lh_w_all, mean_w, &
                                                lh_rt_mc_all, mean_rt_mc )
    covar_w_thl_mc = compute_sample_covariance( nzt, num_samples, ngrdcol, &
                                                 lh_sample_point_weights, lh_w_all, mean_w, &
                                                 lh_thlm_mc_all, mean_thl_mc )
    covar_thl_rt_mc = compute_sample_covariance( nzt, num_samples, ngrdcol, &
                                                  lh_sample_point_weights, lh_thl_all, mean_thl, &
                                                  lh_rt_mc_all, mean_rt_mc )
    covar_rt_thl_mc = compute_sample_covariance( nzt, num_samples, ngrdcol, &
                                                  lh_sample_point_weights, lh_rt_all, mean_rt, &
                                                  lh_thlm_mc_all, mean_thl_mc )

    ! Compute the microphysical variance and covariance tendencies
    lh_rtp2_mc_zt = two * covar_rt_rt_mc
    lh_thlp2_mc_zt = two * covar_thl_thl_mc
    lh_wprtp_mc_zt   = covar_w_rt_mc
    lh_wpthlp_mc_zt  = covar_w_thl_mc
    lh_rtpthlp_mc_zt = covar_thl_rt_mc + covar_rt_thl_mc

    if ( .not. l_lh_instant_var_covar_src ) then
      ! Variances and covariances for timestep-dependent terms
      ! NOTE: these terms arise in rtp2 and thlp2 when rtm and thlm are
      ! explicitly integrated forward in time. These terms are not included
      ! in KK upscaled, so using these terms causes non-convergence with KK
      ! upscaled.
      var_rt_mc = compute_sample_variance( nzt, num_samples, ngrdcol, &
                                            lh_rt_mc_all, lh_sample_point_weights, mean_rt_mc )
      var_thl_mc = compute_sample_variance( nzt, num_samples, ngrdcol, &
                                             lh_thlm_mc_all, lh_sample_point_weights, mean_thl_mc )
      covar_rt_mc_thl_mc = compute_sample_covariance( nzt, num_samples, ngrdcol, &
                                                        lh_sample_point_weights, lh_rt_mc_all, &
                                                        mean_rt_mc, &
                                                        lh_thlm_mc_all, mean_thl_mc )

      ! Add timestep-dependent terms
      lh_rtp2_mc_zt = lh_rtp2_mc_zt + dt * var_rt_mc
      lh_thlp2_mc_zt = lh_thlp2_mc_zt + dt * var_thl_mc
      lh_rtpthlp_mc_zt = lh_rtpthlp_mc_zt + dt * covar_rt_mc_thl_mc
    endif

  end subroutine lh_microphys_var_covar_driver_api
  !-----------------------------------------------------------------------
end module lh_microphys_var_covar_module
