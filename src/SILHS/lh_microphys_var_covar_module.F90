!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module lh_microphys_var_covar_module

  implicit none

  public :: lh_microphys_var_covar_driver

  private ! Default scope

  contains

  !-----------------------------------------------------------------------
  subroutine lh_microphys_var_covar_driver &
             ( nz, num_samples, dt, lh_sample_point_weights, &
               lh_rt_all, lh_thl_all, lh_w_all, &
               lh_rcm_mc_all, lh_rvm_mc_all, lh_thlm_mc_all, &
               lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, &
               lh_wpthlp_mc, lh_rtpthlp_mc )

  ! Description:
  !   Computes the effect of microphysics on gridbox variances and covariances

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd ! Constant

    use stats_variables, only: &
      l_stats_samp, & ! Variable(s)
      ilh_rtp2_mc, ilh_thlp2_mc, ilh_wprtp_mc, ilh_wpthlp_mc, ilh_rtpthlp_mc, &
      stats_zm

    use stats_type_utilities, only: &
      stat_update_var ! Procedure

    use math_utilities, only: &
      compute_sample_mean,        & ! Procedure(s)
      compute_sample_variance,    &
      compute_sample_covariance

    use grid_class, only: &
      zt2zm

    use constants_clubb, only: &
      two    ! Constant

    implicit none

    ! Input Variables!
    integer, intent(in) :: &
      nz,           &                  ! Number of vertical levels
      num_samples                      ! Number of SILHS sample points

    real( kind = core_rknd ), intent(in) :: &
      dt                               ! Model time step                             [s]

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights          ! Weight of SILHS sample points

    real( kind = core_rknd ), dimension(nz,num_samples), intent(in) :: &
      lh_rt_all, &                     ! SILHS samples of total water                [kg/kg]
      lh_thl_all, &                    ! SILHS samples of potential temperature      [K]
      lh_w_all, &                      ! SILHS samples of vertical velocity          [m/s]
      lh_rcm_mc_all, &                 ! SILHS microphys. tendency of rcm            [kg/kg/s]
      lh_rvm_mc_all, &                 ! SILHS microphys. tendency of rvm            [kg/kg/s]
      lh_thlm_mc_all                   ! SILHS microphys. tendency of thlm           [K/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_rtp2_mc,     &                ! SILHS microphys. est. tendency of <rt'^2>   [(kg/kg)^2/s]
      lh_thlp2_mc,    &                ! SILHS microphys. est. tendency of <thl'^2>  [K^2/s]
      lh_wprtp_mc,    &                ! SILHS microphys. est. tendency of <w'rt'>   [m*(kg/kg)/s^2]
      lh_wpthlp_mc,   &                ! SILHS microphys. est. tendency of <w'thl'>  [m*K/s^2]
      lh_rtpthlp_mc                    ! SILHS microphys. est. tendency of <rt'thl'> [K*(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      lh_rt_mc_all

    real( kind = core_rknd ), dimension(nz) :: &
      mean_rt,         &
      mean_rt_mc,      &
      covar_rt_rtmc,   &
      mean_thl,        &
      mean_thl_mc,     &
      covar_thl_thlmc, &
      mean_w,          &
      covar_w_rtmc,    &
      covar_w_thlmc,   &
      covar_thl_rtmc,  &
      covar_rt_thlmc

    ! For timestep-dependent terms
    real( kind = core_rknd ), dimension(nz) :: &
      var_rt_mc, &
      var_thl_mc, &
      covar_rtmc_thlmc

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    lh_rt_mc_all = lh_rcm_mc_all + lh_rvm_mc_all

    ! Calculate means, variances, and covariances needed for the tendency terms
    mean_rt = compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_rt_all )
    mean_rt_mc = compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_rt_mc_all )
    covar_rt_rtmc = compute_sample_covariance( nz, num_samples, lh_sample_point_weights, &
                                               lh_rt_all, mean_rt, lh_rt_mc_all, mean_rt_mc )
    mean_thl = compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_thl_all )
    mean_thl_mc = compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_thlm_mc_all )
    covar_thl_thlmc = compute_sample_covariance( nz, num_samples, lh_sample_point_weights, &
                                                 lh_thl_all, mean_thl, lh_thlm_mc_all, mean_thl_mc )
    mean_w = compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_w_all )
    covar_w_rtmc = compute_sample_covariance( nz, num_samples, lh_sample_point_weights, &
                                              lh_w_all, mean_w, lh_rt_mc_all, mean_rt_mc )
    covar_w_thlmc = compute_sample_covariance( nz, num_samples, lh_sample_point_weights, &
                                              lh_w_all, mean_w, lh_thlm_mc_all, mean_thl_mc )
    covar_thl_rtmc = compute_sample_covariance( nz, num_samples, lh_sample_point_weights, &
                                                lh_thl_all, mean_thl, lh_rt_mc_all, mean_rt_mc )
    covar_rt_thlmc = compute_sample_covariance( nz, num_samples, lh_sample_point_weights, &
                                                lh_rt_all, mean_rt, lh_thlm_mc_all, mean_thl_mc )

    ! Variances and covariances for time-dependent terms
    var_rt_mc  = compute_sample_variance( nz, num_samples, lh_rt_mc_all, lh_sample_point_weights, &
                                          mean_rt_mc  )
    var_thl_mc = compute_sample_variance &
                                ( nz, num_samples, lh_thlm_mc_all, lh_sample_point_weights, &
                                  mean_thl_mc )
    covar_rtmc_thlmc = compute_sample_covariance &
                              ( nz, num_samples, lh_sample_point_weights, &
                                lh_rt_mc_all, mean_rt_mc, lh_thlm_mc_all, mean_thl_mc )

    ! Compute the microphysical variance and covariance tendencies
    lh_rtp2_mc    = zt2zm( two*covar_rt_rtmc   + dt*var_rt_mc  )
    lh_thlp2_mc   = zt2zm( two*covar_thl_thlmc + dt*var_thl_mc )
    lh_wprtp_mc   = zt2zm( covar_w_rtmc  )
    lh_wpthlp_mc  = zt2zm( covar_w_thlmc )
    lh_rtpthlp_mc = zt2zm( covar_thl_rtmc + covar_rt_thlmc + dt*covar_rtmc_thlmc )


    ! Stats sampling
    if ( l_stats_samp ) then
      call stat_update_var( ilh_rtp2_mc, lh_rtp2_mc, stats_zm )
      call stat_update_var( ilh_thlp2_mc, lh_thlp2_mc, stats_zm )
      call stat_update_var( ilh_wprtp_mc, lh_wprtp_mc, stats_zm )
      call stat_update_var( ilh_wpthlp_mc, lh_wpthlp_mc, stats_zm )
      call stat_update_var( ilh_rtpthlp_mc, lh_rtpthlp_mc, stats_zm )
    end if

    return
  end subroutine lh_microphys_var_covar_driver
  !-----------------------------------------------------------------------
end module lh_microphys_var_covar_module
