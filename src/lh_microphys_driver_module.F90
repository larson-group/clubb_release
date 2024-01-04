!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module lh_microphys_driver_module

  implicit none

  private ! Default scope

#ifdef SILHS
  public :: lh_microphys_driver

contains

  !=============================================================================
  subroutine lh_microphys_driver( &
               gr, dt, nz, num_samples, &
               pdf_dim, hydromet_dim, hm_metadata, &
               X_nl_all_levs, lh_sample_point_weights, &
               pdf_params, precip_fracs, p_in_Pa, exner, rho, &
               rcm, delta_zt, cloud_frac, &
               hydromet, X_mixt_comp_all_levs,  &
               lh_rt_clipped, lh_thl_clipped, & ! In
               lh_rc_clipped, lh_rv_clipped, & ! In
               lh_Nc_clipped, & ! In
               l_lh_importance_sampling, &
               l_lh_instant_var_covar_src, &
               saturation_formula, &
               stats_metadata, &
               stats_zt, stats_zm, stats_sfc, stats_lh_zt, &
               lh_hydromet_mc, lh_hydromet_vel, lh_Ncm_mc, &
               lh_rcm_mc, lh_rvm_mc, lh_thlm_mc, &
               lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, &
               lh_wpthlp_mc, lh_rtpthlp_mc, &
               lh_AKm, AKm, AKstd, AKstd_cld, &
               lh_rcm_avg, AKm_rcm, AKm_rcc, &
               microphys_sub )

    ! Description:
    !   Computes an estimate of the change due to microphysics given a set of
    !   subcolumns of thlm, rtm, et cetera from the subcolumn generator.
    !
    ! References:
    !   None
    !---------------------------------------------------------------------------

    use grid_class, only: &
      grid ! Type

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use hydromet_pdf_parameter_module, only: &
      precipitation_fractions ! Type

    use est_kessler_microphys_module, only: &
      est_kessler_microphys

    use clubb_precision, only: &
      core_rknd
      
    use error_code, only: &
       clubb_at_least_debug_level  ! Procedure

    use estimate_scm_microphys_module, only: &
      est_single_column_tndcy

    use stats_type, only: &
      stats ! Type

    use stats_variables, only: &
      stats_metadata_type

    use corr_varnce_module, only: &
      hm_metadata_type

    implicit none

    ! Interface block
#include "microphys_interface.inc"

    ! Input Variables
    type(grid), target, intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in) :: &
      dt ! Model timestep       [s]

    integer, intent(in) :: &
      num_samples,  & ! Number of calls to microphysics per timestep (normally=2)
      nz,           & ! Number of vertical model levels
      pdf_dim,      & ! Number of variables to sample
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(num_samples,nz,pdf_dim) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(in), dimension(num_samples,nz) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(in), dimension(num_samples,nz) :: &
      lh_sample_point_weights ! Weight given the individual sample points

    type(pdf_parameter), intent(in) :: & 
      pdf_params ! PDF parameters       [units vary]

    type(precipitation_fractions), intent(in) :: &
      precip_fracs           ! Precipitation fractions      [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac,  & ! Cloud fraction               [-]
      delta_zt,    & ! Change in meters with height [m]
      rcm,         & ! Liquid water mixing ratio    [kg/kg]
      p_in_Pa,     & ! Pressure                     [Pa]
      exner,       & ! Exner function               [-]
      rho            ! Density on thermo. grid      [kg/m^3]
      
    real( kind = core_rknd ), dimension(num_samples,nz), intent(in) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points

    logical, intent(in) :: &
      l_lh_importance_sampling, & ! Do importance sampling (SILHS) [-]
      l_lh_instant_var_covar_src  ! Produce instantaneous var/covar tendencies [-]

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    type(stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      lh_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_Ncm_mc,     & ! LH estimate of time tndcy. of cloud droplet conc.     [num/kg/s]
      lh_rcm_mc,     & ! LH estimate of time tndcy. of liq. water mixing ratio [kg/kg/s]
      lh_rvm_mc,     & ! LH estimate of time tndcy. of vapor water mix. ratio  [kg/kg/s]
      lh_thlm_mc,    & ! LH estimate of time tndcy. of liquid potential temp.  [K/s]
      lh_rtp2_mc,    & ! LH microphysics tendency for <rt'^2>                  [(kg/kg)^2/s]
      lh_thlp2_mc,   & ! LH microphysics tendency for <thl'^2>                 [K^2/s]
      lh_wprtp_mc,   & ! LH microphysics tendency for <w'rt'>                  [m*(kg/kg)/s^2]
      lh_wpthlp_mc,  & ! LH microphysics tendency for <w'thl'>                 [m*K/s^2]
      lh_rtpthlp_mc    ! LH microphysics tendency for <rt'thl'>                [K*(kg/kg)/s]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_AKm,     & ! Kessler ac estimate                 [kg/kg/s]
      AKm,        & ! Exact Kessler ac                    [kg/kg/s]
      AKstd,      & ! St dev of exact Kessler ac          [kg/kg/s]
      AKstd_cld,  & ! Stdev of exact w/in cloud ac        [kg/kg/s]
      lh_rcm_avg, & ! Monte Carlo rcm estimate            [kg/kg]
      AKm_rcm,    & ! Kessler ac based on rcm             [kg/kg/s]
      AKm_rcc       ! Kessler ac based on rcm/cloud_frac  [kg/kg/s]

    ! ---- Begin Code ----

    ! Perform LH and analytic microphysical calculations
    ! As a test of SILHS, compute an estimate of Kessler microphysics
    if ( clubb_at_least_debug_level( 2 ) ) then
       call est_kessler_microphys &
            ( nz, num_samples, pdf_dim, &                        ! Intent(in)
              X_nl_all_levs, pdf_params, rcm, cloud_frac, &      ! Intent(in)
              X_mixt_comp_all_levs, lh_sample_point_weights, &   ! Intent(in)
              l_lh_importance_sampling, &                        ! Intent(in)
              lh_AKm, AKm, AKstd, AKstd_cld, &                   ! Intent(out)
              AKm_rcm, AKm_rcc, lh_rcm_avg )                     ! Intent(out)
    end if

    ! Call the latin hypercube microphysics driver for microphys_sub
    call est_single_column_tndcy( &
           gr, dt, nz, num_samples, &                                  ! Intent(in)
           pdf_dim, hydromet_dim, hm_metadata, &                    ! Intent(in)
           X_nl_all_levs, X_mixt_comp_all_levs, &                      ! Intent(in)
           lh_sample_point_weights, pdf_params, precip_fracs, &        ! Intent(in)
           p_in_Pa, exner, rho, &                                      ! Intent(in)
           delta_zt, hydromet, rcm, &                                  ! Intent(in)
           lh_rt_clipped, lh_thl_clipped, &                            ! Intent(in)
           lh_rc_clipped, lh_rv_clipped, &                             ! Intent(in)
           lh_Nc_clipped, &                                            ! Intent(in)
           l_lh_instant_var_covar_src, &                               ! Intent(in)
           saturation_formula, &                                       ! Intent(in)
           stats_metadata, &                                           ! Intent(in) 
           stats_zt, stats_zm, stats_sfc, stats_lh_zt, &               ! intent(inout)
           lh_hydromet_mc, lh_hydromet_vel, lh_Ncm_mc, &               ! Intent(out)
           lh_rvm_mc, lh_rcm_mc, lh_thlm_mc, &                         ! Intent(out)
           lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, &                     ! Intent(out)
           lh_wpthlp_mc, lh_rtpthlp_mc, &                              ! Intent(out)
           microphys_sub )                                             ! Intent(Procedure)

    return
  end subroutine lh_microphys_driver
#endif /* SILHS */

end module lh_microphys_driver_module
