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
  subroutine lh_microphys_driver &
             ( dt, nz, num_samples, d_variables, &
               X_nl_all_levs, lh_sample_point_weights, &
               pdf_params, hydromet_pdf_params, p_in_Pa, exner, rho, &
               rcm, delta_zt, cloud_frac, &
               hydromet, X_mixt_comp_all_levs,  &
               lh_clipped_vars, &
               lh_hydromet_mc, lh_hydromet_vel, lh_Ncm_mc, &
               lh_rcm_mc, lh_rvm_mc, lh_thlm_mc, &
               lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, &
               lh_wpthlp_mc, lh_rtpthlp_mc, &
               microphys_sub )

    ! Description:
    !   Computes an estimate of the change due to microphysics given a set of
    !   subcolumns of thlm, rtm, et cetera from the subcolumn generator.
    !
    ! References:
    !   None
    !---------------------------------------------------------------------------

    use variables_diagnostic_module, only: & 
      lh_AKm,  & 
      AKm, & 
      AKstd, & 
      AKstd_cld, & 
      AKm_rcm, & 
      AKm_rcc, & 
      lh_rcm_avg

    use parameters_model, only: hydromet_dim ! Variable

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter ! Type

    use est_kessler_microphys_module, only: &
      est_kessler_microphys

    use clubb_precision, only: &
      core_rknd
      
    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use estimate_scm_microphys_module, only: &
      est_single_column_tndcy

    use latin_hypercube_driver_module, only: &
      lh_clipped_variables_type

    implicit none

    ! Interface block
#include "microphys_interface.inc"

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      dt ! Model timestep       [s]

    integer, intent(in) :: &
      d_variables,     & ! Number of variables to sample
      num_samples,   & ! Number of calls to microphysics per timestep (normally=2)
      nz               ! Number of vertical model levels

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(in), dimension(nz,num_samples) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(in), dimension(num_samples) :: &
      lh_sample_point_weights ! Weight given the individual sample points

    type(pdf_parameter), dimension(nz), intent(in) :: & 
      pdf_params ! PDF parameters       [units vary]

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac,  & ! Cloud fraction               [-]
      delta_zt,    & ! Change in meters with height [m]
      rcm,         & ! Liquid water mixing ratio    [kg/kg]
      p_in_Pa,     & ! Pressure                     [Pa]
      exner,       & ! Exner function               [-]
      rho            ! Density on thermo. grid      [kg/m^3]

    type(lh_clipped_variables_type), dimension(nz,num_samples), intent(in) :: &
      lh_clipped_vars

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

    ! ---- Begin Code ----

    ! Perform LH and analytic microphysical calculations
    ! As a test of SILHS, compute an estimate of Kessler microphysics
    if ( clubb_at_least_debug_level( 2 ) ) then
       call est_kessler_microphys &
            ( nz, num_samples, d_variables, &                    ! Intent(in)
              X_nl_all_levs, pdf_params, rcm, cloud_frac, &      ! Intent(in)
              X_mixt_comp_all_levs, lh_sample_point_weights, &   ! Intent(in)
              lh_AKm, AKm, AKstd, AKstd_cld, &                   ! Intent(out)
              AKm_rcm, AKm_rcc, lh_rcm_avg )                     ! Intent(out)
    end if

    ! Call the latin hypercube microphysics driver for microphys_sub
    call est_single_column_tndcy &
         ( dt, nz, num_samples, d_variables, &                         ! Intent(in)
           X_nl_all_levs, X_mixt_comp_all_levs, &                      ! Intent(in)
           lh_sample_point_weights, pdf_params, hydromet_pdf_params, & ! Intent(in)
           p_in_Pa, exner, rho, &                                      ! Intent(in)
           delta_zt, hydromet, rcm, &                                  ! Intent(in)
           lh_clipped_vars, &                                          ! Intent(in)
           lh_hydromet_mc, lh_hydromet_vel, lh_Ncm_mc, &               ! Intent(out)
           lh_rvm_mc, lh_rcm_mc, lh_thlm_mc, &                         ! Intent(out)
           lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, &                     ! Intent(out)
           lh_wpthlp_mc, lh_rtpthlp_mc, &                              ! Intent(out)
           microphys_sub )                                             ! Intent(Procedure)

    return
  end subroutine lh_microphys_driver
#endif /* SILHS */

end module lh_microphys_driver_module
