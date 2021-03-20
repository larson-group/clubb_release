!----------------------------------------------------------------------
! $Id$
module cloud_feedback

  ! Description:
  !   Contains subroutines for the Cloud Feedback cases.
  !
  ! References:
  !   http://atmgcm.msrc.sunysb.edu/cfmip_figs/Case_specification.html
  !----------------------------------------------------------------------

  implicit none

  private ! Default Scope

  public :: cloud_feedback_sfclyr

  contains

  !----------------------------------------------------------------------
  subroutine cloud_feedback_sfclyr( time, sfctype,                &
                                    thlm_sfc, rtm_sfc, lowest_level, &
                                    ubar, p_sfc, T_sfc,                &
                                    wpthlp_sfc, wprtp_sfc, ustar )

  ! Description:
  !   Sets up surface information for the cloud feedback case

  ! References: 
  !   http://cfmip.metoffice.com/
  !   http://atmgcm.msrc.sunysb.edu/cfmip_figs/Case_specification.html
  !----------------------------------------------------------------------

  use constants_clubb, only: p0, kappa ! Variable(s)

  use saturation, only: sat_mixrat_liq ! Variable(s)

  use clubb_precision, only: time_precision, core_rknd ! Variable(s)

  use sfc_flux, only: compute_wprtp_sfc, compute_wpthlp_sfc

  use time_dependent_input, only: time_sfc_given, T_sfc_given, &! Variable(s)
                                  time_select                   ! Procedure(s)

  use interpolation, only: linear_interp_factor ! Procedure(s)

  implicit none

  intrinsic :: max, sqrt

  integer, intent(in) :: sfctype

  real(time_precision), intent(in) :: &
    time ! The current time [s]

  real( kind = core_rknd ), intent(in) ::  & 
    thlm_sfc,  & ! thlm at (2)         [m/s]
    rtm_sfc,   & ! rtm at (2)          [kg/kg]
    p_sfc,      & ! Surface pressure    [Pa]
    ubar,      & ! This is root (u^2 + v^2), per ATEX and RICO spec.
    lowest_level ! This is z at the lowest above-ground model level.  [m]

  ! Output variables
  real( kind = core_rknd ), intent(out) ::  & 
    T_sfc,      & ! Temperature         [K]
    wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
    wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
    ustar           ! surface friction velocity [m/s]

  ! Constants
  real( kind = core_rknd ), parameter :: & 
  !  rho_sfc_flux = 1._core_rknd0, &
  !  C_10    = 0._core_rknd0013,    & ! Drag coefficient, defined by ATEX specification
    C_h_20  = 0.001094_core_rknd,  & ! Drag coefficient, defined by RICO 3D specification
    C_q_20  = 0.001133_core_rknd,  & ! Drag coefficient, defined by RICO 3D specification
    z0      = 0.00015_core_rknd      ! Roughness length, defined by ATEX specification

  ! Internal variables
  real( kind = core_rknd ) :: &
    Ch,   &                ! This is C_h_20 scaled to the height of the lowest model level.
    Cq,   &                ! This is C_q_20 scaled to the height of the lowest model level.
    exner_sfc, & ! Value of exner at the surface [-]
    time_frac ! The time fraction used for interpolation
   
  integer :: &
    before_time, after_time ! The times used for interpolation

  real( kind = core_rknd ), parameter :: &
    standard_flux_alt = 20._core_rknd ! default height at which the surface flux is computed [m]
 
  !--------------BEGIN CODE---------------------

  call time_select( time, size(time_sfc_given), time_sfc_given, &
                    before_time, after_time, time_frac )

  T_sfc = linear_interp_factor( time_frac, T_sfc_given(after_time), &
                                T_sfc_given(before_time) )

  ! Calculate exner_sfc based on p_sfc.
  exner_sfc = ( p_sfc / p0 )**kappa

  ! Just set ustar = 0.3
  ustar = 0.3_core_rknd

  !--------------------------------------------------------------------------------
  ! Email way
  ! wprtp = value_from_forcings_file_in_W_m**2 / ( rho_sfc_flux * Lv )
  ! wpthlp = value_from_forcings_file_in_W_m**2 / ( rho_sfc_flux * Cp )

  !lhflx(1) = 0.001_core_rknd * ubar * rho_sfc_flux * Lv * ( sat_mixrat_liq( p_sfc, T_sfc ) - & 
  !                                                sat_mixrat_liq( p_sfc, T_in_K ) * 0.8_core_rknd )
  !shflx(1) = 0.001_core_rknd * ubar * rho_sfc_flux * Cp * ( T_sfc - T_in_K )

  ! If this is the S6 case, fudge the values of the fluxes using values from the forcings
  !if ( runtype == "cloud_feedback_s6" .or. runtype == "cloud_feedback_s6_p2k" ) then
  !    wprtp_sfc = lhflx(1) / ( 1.0_core_rknd * Lv )
  !    wpthlp_sfc = shflx(1) / ( 1.0_core_rknd * Cp )
  !else
  !    wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, sat_mixrat_liq( p_sfc, T_sfc ) )
  !    wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, & 
  !                                     T_sfc, exner_sfc )
  !end if

  ! (Stevens, et al. 2000, eq 3)
  ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Ch   = C_h_20 * ((log(standard_flux_alt/z0))/(log(lowest_level/z0))) * & 
         ((log(standard_flux_alt/z0))/(log(lowest_level/z0)))
  ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Cq   = C_q_20 * ((log(standard_flux_alt/z0))/(log(lowest_level/z0))) * & 
         ((log(standard_flux_alt/z0))/(log(lowest_level/z0)))
 
  if ( sfctype == 1 ) then
    wprtp_sfc = compute_wprtp_sfc( Cq, ubar, rtm_sfc, sat_mixrat_liq( p_sfc, T_sfc ) )
    wpthlp_sfc = compute_wpthlp_sfc( Ch, ubar, thlm_sfc, & 
                                     T_sfc, exner_sfc )
  else
    error stop "Invalid value for sfctype."
  end if


  return
  end subroutine cloud_feedback_sfclyr

end module cloud_feedback
