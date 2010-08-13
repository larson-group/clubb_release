!----------------------------------------------------------------------
! $Id$
module cloud_feedback

  !       Description:
  !       Contains subroutines for the Cloud Feedback cases.
  !----------------------------------------------------------------------

  implicit none

  private ! Default Scope

  public :: cloud_feedback_sfclyr

  contains

  !----------------------------------------------------------------------
  subroutine cloud_feedback_sfclyr( runtype, sfctype,                &
                                    thlm_sfc, rtm_sfc, lowest_level, &
                                    ubar, p_sfc, T_sfc,                &
                                    wpthlp_sfc, wprtp_sfc, ustar )

  !       Description:
  !       Sets up surface information for the cloud feedback case

  !       References:
  !----------------------------------------------------------------------

  use constants_clubb, only: pi, grav, Lv, Cp, p0, kappa ! Variable(s)

  use saturation, only: sat_mixrat_liq ! Variable(s)

  use stats_precision, only: time_precision ! Variable(s)

  use surface_flux, only: compute_wprtp_sfc, compute_wpthlp_sfc

  implicit none

  intrinsic :: max, sqrt

  ! Input Variables
  character(len=50), intent(in) :: runtype ! The case that is being run

  integer, intent(in) :: sfctype

  real, intent(in) ::  & 
    thlm_sfc,  & ! thlm at (2)         [m/s]
    rtm_sfc,   & ! rtm at (2)          [kg/kg]
    T_sfc,      & ! Temperature         [K]
    p_sfc,      & ! Surface pressure    [Pa]
    ubar,      & ! This is root (u^2 + v^2), per ATEX and RICO spec.
    lowest_level ! This is z at the lowest above-ground model level.  [m]

  ! Output variables
  real, intent(out) ::  & 
    wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
    wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
    ustar           ! surface friction velocity [m/s]

  ! Constants
  real, parameter :: & 
  !  rho_sfc_flux = 1.0, &
  !  C_10    = 0.0013,    & ! Drag coefficient, defined by ATEX specification
    C_h_20  = 0.001094,  & ! Drag coefficient, defined by RICO 3D specification
    C_q_20  = 0.001133,  & ! Drag coefficient, defined by RICO 3D specification
    z0      = 0.00015      ! Roughness length, defined by ATEX specification

  ! Internal variables
  real :: &
    Ch,   &                ! This is C_h_20 scaled to the height of the lowest model level.
    Cq,   &                ! This is C_q_20 scaled to the height of the lowest model level.
    exner_sfc ! Value of exner at the surface [-]
    
  !--------------BEGIN CODE---------------------

  ! Calculate exner_sfc based on p_sfc.
  exner_sfc = ( p_sfc / p0 )**kappa

  ! Just set ustar = 0.3
  ustar = 0.3

  ! Get rid of a compiler warning
  if( runtype == "anything" .or. thlm_sfc == 1 .or. rtm_sfc == 1 .or. &
          exner_sfc == 1 .or. p_sfc == 1 .or. T_sfc == 1) then
      ustar = 0.3
  end if

  !--------------------------------------------------------------------------------
  ! Email way
  ! wprtp = value_from_forcings_file_in_W_m**2 / ( rho_sfc_flux * Lv )
  ! wpthlp = value_from_forcings_file_in_W_m**2 / ( rho_sfc_flux * Cp )

  !lhflx(1) = 0.001 * ubar * rho_sfc_flux * Lv * ( sat_mixrat_liq( p_sfc, T_sfc ) - & 
  !                                                sat_mixrat_liq( p_sfc, T_in_K ) * 0.8 )
  !shflx(1) = 0.001 * ubar * rho_sfc_flux * Cp * ( T_sfc - T_in_K )

  ! If this is the S6 case, fudge the values of the fluxes using values from the forcings
  !if ( runtype == "cloud_feedback_s6" .or. runtype == "cloud_feedback_s6_p2k" ) then
  !    wprtp_sfc = lhflx(1) / ( 1.0 * Lv )
  !    wpthlp_sfc = shflx(1) / ( 1.0 * Cp )
  !else
  !    wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, sat_mixrat_liq( p_sfc, T_sfc ) )
  !    wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, & 
  !                                     T_sfc, exner_sfc )
  !end if

  !
  ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Ch   = C_h_20 * ((log(20/z0))/(log(lowest_level/z0))) * & 
         ((log(20/z0))/(log(lowest_level/z0)))
  ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Cq   = C_q_20 * ((log(20/z0))/(log(lowest_level/z0))) * & 
         ((log(20/z0))/(log(lowest_level/z0)))
 
  if ( sfctype == 1 ) then
    wprtp_sfc = compute_wprtp_sfc( Cq, ubar, rtm_sfc, sat_mixrat_liq( p_sfc, T_sfc ) )
    wpthlp_sfc = compute_wpthlp_sfc( Ch, ubar, thlm_sfc, & 
                                     T_sfc, exner_sfc )

  end if


  return
  end subroutine cloud_feedback_sfclyr

end module cloud_feedback
