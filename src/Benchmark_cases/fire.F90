!----------------------------------------------------------------------
! $Id$
module fire

  ! Description:
  !   Contains subroutines for the GCSS FIRE case.
  !
  ! References:
  !   Moeng, C.-H., and coauthors, 1996: Simulation of a
  !   stratocumulus-topped PBL: Intercomparison among different 
  !   numerical codes. Bull. Amer. Meteor. Soc., 77, 261-278.
  !   ftp://eos.atmos.washington.edu/pub/breth/papers/1996/GCSS1-Moeng.pdf
  !----------------------------------------------------------------------

  implicit none

  public :: fire_sfclyr

  private ! Default Scope

  contains

  !======================================================================
  subroutine fire_sfclyr( time, ubar, p_sfc, & 
                          thlm_sfc, rtm_sfc, exner_sfc, & 
                          saturation_formula, &
                          wpthlp_sfc, wprtp_sfc, ustar, T_sfc )
                                          
  ! Description:
  !   This subroutine computes surface fluxes of heat and moisture 
  !   using aerodynamic formulas.

  ! References:
  !   Moeng, C.-H., and coauthors, 1996: Simulation of a
  !   stratocumulus-topped PBL: Intercomparison among different 
  !   numerical codes. Bull. Amer. Meteor. Soc., 77, 261-278.
  !   ftp://eos.atmos.washington.edu/pub/breth/papers/1996/GCSS1-Moeng.pdf
  !------------------------------------------------------------------------

  use saturation, only: sat_mixrat_liq ! Procedure(s)
 
  use sfc_flux, only: compute_wprtp_sfc, compute_wpthlp_sfc

  use clubb_precision, only: time_precision, core_rknd ! Variable(s)

  use interpolation, only: linear_interp_factor ! Procedure(s)

  use time_dependent_input, only: time_sfc_given, T_sfc_given, & ! Variable(s)
                                  time_select ! Procedure(s)

  implicit none

  ! Input Variables
  real(time_precision), intent(in) :: &
    time   ! current time [s]

  real( kind = core_rknd ), intent(in) ::  & 
    ubar,    & ! mean sfc wind speed                           [m/s]
    p_sfc,    & ! Surface pressure                              [Pa]
    thlm_sfc,& ! theta_l at first model layer                  [K]
    rtm_sfc, & ! Total water mixing ratio at first model layer [kg/kg]
    exner_sfc

  integer, intent(in) :: &
    saturation_formula ! Integer that stores the saturation formula to be used

  ! Output Variables
  real( kind = core_rknd ), intent(out) ::  & 
    wpthlp_sfc, &   ! surface thetal flux        [K m/s]
    wprtp_sfc,  &   ! surface moisture flux      [kg/kg m/s]
    ustar,      &
    T_sfc           ! Surface temperature        [K]

  ! Local Variable
  real( kind = core_rknd ) :: & 
    Cz, &  ! Coefficient
    time_frac ! the time fraction used for interpolation

  integer :: &
    before_time, after_time ! the time indexes used for interpolation

  !--------------BEGIN CODE---------------

  ! Interpolate variables from time_dependent_input

  call time_select( time, size(time_sfc_given), time_sfc_given, &
                    before_time, after_time, time_frac )

  T_sfc = linear_interp_factor( time_frac, T_sfc_given(after_time), &
                                    T_sfc_given(before_time) )

  ! Compute wpthlp_sfc and wprtp_sfc

  Cz = 0.0013_core_rknd

  ustar = 0.3_core_rknd

  wpthlp_sfc = compute_wpthlp_sfc ( Cz, ubar, thlm_sfc, T_sfc, exner_sfc )
  wprtp_sfc = compute_wprtp_sfc( Cz, ubar, rtm_sfc, sat_mixrat_liq( p_sfc, T_sfc, saturation_formula ) )

  return
  end subroutine fire_sfclyr

!----------------------------------------------------------------------
end module fire
