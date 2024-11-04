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
  subroutine fire_sfclyr( ngrdcol, time, ubar, p_sfc, & 
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
  integer, intent(in) :: &
    ngrdcol

  real(time_precision), intent(in) :: &
    time   ! current time [s]

  real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  & 
    ubar,    & ! mean sfc wind speed                           [m/s]
    p_sfc,    & ! Surface pressure                              [Pa]
    thlm_sfc,& ! theta_l at first model layer                  [K]
    rtm_sfc, & ! Total water mixing ratio at first model layer [kg/kg]
    exner_sfc

  integer, intent(in) :: &
    saturation_formula ! Integer that stores the saturation formula to be used

  ! Output Variables
  real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
    wpthlp_sfc, &   ! surface thetal flux        [K m/s]
    wprtp_sfc,  &   ! surface moisture flux      [kg/kg m/s]
    ustar,      &
    T_sfc           ! Surface temperature        [K]

  ! Local Variable
  real( kind = core_rknd ), dimension(ngrdcol) :: & 
    rsat, &
    Cz   ! Coefficient

  real( kind = core_rknd ) :: & 
    time_frac, & ! the time fraction used for interpolation
    T_sfc_interp

  integer :: &
    before_time, after_time, & ! the time indexes used for interpolation
    i

  !--------------BEGIN CODE---------------

  !$acc enter data create( rsat, Cz )

  ! Interpolate variables from time_dependent_input

  call time_select( time, size(time_sfc_given), time_sfc_given, &
                    before_time, after_time, time_frac )

  T_sfc_interp = linear_interp_factor( time_frac, T_sfc_given(after_time), &
                                       T_sfc_given(before_time) )

  !$acc parallel loop gang vector default(present)
  do i = 1, ngrdcol
    T_sfc(i) = T_sfc_interp
    Cz(i) = 0.0013_core_rknd
    ustar(i) = 0.3_core_rknd

    rsat(i) = sat_mixrat_liq( p_sfc(i), T_sfc(i), saturation_formula )
  end do

  ! Compute wpthlp_sfc and wprtp_sfc
  call compute_wpthlp_sfc( ngrdcol, Cz, ubar, thlm_sfc, T_sfc, exner_sfc, &
                           wpthlp_sfc ) 

  call compute_wprtp_sfc( ngrdcol, Cz, ubar, rtm_sfc, rsat, &
                          wprtp_sfc )

  !$acc exit data delete( rsat, Cz )

  return

  end subroutine fire_sfclyr

!----------------------------------------------------------------------
end module fire
