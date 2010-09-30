!----------------------------------------------------------------------
! $Id$
module fire

  !       Description:
  !       Contains subroutines for the GCSS FIRE case.
  !----------------------------------------------------------------------

  implicit none

  public :: fire_sfclyr

  private ! Default Scope

  contains

  !======================================================================
  subroutine fire_sfclyr( ubar, T_sfc, p_sfc, & 
                          thlm_sfc, rtm_sfc, exner_sfc, & 
                          wpthlp_sfc, wprtp_sfc, ustar )
                                          
  !       Description:
  !       This subroutine computes surface fluxes of heat and moisture 
  !       using aerodynamic formulas.

  !       References:
  !       None
  !------------------------------------------------------------------------

  use saturation, only: sat_mixrat_liq ! Procedure(s)
  use surface_flux, only: compute_wprtp_sfc, compute_wpthlp_sfc

  implicit none

  ! Input Variables
  real, intent(in) ::  & 
    ubar,    & ! mean sfc wind speed                           [m/s]
    T_sfc,    & ! Surface temperature                           [K]
    p_sfc,    & ! Surface pressure                              [Pa]
    thlm_sfc,& ! theta_l at first model layer                  [K]
    rtm_sfc, & ! Total water mixing ratio at first model layer [kg/kg]
    exner_sfc

  ! Output Variables
  real, intent(out) ::  & 
    wpthlp_sfc, & ! surface thetal flux        [K m/s]
    wprtp_sfc, &     ! surface moisture flux      [kg/kg m/s]
    ustar
    
  ! Local Variable
  real :: & 
    Cz  ! Coefficient

  !--------------BEGIN CODE---------------

  Cz = 0.0013

  ustar = 0.3

  wpthlp_sfc = compute_wpthlp_sfc ( Cz, ubar, thlm_sfc, T_sfc, exner_sfc )
  wprtp_sfc = compute_wprtp_sfc( Cz, ubar, rtm_sfc, sat_mixrat_liq( p_sfc, T_sfc ) )

  return
  end subroutine fire_sfclyr

!----------------------------------------------------------------------
end module fire
