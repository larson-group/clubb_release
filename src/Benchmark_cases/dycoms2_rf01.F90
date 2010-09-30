!----------------------------------------------------------------------
! $Id$
module dycoms2_rf01

!       Description:
!       Contains subroutines for the DYCOMS II RF01 case.
!----------------------------------------------------------------------
  implicit none

  public :: dycoms2_rf01_tndcy, dycoms2_rf01_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine dycoms2_rf01_tndcy( thlm_forcing, rtm_forcing, &
                                 sclrm_forcing, edsclrm_forcing )
! Description:
!   Subroutine to set theta and water tendencies for DYCOMS RF01 case.

! References:
!----------------------------------------------------------------------

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variables(s)

    use grid_class, only: gr

    implicit none

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) ::  & 
      thlm_forcing,  & ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing      ! Total water mixing ratio tendency            [kg/kg/s]

    real, intent(out), dimension(gr%nnzp, sclr_dim) :: & 
      sclrm_forcing   ! Passive scalar tendency         [units/s]

    real, intent(out), dimension(gr%nnzp, edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar tendency    [units/s]

    thlm_forcing = 0.
    rtm_forcing  = 0.

    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

    return
  end subroutine dycoms2_rf01_tndcy
  
  !======================================================================
  subroutine dycoms2_rf01_sfclyr( sfctype, T_sfc, p_sfc,  & 
                                    exner_sfc, ubar, & 
                                    thlm_sfc, rtm_sfc, rho_zm_sfc, &
                                    wpthlp_sfc, wprtp_sfc, ustar )
  ! Description:
  !   This subroutine computes surface fluxes of
  !   heat and moisture according to GCSS DYCOMS II RF 01 specifications

  ! References:
  !   None
  !----------------------------------------------------------------------
  use constants_clubb, only: Cp, fstderr, Lv ! Variable(s)

  use saturation, only: sat_mixrat_liq ! Variable(s)

  use surface_flux, only: compute_wpthlp_sfc, compute_wprtp_sfc ! Procedure(s)

  implicit none

  ! Input variables
  integer, intent(in) :: &
    sfctype
  real, intent(in) ::  &
    T_sfc,      & ! Surface temperature                           [K]
    p_sfc,      & ! Surface pressure                              [Pa]
    exner_sfc, & ! Exner function                                [-]
    ubar,      & ! mean sfc wind speed                           [m/s]
    thlm_sfc,  & ! theta_l at first model layer                  [K]
    rtm_sfc,   & ! Total water mixing ratio at first model layer [kg/kg]
    rho_zm_sfc   ! Density at the surface                        [kg/m^3]

  ! Output variables
  real, intent(out) ::  & 
    wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
    wprtp_sfc, &      ! w'rt' surface flux        [(m kg)/(kg s)]
    ustar
    
  ! Local Variable
  real :: & 
    Cd  ! Coefficient

  !-----------------BEGIN CODE-----------------------

  Cd = 0.0011

  ustar = 0.25

  ! Compute heat and moisture fluxes
  if ( sfctype == 0 ) then

    wpthlp_sfc =  15.0 / ( rho_zm_sfc * Cp )
    wprtp_sfc  = 115.0 / ( rho_zm_sfc * Lv )

  else if ( sfctype == 1 ) then

    wpthlp_sfc = compute_wpthlp_sfc( Cd, ubar, thlm_sfc, T_sfc, exner_sfc )
    wprtp_sfc = compute_wprtp_sfc( Cd, ubar, rtm_sfc, sat_mixrat_liq( p_sfc, T_sfc ) )

  else  ! Undefined value for sfctype

    write(fstderr,*) "Invalid sfctype value = ", sfctype
    stop

  end if ! sfctype
  return
  end subroutine dycoms2_rf01_sfclyr

!----------------------------------------------------------------------
end module dycoms2_rf01
