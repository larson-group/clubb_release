!----------------------------------------------------------------------
! $Id$
module dycoms2_rf01

!       Description:
!       Contains subroutines for the DYCOMS II RF01 case.
!
!       References:
!       <http://www.atmos.ucla.edu/~bstevens/dycoms/rf01/rf01.html>
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
!   <http://www.atmos.ucla.edu/~bstevens/dycoms/rf01/rf01.html>
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
  subroutine dycoms2_rf01_sfclyr( time, sfctype, p_sfc,  & 
                                    exner_sfc, ubar, & 
                                    thlm_sfc, rtm_sfc, rho_sfc, &
                                    wpthlp_sfc, wprtp_sfc, ustar, T_sfc )
  ! Description:
  !   This subroutine computes surface fluxes of
  !   heat and moisture according to GCSS DYCOMS II RF 01 specifications

  ! References:
  !   <http://www.atmos.ucla.edu/~bstevens/dycoms/rf01/rf01.html>
  !----------------------------------------------------------------------
  use constants_clubb, only: fstderr ! Variable(s)

  use saturation, only: sat_mixrat_liq ! Variable(s)

  use surface_flux, only: compute_wpthlp_sfc, compute_wprtp_sfc, &
                          convert_SH_to_km_s, convert_LH_to_m_s ! Procedure(s)

  use time_dependent_input, only: SH_given, LH_given, time_sfc_given,& ! Variable(s)
                                  T_sfc_given, &
                                  time_select ! Procedure(s)

  use stats_precision, only: time_precision ! Variable(s)
  
  use interpolation, only: factor_interp ! Procedure(s)

  implicit none

  ! Input variables
  real(time_precision), intent(in) :: &
    time ! The current time [s]
  integer, intent(in) :: &
    sfctype
  real, intent(in) ::  &
    p_sfc,      & ! Surface pressure                              [Pa]
    exner_sfc, & ! Exner function                                [-]
    ubar,      & ! mean sfc wind speed                           [m/s]
    thlm_sfc,  & ! theta_l at first model layer                  [K]
    rtm_sfc,   & ! Total water mixing ratio at first model layer [kg/kg]
    rho_sfc   ! Density at the surface                        [kg/m^3]

  ! Output variables
  real, intent(out) ::  & 
    wpthlp_sfc,  &  ! w'theta_l' surface flux   [(m K)/s]
    wprtp_sfc,   &  ! w'rt' surface flux        [(m kg)/(kg s)]
    ustar,       &
    T_sfc           ! Surface temperature       [K]
    
  ! Local Variable
  real, parameter :: & 
    Cd = 0.0011   ! Coefficient
    
  real :: &
    SH, &  ! Sensible heat flux
    LH, &  ! Latent heat flux
    time_frac ! The time fraction used for interpolation

  integer :: &
    before_time, after_time ! The times used for interpolation

  !-----------------BEGIN CODE-----------------------

  ustar = 0.25

  call time_select( time, size(time_sfc_given), time_sfc_given, &
                    before_time, after_time, time_frac )

  SH = factor_interp( time_frac, SH_given(after_time), SH_given(before_time) )
  LH = factor_interp( time_frac, LH_given(after_time), LH_given(before_time) )
  T_sfc = factor_interp( time_frac, T_sfc_given(after_time), &
                                    T_sfc_given(before_time) )

  ! Compute heat and moisture fluxes
  if ( sfctype == 0 ) then

    wpthlp_sfc = convert_SH_to_km_s( SH, rho_sfc )
    wprtp_sfc  = convert_LH_to_m_s( LH, rho_sfc )

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
