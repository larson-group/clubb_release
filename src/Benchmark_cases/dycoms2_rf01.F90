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
  subroutine dycoms2_rf01_tndcy( gr, sclr_dim, edsclr_dim, sclr_idx, &
                                 thlm_forcing, rtm_forcing, &
                                 sclrm_forcing, edsclrm_forcing )
! Description:
!   Subroutine to set theta and water tendencies for DYCOMS RF01 case.

! References:
!   <http://www.atmos.ucla.edu/~bstevens/dycoms/rf01/rf01.html>
!----------------------------------------------------------------------

    use array_index, only: &
      sclr_idx_type

    use grid_class, only: &
      grid

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    !--------------------- Input Variables ---------------------
    type(grid), intent(in) :: &
      gr

    integer, intent(in) :: &
      sclr_dim, & 
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(gr%nzt) ::  & 
      thlm_forcing,  & ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing      ! Total water mixing ratio tendency            [kg/kg/s]

    real( kind = core_rknd ), intent(out), dimension(gr%nzt, sclr_dim) :: & 
      sclrm_forcing   ! Passive scalar tendency         [units/s]

    real( kind = core_rknd ), intent(out), dimension(gr%nzt, edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar tendency    [units/s]

    thlm_forcing = 0._core_rknd
    rtm_forcing  = 0._core_rknd


    ! Test scalars with thetal and rt if desired
    if ( sclr_idx%iisclr_thl > 0 ) sclrm_forcing(:,sclr_idx%iisclr_thl) = thlm_forcing
    if ( sclr_idx%iisclr_rt  > 0 ) sclrm_forcing(:,sclr_idx%iisclr_rt)  = rtm_forcing

    if ( sclr_idx%iiedsclr_thl > 0 ) edsclrm_forcing(:,sclr_idx%iiedsclr_thl) = thlm_forcing
    if ( sclr_idx%iiedsclr_rt  > 0 ) edsclrm_forcing(:,sclr_idx%iiedsclr_rt)  = rtm_forcing

    return
  end subroutine dycoms2_rf01_tndcy
  
  !======================================================================
  subroutine dycoms2_rf01_sfclyr( time, sfctype, p_sfc,  & 
                                  exner_sfc, ubar, & 
                                  thlm_sfc, rtm_sfc, rho_sfc, &
                                  saturation_formula, &
                                  wpthlp_sfc, wprtp_sfc, ustar, T_sfc )
  ! Description:
  !   This subroutine computes surface fluxes of
  !   heat and moisture according to GCSS DYCOMS II RF 01 specifications

  ! References:
  !   <http://www.atmos.ucla.edu/~bstevens/dycoms/rf01/rf01.html>
  !----------------------------------------------------------------------
  use constants_clubb, only: fstderr ! Variable(s)

  use saturation, only: sat_mixrat_liq ! Variable(s)

  use sfc_flux, only: compute_wpthlp_sfc, compute_wprtp_sfc, &
                          convert_sens_ht_to_km_s, convert_latent_ht_to_m_s ! Procedure(s)

  use time_dependent_input, only: sens_ht_given, latent_ht_given, time_sfc_given,& ! Variable(s)
                                  T_sfc_given, &
                                  time_select ! Procedure(s)

  use clubb_precision, only: time_precision, core_rknd ! Variable(s)
  
  use interpolation, only: linear_interp_factor ! Procedure(s)

  implicit none

  ! Input variables
  real(time_precision), intent(in) :: &
    time ! The current time [s]
  integer, intent(in) :: &
    sfctype
  real( kind = core_rknd ), intent(in) ::  &
    p_sfc,      & ! Surface pressure                              [Pa]
    exner_sfc, & ! Exner function                                [-]
    ubar,      & ! mean sfc wind speed                           [m/s]
    thlm_sfc,  & ! theta_l at first model layer                  [K]
    rtm_sfc,   & ! Total water mixing ratio at first model layer [kg/kg]
    rho_sfc   ! Density at the surface                        [kg/m^3]

  integer, intent(in) :: &
    saturation_formula ! Integer that stores the saturation formula to be used

  ! Output variables
  real( kind = core_rknd ), intent(out) ::  & 
    wpthlp_sfc,  &  ! w'theta_l' surface flux   [(m K)/s]
    wprtp_sfc,   &  ! w'rt' surface flux        [(m kg)/(kg s)]
    ustar,       &
    T_sfc           ! Surface temperature       [K]
    
  ! Local Variable
  real( kind = core_rknd ), parameter :: & 
    Cd = 0.0011_core_rknd   ! Coefficient
    
  real( kind = core_rknd ) :: &
    sens_ht, &  ! Sensible heat flux
    latent_ht, &  ! Latent heat flux
    time_frac ! The time fraction used for interpolation

  integer :: &
    before_time, after_time ! The times used for interpolation

  !-----------------BEGIN CODE-----------------------

  ustar = 0.25_core_rknd

  call time_select( time, size(time_sfc_given), time_sfc_given, &
                    before_time, after_time, time_frac )

  sens_ht = linear_interp_factor( time_frac, sens_ht_given(after_time), sens_ht_given(before_time) )
  latent_ht = linear_interp_factor( time_frac, latent_ht_given(after_time), &
                                    latent_ht_given(before_time) )
  T_sfc = linear_interp_factor( time_frac, T_sfc_given(after_time), &
                                    T_sfc_given(before_time) )

  ! Compute heat and moisture fluxes
  if ( sfctype == 0 ) then

    wpthlp_sfc = convert_sens_ht_to_km_s( sens_ht, rho_sfc )
    wprtp_sfc  = convert_latent_ht_to_m_s( latent_ht, rho_sfc )

  else if ( sfctype == 1 ) then

    wpthlp_sfc = compute_wpthlp_sfc( Cd, ubar, thlm_sfc, T_sfc, exner_sfc )
    wprtp_sfc = compute_wprtp_sfc( Cd, ubar, rtm_sfc, sat_mixrat_liq( p_sfc, &
                                                                      T_sfc, saturation_formula ) )

  else  ! Undefined value for sfctype

    write(fstderr,*) "Invalid sfctype value = ", sfctype
    error stop

  end if ! sfctype
  return
  end subroutine dycoms2_rf01_sfclyr

!----------------------------------------------------------------------
end module dycoms2_rf01
