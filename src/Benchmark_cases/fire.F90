!----------------------------------------------------------------------
! $Id$
module fire

  !       Description:
  !       Contains subroutines for the GCSS FIRE case.
  !----------------------------------------------------------------------

  implicit none

  public :: fire_tndcy, fire_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine fire_tndcy & 
             ( rho, rcm, exner,  & 
               Frad, radht,  & 
               thlm_forcing, rtm_forcing, & 
               sclrm_forcing, edsclrm_forcing )
  !       Description:
  !       Subroutine to large-scale subsidence for FIRE case. Calls
  !       cloud_rad for computing radiation

  !       References:
  !       None
  !----------------------------------------------------------------------

  use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

  use parameters_radiation, only: rad_scheme  ! Variable(s)

  use grid_class, only: gr ! Variable(s)

  use grid_class, only: zt2zm ! Procedure(s)

  use cloud_rad_module, only: cloud_rad ! Procedure(s)

  use stats_precision, only: time_precision ! Variable(s)

  use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)
   
  use stats_type, only: stat_update_var ! Procedure(s)

  use stats_variables, only: zt, iradht_LW, l_stats_samp ! Variable(s)

  implicit none

  ! Input Variables
  real, intent(in), dimension(gr%nnzp) :: & 
    rho,   & ! Density                         [kg/m^3]
    rcm,   & ! Liquid water mixing ratio       [kg/kg]
    exner    ! Exner function                  [-]

  ! Output Variables
  real, intent(out), dimension(gr%nnzp) :: & 
    Frad,         & ! Radiative flux                   [W/m^2]
    radht,        & ! Radiative heating rate           [K/s]
    thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
    rtm_forcing     ! Total water mixing ratio tendency [kg/kg/s]

  real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
    sclrm_forcing ! Passive scalar tendency [units/s]

  real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
    edsclrm_forcing ! Passive scalar tendency [units/s]


  ! Radiative theta-l tendency is computed interactively elsewhere

  thlm_forcing = 0.0

  ! Large scale advective moisture tendency

  rtm_forcing = 0.0

  ! Use cloud_rad to compute radiation
  if ( trim( rad_scheme ) == "simplified" ) then

    call cloud_rad( rho, rcm, exner, Frad, radht, thlm_forcing )

    if ( l_stats_samp ) then
      call stat_update_var( iradht_LW, radht, zt )
    end if

  end if

  ! Test scalars with thetal and rt if desired
  if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
  if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

  if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
  if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

  return
  end subroutine fire_tndcy

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
