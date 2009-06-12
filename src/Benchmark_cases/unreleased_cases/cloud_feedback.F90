!----------------------------------------------------------------------
! $Id$
module cloud_feedback

!       Description:
!       Contains subroutines for the Cloud Feedback cases.
!----------------------------------------------------------------------

implicit none

private ! Default Scope

public :: cloud_feedback_tndcy, cloud_feedback_sfclyr, cloud_feedback_init

private :: divT, divq, press, ndiv, lhflx, shflx

! Constant Parameters
integer, parameter :: &  
  ndiv = 26,          & ! Number of values in the forcing files
  per_line = 1          ! Number of values per line in the forcing files

! Forcing arrays
real, dimension(ndiv) :: divT  ! Horizontal large scale temp. forcing
real, dimension(ndiv) :: divq  ! Horizontal large scale water vapor forcing
real, dimension(ndiv) :: press ! Pressure levels
real, dimension(1)    :: lhflx ! Surface latent heat flux
real, dimension(1)    :: shflx ! Surface sensible heat flux

contains

!----------------------------------------------------------------------
subroutine cloud_feedback_tndcy( exner, p_in_Pa, & 
                                 thlm_forcing, rtm_forcing, & 
                                 sclrm_forcing, edsclrm_forcing )
!       Description:
!       Sets up the tendency information for the cloud feeback case

!       References:
!----------------------------------------------------------------------

use grid_class, only: gr !  Variable(s)

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

use interpolation, only: zlinterp_fnc ! Procedure(s)

use stats_precision, only: time_precision ! Variable(s)

use array_index, only:  & 
    iisclr_thl, iisclr_rt ! Variable(s)

implicit none

! Input
real, intent(in), dimension(gr%nnzp) :: &
  exner,        & ! Exner function                [-]
  p_in_Pa         ! Pressure                      [Pa]

! Output Variables
real, intent(out), dimension(gr%nnzp) :: & 
  thlm_forcing, & ! Liquid water potential temperature tendency  [K/s]
  rtm_forcing     ! Total water mixing ratio tendency            [kg/kg/s]

real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
  sclrm_forcing ! Passive scalar forcing [units vary]

real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
  edsclrm_forcing ! Passive eddy-scalar forcing [units vary]

! Horizontal large scale temp. forcing
thlm_forcing = zlinterp_fnc( gr%nnzp, ndiv, -p_in_Pa, -press, divT ) / exner

! Large scale advective moisture tendency
rtm_forcing = zlinterp_fnc( gr%nnzp, ndiv, -p_in_Pa, -press, divq )

! Test scalars with thetal and rt if desired
if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

return
end subroutine cloud_feedback_tndcy

!----------------------------------------------------------------------
subroutine cloud_feedback_sfclyr( runtype, time, p_in_Pa, rho0, & 
                                  lowestlevel, thlm_sfc, rtm_sfc, um_sfc, vm_sfc,  &
                                  exner_sfc, psfc, Tsfc, rcm, & 
                                  upwp_sfc, vpwp_sfc, & 
                                  wpthlp_sfc, wprtp_sfc, ustar, & 
                                  wpsclrp_sfc, wpedsclrp_sfc )

!       Description:
!       Sets up surface information for the cloud feedback case

!       References:
!----------------------------------------------------------------------

use constants, only: pi, grav, Lv, Cp ! Variable(s)

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use saturation, only: sat_mixrat_liq ! Variable(s)

use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

use stats_precision, only: time_precision ! Variable(s)

use array_index, only:  & 
    iisclr_thl, iisclr_rt ! Variable(s)

use surface_flux, only: &
    compute_ubar, compute_momentum_flux, compute_wprtp_sfc, compute_wpthlp_sfc

!use T_in_K_mod, only: &
!    thlm2T_in_K ! Procedure

implicit none

intrinsic :: max, sqrt

! Input Variables
real(kind=time_precision), intent(in) ::  & 
  time      ! Current time        [s] 

character(len=50), intent(in) :: runtype ! The case that is being run

real, intent(in) ::  & 
  p_in_Pa,   & ! Pressure            [Pa] 
  rho0,      & ! Density at zm=1     [kg/m^3] 
  thlm_sfc,  & ! thlm at (2)         [m/s]
  rtm_sfc,   & ! rtm at (2)          [kg/kg]
  Tsfc,      & ! Temperature         [K]
  psfc,      &
  rcm,       & ! Cloud water mixing ratio      [kg/kg]
  exner_sfc, & ! Exner function      [-]
  lowestlevel,   & ! This is z at the lowest above-ground model level.  [m]
  um_sfc,    & ! um at (2)           [m/s]
  vm_sfc       ! vm at (2)           [m/s]

! Output variables
real, intent(out) ::  & 
  upwp_sfc,     & ! u'w' at (1)      [m^2/s^2]
  vpwp_sfc,     & ! v'w'at (1)       [m^2/s^2]
  wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
  wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
  ustar           ! surface friction velocity [m/s]

real, intent(out), dimension(sclr_dim) ::  & 
  wpsclrp_sfc     ! Passive scalar surface flux      [units m/s] 

real, intent(out), dimension(edsclr_dim) ::  & 
  wpedsclrp_sfc   ! Passive eddy-scalar surface flux [units m/s]

! Constants
real, parameter :: & 
  C_10    = 0.0013, &     ! Drag coefficient, defined by ATEX specification
!  z0      = 0.00015,   & ! Roughness length, defined by ATEX specification
  rho_sfc_flux = 1.0

! Internal variables
real :: & 
  ubar 
!  T_in_K

!T_in_K = thlm2T_in_K( thlm_sfc, exner_sfc, rcm )

ubar = compute_ubar( um_sfc, vm_sfc )

! Just set ustar = 0.3
ustar = 0.3

!--------------------------------------------------------------------------------
! Old Style

!wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, & 
!                                 Tsfc, exner_sfc )
!wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, 0.8 * sat_mixrat_liq( psfc, Tsfc ) )

!call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
!                            upwp_sfc, vpwp_sfc )

!--------------------------------------------------------------------------------
! Email way
! wprtp = value_from_forcings_file_in_W_m**2 / ( rho_sfc_flux * Lv )
! wpthlp = value_from_forcings_file_in_W_m**2 / ( rho_sfc_flux * Cp )

! S11 is not affected much by using the calculated surface fluxes other then cf is lower (with 0.8)
! S6 crashes when when compute_wpthlp_sfc and compute_wprtp_sfc are used (with 0.8)
! S6 and S12 crash when the calculated surface fluxes are used to calculate wpthlp (0.8)
! S6 crashes with calculated surface fluxes (without 0.8)

!lhflx(1) = 0.001 * ubar * rho_sfc_flux * Lv * ( sat_mixrat_liq( psfc, Tsfc ) - & 
!                                                sat_mixrat_liq( psfc, T_in_K ) * 0.8 )
!shflx(1) = 0.001 * ubar * rho_sfc_flux * Cp * ( Tsfc - T_in_K )

! If this is the S6 case, fudge the values of the fluxes using values from the forcings
if ( runtype == "cloud_feedback_s6" .or. runtype == "cloud_feedback_s6_p2k" ) then
    wprtp_sfc = lhflx(1) / ( rho_sfc_flux * Lv )
    wpthlp_sfc = shflx(1) / ( rho_sfc_flux * Cp )
else
    wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, sat_mixrat_liq( psfc, Tsfc ) )
    wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, & 
                                     Tsfc, exner_sfc )
end if

call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                            upwp_sfc, vpwp_sfc )

! Let passive scalars be equal to rt and theta_l for now
if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc

return
end subroutine cloud_feedback_sfclyr

!----------------------------------------------------------------
subroutine cloud_feedback_init( iunit, file_path )
!
!       Description:
!       This subroutine initializes the module by reading in forcing
!       data used in the tndcy subroutine.
!----------------------------------------------------------------

  use file_functions, only: file_read_1d ! Procedure(s)

  implicit none

  integer, intent(in) :: iunit ! File unit number

  character(len=*), intent(in) :: &
    file_path ! Path to the forcing files

  call file_read_1d( iunit, & 
    file_path//'cloud_feedback_divT.dat', & 
    ndiv, per_line, divT )

  call file_read_1d( iunit, & 
    file_path//'cloud_feedback_divq.dat', & 
    ndiv, per_line, divq )

  call file_read_1d( iunit, & 
    file_path//'cloud_feedback_press.dat', & 
    ndiv, per_line, press )

  call file_read_1d( iunit, & 
    file_path//'cloud_feedback_press.dat', & 
    ndiv, per_line, press )

  call file_read_1d( iunit, & 
    file_path//'cloud_feedback_shflx.dat', & 
    1, 1, shflx )

  call file_read_1d( iunit, & 
    file_path//'cloud_feedback_lhflx.dat', & 
    1, 1, lhflx )

  return 
end subroutine cloud_feedback_init

end module cloud_feedback
