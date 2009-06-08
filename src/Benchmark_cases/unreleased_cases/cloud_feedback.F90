!----------------------------------------------------------------------
! $Id:$
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
subroutine cloud_feedback_tndcy( time, rcm, exner, & 
                                 p_in_Pa, & 
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
real(kind=time_precision), intent(in) :: time ! Model time [s]

real, intent(in), dimension(gr%nnzp) :: &
  rcm,          & ! Liquid water mixing ratio     [kg/kg]
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
subroutine cloud_feedback_sfclyr( time, p_in_Pa, rho0, lowestlevel, & 
                                  thlm_sfc, rtm_sfc, um_sfc, vm_sfc,  &
                                  exner_sfc, psfc, Tsfc, & 
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

implicit none

intrinsic :: max, sqrt

! Input Variables
real(kind=time_precision), intent(in) ::  & 
  time      ! Current time        [s] 

real, intent(in) ::  & 
  p_in_Pa,   & ! Pressure            [Pa] 
  rho0,      & ! Density at zm=1     [kg/m^3] 
  thlm_sfc,  & ! thlm at (2)         [m/s]
  rtm_sfc,   & ! rtm at (2)          [kg/kg]
  Tsfc,      & ! Temperature         [K]
  psfc,      &
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
  C_10    = 0.0013       ! Drag coefficient, defined by ATEX specification
!  C_m_20  = 0.001229,  & ! Drag coefficient, defined by RICO 3D specification
!  C_h_20  = 0.001094,  & ! Drag coefficient, defined by RICO 3D specification
!  C_q_20  = 0.001133,  & ! Drag coefficient, defined by RICO 3D specification
!  z0      = 0.00015      ! Roughness length, defined by ATEX specification
!  rho_sfc_flux = 1.0

! Internal variables
real :: & 
  ubar
!  Cz,   & ! This is C_10 scaled to the height of the lowest model level.
!  Cm,   & ! This is C_m_20 scaled to the height of the lowest model level.
!  Ch,   & ! This is C_h_20 scaled to the height of the lowest model level.
!  Cq      ! This is C_q_20 scaled to the height of the lowest model level.

ubar = compute_ubar( um_sfc, vm_sfc )

! Just set ustar = 0.3
ustar = 0.3

! Modification in case lowest model level isn't at 10 m, from ATEX specification
!Cz   = C_10 * ((log(10/z0))/(log(lowestlevel/z0))) * & 
!       ((log(10/z0))/(log(lowestlevel/z0)))         
! Modification in case lowest model level isn't at 10 m, from ATEX specification
!Cm   = C_m_20 * ((log(20/z0))/(log(lowestlevel/z0))) * & 
!       ((log(20/z0))/(log(lowestlevel/z0)))             
! Modification in case lowest model level isn't at 10 m, from ATEX specification
!Ch   = C_h_20 * ((log(20/z0))/(log(lowestlevel/z0))) * & 
!       ((log(20/z0))/(log(lowestlevel/z0)))          
! Modification in case lowest model level isn't at 10 m, from ATEX specification
!Cq   = C_q_20 * ((log(20/z0))/(log(lowestlevel/z0))) * & 
!       ((log(20/z0))/(log(lowestlevel/z0)))

!--------------------------------------------------------------------------------
! ATEX Style
!wpthlp_sfc = compute_wpthlp_sfc( Cz, ubar, thlm_sfc, Tsfc, exner_sfc )
!wprtp_sfc = compute_wprtp_sfc( Cz, ubar, rtm_sfc, sat_mixrat_liq( psfc,Tsfc ) )
!call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
!                            upwp_sfc, vpwp_sfc )

!--------------------------------------------------------------------------------
! Rico Style
!wpthlp_sfc = compute_wpthlp_sfc( Ch, ubar, thlm_sfc, Tsfc, exner_sfc )
!wprtp_sfc  = compute_wprtp_sfc( Cq, ubar, rtm_sfc, 0.8 * sat_mixrat_liq( psfc,Tsfc ) )
!upwp_sfc   = -um_sfc * Cm * ubar  ! m^2 s^-2
!vpwp_sfc   = -vm_sfc * Cm * ubar  ! m^2 s^-2

!--------------------------------------------------------------------------------
! Old Style

wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, & 
                                 Tsfc, exner_sfc )
wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, sat_mixrat_liq( psfc, Tsfc ) )

call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                            upwp_sfc, vpwp_sfc )

!--------------------------------------------------------------------------------
! Email way
! wprtp = value_from_forcings_file_in_W_m**2 / ( rho_sfc_flux * Lv )
! wpthlp = value_from_forcings_file_in_W_m**2 / ( rho_sfc_flux * Cp )

!lhflx = 0.001 * ubar * rho_sfc_flux * Lv * ( sat_mixrat_liq( psfc, Tsfc ) - & 
!                                              sat_mixrat_liq( psfc, T_in_K ) * 0.8 )
!shflx = 0.001 * ubar * rho_sfc_flux * Cp * [ Tsfc - T_in_K ]

!wprtp_sfc = lhflx / ( rho_sfc_flux * Lv )
!wpthlp_sfc = shflx / ( rho_sfc_flux * Cp )

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
