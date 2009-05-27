!----------------------------------------------------------------------
! $Id$
module fire

!       Description:
!       Contains subroutines for the GCSS FIRE case.
!----------------------------------------------------------------------

implicit none

public :: fire_tndcy, sfc_momentum_fluxes, sfc_thermo_fluxes

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

use atex_cloud_rad, only: cloud_rad ! Procedure(s)

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

!------------------------------------------------------------------------
subroutine sfc_momentum_fluxes( um_sfc, vm_sfc, &
                                upwp_sfc, vpwp_sfc, ustar )

!       Description:
!       This subroutine computes surface momentum fluxes using aerodynamic
!       formulas.

!       References:
!       None
!------------------------------------------------------------------------

use surface_flux, only: compute_ubar, compute_momentum_flux

implicit none

! External
intrinsic :: sqrt

! Input variables
real, intent(in) ::  & 
  um_sfc,  & ! u wind first level above ground    [m/s]
  vm_sfc     ! v wind first level above ground    [m/s]

! Output Variables
real, intent(out) ::  & 
  upwp_sfc, & ! sfc u momentum flux (m^2/s^2)
  vpwp_sfc, & ! sfc v momentum flux (m^2/s^2)
  ustar       ! surface friction velocity [m/s]

! Local Variables
real :: ubar ! total wind speed above ground

! Declare the value of ustar
ustar = 0.3

! Computes fluxes

ubar = compute_ubar( um_sfc, vm_sfc )

call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                            upwp_sfc, vpwp_sfc )

return
end subroutine sfc_momentum_fluxes

!------------------------------------------------------------------------
subroutine sfc_thermo_fluxes( um_sfc, vm_sfc, &
                              Tsfc, psfc, & 
                              thlair, rtair, exner_sfc, & 
                              wpthlp_sfc, wprtp_sfc, & 
                              wpsclrp_sfc, & 
                              wpedsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of heat and moisture 
!       using aerodynamic formulas.

!       References:
!       None
!------------------------------------------------------------------------

use constants, only: kappa, p0 ! Variable(s)

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use saturation, only: sat_mixrat_liq ! Procedure(s)

use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl

use surface_flux, only: compute_ubar, compute_wprtp_sfc, compute_wpthlp_sfc

implicit none

! External
intrinsic :: sqrt

! Parameter
real, parameter :: C = 1.3e-3

! Input Variables
real, intent(in) ::  & 
  um_sfc,  & ! u wind                        [m/s]
  vm_sfc,  & ! u wind                        [m/s]
  Tsfc,    & ! Surface temperature           [K]
  psfc,    & ! Surface pressure              [Pa]
  thlair,  & ! theta_l at first model layer  [K]
  rtair,   & ! rt at first model layer       [kg/kg]
  exner_sfc

! Output Variables
real, intent(out) ::  & 
  wpthlp_sfc, & ! surface thetal flux        [K m/s]
  wprtp_sfc     ! surface moisture flux      [kg/kg m/s]

real, intent(out), dimension(sclr_dim) ::  & 
  wpsclrp_sfc       ! scalar surface flux            [units m/s]

real, intent(out), dimension(edsclr_dim) ::  & 
  wpedsclrp_sfc     ! eddy-scalar surface flux       [units m/s]

! Local Variables
real :: ubar  ! Total wind speed above ground

! Compute fluxes
ubar = compute_ubar( um_sfc, vm_sfc )

wpthlp_sfc = compute_wpthlp_sfc ( C, ubar, thlair, Tsfc, exner_sfc )
wprtp_sfc = compute_wprtp_sfc( C, ubar, rtair, sat_mixrat_liq( psfc, Tsfc ) )

! Let passive scalars be equal to rt and theta_l for now
if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc

return
end subroutine sfc_thermo_fluxes

end module fire
