!----------------------------------------------------------------------
! $Id$
module cobra
  !       Description:
  !       Contains subroutines for the COBRA CO2 case.
  !----------------------------------------------------------------------

  implicit none

  private ! Default Scope

  public :: cobra_sfclyr

  contains

  !-----------------------------------------------------------------------
  subroutine cobra_sfclyr( sclr_dim, edsclr_dim, sclr_idx, &
                           time, z, rho_sfc, thlm_sfc, ubar, & 
                           wpthlp_sfc, wprtp_sfc, ustar, & 
                           wpsclrp_sfc, wpedsclrp_sfc, T_sfc )

  !       Description:
  !       This subroutine computes surface fluxes of horizontal momentum,
  !       heat and moisture according to the format used for the GCSS ARM 
  !       case.

  !       Notes:
  !       The data has been altered so it can be used for the COBRA CO2 
  !       case.
  !-----------------------------------------------------------------------

  use constants_clubb, only: &
    grav ! Variable(s)

  use interpolation, only: &
    linear_interp_factor

  use clubb_precision, only: &
    time_precision, & ! Variable(s)
    core_rknd 

  use diag_ustar_module, only: &
    diag_ustar ! Variable(s)

  use sfc_flux, only: &
    convert_sens_ht_to_km_s, & ! Procedure(s)
    convert_latent_ht_to_m_s 

  use time_dependent_input, only: &
    latent_ht_given, & ! Variable (s)
    sens_ht_given, &
    time_sfc_given, &
    CO2_sfc_given, &
    T_sfc_given, &
    time_select ! Procedure(s)

  use array_index, only: &
    sclr_idx_type

  implicit none

  intrinsic :: sqrt, max

  ! Parameter Constants
  real( kind = core_rknd ), parameter :: & 
    M_da  = 0.02897_core_rknd  ! Molecular weight of dry air.
  integer, parameter :: &
    ntimes = 49

  !--------------------- Input Variables ---------------------
  integer, intent(in) :: &
    sclr_dim, & 
    edsclr_dim

  type (sclr_idx_type), intent(in) :: &
    sclr_idx

  real(kind=time_precision), intent(in) ::  & 
    time      ! Current time                [s]

  real( kind = core_rknd ), intent(in) :: & 
    z,         & ! Elevation at zt=1           [m]
    rho_sfc,       & ! Air density at surface      [kg/m^3]
    thlm_sfc,  & ! Theta_l at zt(1)            [K]
    ubar         ! mean sfc wind speed         [m/s]

  !--------------------- Output Variables ---------------------
  real( kind = core_rknd ), intent(out) ::  & 
    wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
    wprtp_sfc,   & ! w'rt' surface flux        [(m kg)/(kg s)]
    ustar,       & ! surface friction velocity [m/s]
    T_sfc          ! Temperature at the surface [K]

  ! Output variables
  real( kind = core_rknd ), intent(out), dimension(sclr_dim) ::  & 
    wpsclrp_sfc    ! w'sclr' surface flux          [units m/s]

  real( kind = core_rknd ), intent(out), dimension(edsclr_dim) ::  & 
    wpedsclrp_sfc  ! w' edsclr' surface flux       [units m/s]

  !--------------------- Local Variables ---------------------
  integer :: &
    before_time, after_time

  real( kind = core_rknd ) ::  &  
    heat_flx, moisture_flx, &                 ! [W/m^2]
    heat_flx2, moisture_flx2, &
    time_frac, &
    bflx

  real( kind = core_rknd ) :: &
    CO2_flx, &
    CO2_flx2

  ! COBRA roughness height
  ! real, parameter :: z0 = 0.035_core_rknd  ! ARM momentum roughness height
  real( kind = core_rknd ), parameter :: z0 = 1.75_core_rknd   ! momentum roughness height

  !--------------------- Begin Code ---------------------

  ! Default Initialization
  heat_flx = 0.0_core_rknd
  moisture_flx = 0.0_core_rknd
  CO2_flx = 0.0_core_rknd
  CO2_flx2 = 0.0_core_rknd ! Default initialization

  ! Compute heat and moisture fluxes from ARM data in (W/m2)

  ! Use time_select to caluclate the indexes before and after time
  ! and the time fraction necessary for linear_interp_factor
  call time_select( time, ntimes, time_sfc_given, &
                    before_time, after_time, time_frac )

  ! Interpolate fluxes
  heat_flx = linear_interp_factor( time_frac, sens_ht_given(after_time), &
                                   sens_ht_given(before_time) )
  moisture_flx = linear_interp_factor( time_frac, latent_ht_given(after_time), &
                                       latent_ht_given(before_time) )
  CO2_flx = linear_interp_factor( time_frac, CO2_sfc_given(after_time), &
                                  CO2_sfc_given(before_time) )
  T_sfc = linear_interp_factor( time_frac, T_sfc_given(after_time), &
                                T_sfc_given(before_time) )

  ! Convert heat_flx and moisture_flx to natural units
  heat_flx2     = convert_sens_ht_to_km_s( heat_flx, rho_sfc )    ! (K m/s)
  moisture_flx2 = convert_latent_ht_to_m_s( moisture_flx, rho_sfc )! (m/s)

  !       Convert CO2 surface flux to natural units.
  !       The CO2 flux has been given in units of:  umol/(m^2 s).
  !       umol stands for micromoles.  The CO2 concentration in
  !       this code is in units of ppmv, which is also the molar
  !       mixing ratio times 10^6.
  !       The units are:  10^6 * [ mol (CO2) / mol (dry air) ].
  !       w'CO2' = (Flux) * [ M (dry air) / rho (dry air) ];
  !       where M is the molecular weight of dry air.
  CO2_flx2 = CO2_flx * ( M_da / rho_sfc )

  ! Heat flux in units of (m2/s3) (needed by diag_ustar)
  bflx = grav/thlm_sfc * heat_flx2

  ! Compute ustar
  ustar = diag_ustar( z, bflx, ubar, z0 )

  ! Assign fluxes

  wpthlp_sfc = heat_flx2
  wprtp_sfc  = moisture_flx2

  if ( sclr_idx%iisclr_CO2 > 0 ) wpsclrp_sfc(sclr_idx%iisclr_CO2) = CO2_flx2
  if ( sclr_idx%iisclr_thl > 0 ) wpsclrp_sfc(sclr_idx%iisclr_thl) = wpthlp_sfc
  if ( sclr_idx%iisclr_rt  > 0 ) wpsclrp_sfc(sclr_idx%iisclr_rt)  = wprtp_sfc

  if ( sclr_idx%iiedsclr_CO2 > 0 ) wpedsclrp_sfc(sclr_idx%iiedsclr_CO2) = CO2_flx2
  if ( sclr_idx%iiedsclr_thl > 0 ) wpedsclrp_sfc(sclr_idx%iiedsclr_thl) = wpthlp_sfc
  if ( sclr_idx%iiedsclr_rt  > 0 ) wpedsclrp_sfc(sclr_idx%iiedsclr_rt)  = wprtp_sfc

  return
  end subroutine cobra_sfclyr

end module cobra
