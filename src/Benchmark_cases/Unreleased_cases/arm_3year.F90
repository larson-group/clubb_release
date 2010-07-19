!----------------------------------------------------------------------
! $Id$
module arm_3year

!       Description:
!       Contains subroutines for the ARM 3 Year case.
!----------------------------------------------------------------------


  implicit none

  public :: arm_3year_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine arm_3year_sfclyr( time, z, rho0, & 
                               thlm_sfc, ubar,  & 
                               wpthlp_sfc, wprtp_sfc, ustar )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ARM specifications
    !----------------------------------------------------------------------

    use constants_clubb, only: Cp, Lv, grav ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    use surface_flux, only: compute_ht_mostr_flux ! Procedures

    use time_dependent_input, only: time_sfc_given

    implicit none

    intrinsic :: max, sqrt, present

    real, parameter ::  &  
      z0    = 0.035   ! ARM Cu mom. roughness height

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      time      ! Current time        [s]

    real, intent(in) ::  & 
      z,         & ! Height at zt=2      [s] 
      rho0,      & ! Density at zm=1     [kg/m^3] 
      ubar,      & ! mean sfc wind speed [m/s]
      thlm_sfc     ! thlm at (2)         [m/s]

    ! Output variables
    real, intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variables
    real :: bflx, heat_flx, moisture_flx
    
    !-------------BEGIN CODE--------------

    ! Compute heat and moisture fluxes from ARM data in (W/m2)
    call compute_ht_mostr_flux( time, size( time_sfc_given ), &
                                heat_flx, moisture_flx )

    ! Convert W/m^2 into w'thl' w'rt' units
    wpthlp_sfc = heat_flx / ( Cp * rho0 )     ! (K m/s)
    wprtp_sfc  = moisture_flx / ( Lv * rho0 ) ! (kg m/ kg s)

    ! Compute momentum fluxes using ARM Cu formulae

    bflx = grav/thlm_sfc * wpthlp_sfc

    ! Compute ustar
    ustar = diag_ustar( z, bflx, ubar, z0 )

    return
  end subroutine arm_3year_sfclyr

end module arm_3year
