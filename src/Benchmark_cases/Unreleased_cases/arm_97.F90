!----------------------------------------------------------------------
! $Id$
module arm_97

  !       Description:
  !       Contains subroutines for the July 26-30 1997 ARM IOP A case.
  !       
  !       References:
  !       http://www.mmm.ucar.edu/gcss-wg4/gcss/case3.html
  !----------------------------------------------------------------------

  implicit none

  public :: arm_97_sfclyr

  private ! Defualt Scope


  contains

  !----------------------------------------------------------------------
  subroutine arm_97_sfclyr( time, z, rho_sfc, & 
                            thlm_sfc, ubar,  & 
                            wpthlp_sfc, wprtp_sfc, ustar, T_sfc )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ARM specifications
    !
    !       References:
    !       http://www.mmm.ucar.edu/gcss-wg4/gcss/case3.html
    !----------------------------------------------------------------------

    use constants_clubb, only: grav ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    use interpolation, only: factor_interp ! Procedure(s)

    use error_code, only: clubb_debug ! Procedure(s)

    use surface_flux, only: convert_SH_to_km_s, convert_LH_to_m_s ! Procedure(s)

    use time_dependent_input, only: &
      time_select, &
      time_sfc_given, &
      LH_given, &
      SH_given, &
      T_sfc_given, &
      l_t_dependent

    implicit none

    intrinsic :: max, sqrt, present

    real, parameter ::  & 
      z0    = 0.035   ! ARM Cu mom. roughness height

    ! Input Variables
    real(time_precision), intent(in) ::  & 
      time      ! Current time        [s]

    real, intent(in) ::  & 
      z,         & ! Height at zt=2      [s] 
      rho_sfc,      & ! Density at zm=1     [kg/m^3] 
      ubar,      & ! mean sfc wind speed [m/s]
      thlm_sfc     ! thlm at (2)         [m/s]

    ! Output variables
    real, intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar,        & ! surface friction velocity [m/s]
      T_sfc           ! surface temperature [K]

    ! Local variables
    real :: bflx, heat_flx, moisture_flx, time_frac
    integer :: before_time, after_time
    !----------------------------------------------------------------------
    if( l_t_dependent ) then
      ! Default initialization
      heat_flx = 0.0
      moisture_flx = 0.0

      time_frac = -1.0 ! Default initialization

      call time_select( time, size(time_sfc_given), time_sfc_given, &
                                   before_time, after_time, time_frac )

      if( time_frac == -1.0 ) then
        call clubb_debug(1,"times is not sorted in arm_97_tndcy")
      endif


      heat_flx = factor_interp( time_frac, SH_given(after_time), SH_given(before_time) )
      moisture_flx = factor_interp( time_frac, LH_given(after_time), LH_given(before_time) )
      T_sfc = factor_interp( time_frac, T_sfc_given(after_time), T_sfc_given(before_time) )

      ! Convert W/m^2 into w'thl' w'rt' units
      wpthlp_sfc = convert_SH_to_km_s( heat_flx, rho_sfc )     ! (K m/s)
      wprtp_sfc  = convert_LH_to_m_s( moisture_flx, rho_sfc ) ! (kg m/ kg s)

      ! Compute momentum fluxes using ARM Cu formulae

      bflx = grav/thlm_sfc * wpthlp_sfc

      ! Compute ustar
      ustar = diag_ustar( z, bflx, ubar, z0 )

    endif
    return
  end subroutine arm_97_sfclyr
  !----------------------------------------------------------------------
end module arm_97
