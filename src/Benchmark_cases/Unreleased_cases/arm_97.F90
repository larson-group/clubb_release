!----------------------------------------------------------------------
! $Id$
module arm_97

  !       Description:
  !       Contains subroutines for the July 26-30 1997 ARM IOP A case.
  !----------------------------------------------------------------------

  implicit none

  public :: arm_97_sfclyr

  private ! Defualt Scope


  contains

  !----------------------------------------------------------------------
  subroutine arm_97_sfclyr( time, z, rho0, & 
                            thlm_sfc, ubar,  & 
                            wpthlp_sfc, wprtp_sfc, ustar )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ARM specifications
    !----------------------------------------------------------------------

    use constants_clubb, only: Cp, Lv, grav ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    use interpolation, only: factor_interp ! Procedure(s)

    use error_code, only: clubb_debug ! Procedure(s)

    use time_dependent_input, only: &
      time_select, &
      time_sfc_given, &
      LH_given, &
      SH_given, &
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
      rho0,      & ! Density at zm=1     [kg/m^3] 
      ubar,      & ! mean sfc wind speed [m/s]
      thlm_sfc     ! thlm at (2)         [m/s]

    ! Output variables
    real, intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variables
    real :: bflx, heat_flx, moisture_flx, time_frac
    integer :: i1, i2
    !----------------------------------------------------------------------
    if( l_t_dependent ) then
      ! Default initialization
      heat_flx = 0.0
      moisture_flx = 0.0

      time_frac = -1.0 ! Default initialization

      call time_select( time, size(time_sfc_given), time_sfc_given, i1, i2 )

      time_frac = real((time-time_sfc_given(i1))/(time_sfc_given(i2)-time_sfc_given(i1)))

      if( time_frac == -1.0 ) then
        call clubb_debug(1,"times is not sorted in arm_97_tndcy")
      endif


      heat_flx = factor_interp( time_frac, SH_given(i2), SH_given(i1) )
      moisture_flx = factor_interp( time_frac, LH_given(i2), LH_given(i1) )

      ! Convert W/m^2 into w'thl' w'rt' units
      wpthlp_sfc = heat_flx / ( Cp * rho0 )     ! (K m/s)
      wprtp_sfc  = moisture_flx / ( Lv * rho0 ) ! (kg m/ kg s)

      ! Compute momentum fluxes using ARM Cu formulae

      bflx = grav/thlm_sfc * wpthlp_sfc

      ! Compute ustar
      ustar = diag_ustar( z, bflx, ubar, z0 )

    endif
    return
  end subroutine arm_97_sfclyr
  !----------------------------------------------------------------------
end module arm_97
