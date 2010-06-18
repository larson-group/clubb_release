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
                            thlm_sfc, um_sfc, vm_sfc,  & 
                            upwp_sfc, vpwp_sfc, & 
                            wpthlp_sfc, wprtp_sfc, ustar, & 
                            wpsclrp_sfc, wpedsclrp_sfc )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ARM specifications
    !----------------------------------------------------------------------

    use constants_clubb, only: Cp, Lv, grav ! Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    use interpolation, only: factor_interp ! Procedure(s)

    use surface_flux, only: compute_ubar, compute_momentum_flux ! Procedure(s)
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
      um_sfc,    & ! um at (2)           [m/s]
      vm_sfc,    & ! vm at (2)           [m/s]
      thlm_sfc     ! thlm at (2)         [m/s]

    ! Output variables
    real, intent(out) ::  & 
      upwp_sfc,     & ! u'w' at (1)      [m^2/s^2]
      vpwp_sfc,     & ! v'w'at (1)       [m^2/s^2]
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    real, intent(out), dimension(sclr_dim) ::  & 
      wpsclrp_sfc      ! Passive scalar surface flux      [units m/s]

    real, intent(out), dimension(edsclr_dim) ::  & 
      wpedsclrp_sfc    ! Passive eddy-scalar surface flux [units m/s]

    ! Local variables
    !        real :: ubar, ustar, bflx, heat_flx, moisture_flx, time_frac
    real :: ubar, bflx, heat_flx, moisture_flx, time_frac
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

      ! Let passive scalars be equal to rt and theta_l for now
      if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
      if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

      if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
      if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc

      ! Compute momentum fluxes using ARM Cu formulae

      ubar = compute_ubar( um_sfc, vm_sfc )

      bflx = grav/thlm_sfc * wpthlp_sfc

      ! Compute ustar
      ustar = diag_ustar( z, bflx, ubar, z0 )

      call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                                  upwp_sfc, vpwp_sfc )
    endif
    return
  end subroutine arm_97_sfclyr
  !----------------------------------------------------------------------
end module arm_97
