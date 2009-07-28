!----------------------------------------------------------------------
! $Id$
module arm_0003

!       Description:
!       Contains subroutines for the March 2000 IOP ARM case.
!----------------------------------------------------------------------

  implicit none

  public :: arm_0003_sfclyr

  private ! Default Scope


  contains

  !----------------------------------------------------------------------
  subroutine arm_0003_sfclyr( time, z, rho0, & 
                              thlm_sfc, um_sfc, vm_sfc,  & 
                              upwp_sfc, vpwp_sfc, & 
                              wpthlp_sfc, wprtp_sfc, ustar, & 
                              wpsclrp_sfc, wpedsclrp_sfc )
    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS ARM specifications
    !----------------------------------------------------------------------

    use constants, only: Cp, Lv, grav ! Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use diag_ustar_mod, only: diag_ustar ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    use interpolation, only: factor_interp ! Procedure(s)

    use surface_flux, only: compute_ubar, compute_momentum_flux ! Procedure(s)

    use time_dependent_input, only: LH_given, SH_given, time_sfc_given ! Variable(s)

    implicit none

    intrinsic :: max, sqrt

    real, parameter ::  &  
      z0    = 0.035   ! ARM Cu mom. roughness height

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
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
    real :: ubar, bflx, heat_flx, moisture_flx, time_frac
    integer :: i1, i2
    integer :: ntimes
    !----------------------------------------------------------------------

    ! Default Initialization
    heat_flx = 0.0
    moisture_flx = 0.0

    ntimes = size (time_sfc_given)

    if ( time <= time_sfc_given(1) ) then
      heat_flx     = SH_given(1)
      moisture_flx = LH_given(1)
    else if ( time >= time_sfc_given(ntimes) ) then
      heat_flx     = SH_given(ntimes)
      moisture_flx = LH_given(ntimes)
    else
      i1 = 1
      do while ( i1 <= ntimes-1 )
        i2 = i1 + 1
        if ( time >= time_sfc_given(i1) .and. time < time_sfc_given(i2) ) then
          time_frac       = real((time-time_sfc_given(i1))/(time_sfc_given(i2) & 
             - time_sfc_given(i1)))
          heat_flx = factor_interp( time_frac, SH_given(i2), SH_given(i1) )
          moisture_flx = factor_interp( time_frac, LH_given(i2), LH_given(i1) )
          i1           = ntimes
        end if
        i1 = i2
      end do
    end if ! time <= time_sfc_given(1)

    ! Convert W/m^2 into w'thl' w'rt' units
    wpthlp_sfc = heat_flx / ( Cp * rho0 )     ! (K m/s)
    wprtp_sfc  = moisture_flx / ( Lv * rho0 ) ! (kg m/ kg s)

    ! Compute momentum fluxes using ARM Cu formulae

    ubar = compute_ubar( um_sfc, vm_sfc )

    bflx = grav/thlm_sfc * wpthlp_sfc

    ! Compute ustar
    ustar = diag_ustar( z, bflx, ubar, z0 )

    call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                                upwp_sfc, vpwp_sfc )

    ! Let passive scalars be equal to rt and theta_l for testing
    if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
    if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

    if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
    if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc

    return
  end subroutine arm_0003_sfclyr

!----------------------------------------------------------------------
end module arm_0003
