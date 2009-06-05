!----------------------------------------------------------------------
! $Id$
module arm_97

  !       Description:
  !       Contains subroutines for the July 26-30 1997 ARM IOP A case.
  !----------------------------------------------------------------------

  implicit none

  public :: arm_97_tndcy, arm_97_sfclyr

  private ! Defualt Scope


  contains

  !----------------------------------------------------------------------
  subroutine arm_97_tndcy( time, &
                           thlm_forcing, &
                           rtm_forcing, um_hoc_grid, vm_hoc_grid, &
                           sclrm_forcing, edsclrm_forcing )
    !       Description:
    !       Subroutine to set thetal and total water tendencies for ARM 97 case

    !       References:
    !       None
    !----------------------------------------------------------------------

    use grid_class, only: gr ! Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use error_code, only: clubb_debug ! Procedure(s)

    use array_index, only:  & 
        iisclr_thl, iisclr_rt, iiedsclr_thl, iiedsclr_rt ! Variable(s)

    use interpolation, only: factor_interp ! Procedure(s)

    use time_dependant_input, only: &
      time_select,  &
      time_f_given, &
      thlm_f_given, &
      rtm_f_given,  &
      um_given,     &
      vm_given, &
      l_time_dependant


    implicit none

    ! External
    intrinsic :: present

    ! Input Variables
    real(kind=time_precision), intent(in) :: time ! Model time [s]

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) ::  & 
      thlm_forcing,  & ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing,   & ! Total water mixing ratio tendency            [kg/kg/s]
      um_hoc_grid,   & ! Observed wind, for nudging                   [m/s]
      vm_hoc_grid      ! Observed wind, for nudging                   [m/s]

    real, intent(out), dimension(gr%nnzp,sclr_dim) ::  & 
      sclrm_forcing ! Passive scalar tendency [units vary]

    real, intent(out), dimension(gr%nnzp,edsclr_dim) ::  & 
      edsclrm_forcing ! Eddy passive scalar tendency [units vary]

    real :: time_frac
    integer :: i1, i2
    if(l_time_dependant ) then
    !-----------------------------------------------------------------------

    ! Thetal forcing is equal to the LS tendency given here and the
    ! interactive calculation done in BUGSrad

    ! Compute time_frac, where 0 < time_frac < 1.
    !     time_frac measures the fractional time interval
    !     between the current time (time), and the previous
    !     time at which forcing data are available (times(i1))
    !     and the next time (times(i2)).

    time_frac = -1.0 ! Default initialization

    call time_select( time,size(time_f_given), time_f_given, i1, i2 )

    time_frac = real((time-time_f_given(i1))/(time_f_given(i2)-time_f_given(i1)))

    if( time_frac == -1.0 ) then
      call clubb_debug(1,"times is not sorted in arm_97_tndcy")
    endif

    ! Interpolate LS thetal tendency to the CLUBB grid
    ! Time
    thlm_forcing = factor_interp( time_frac, thlm_f_given(:,i2), thlm_f_given(  :,i1 ) )
    ! Interpolate LS rt tendency to the HOC grid

    ! Time
    rtm_forcing = factor_interp( time_frac, rtm_f_given(:,i2), rtm_f_given(:,i1) )
    
    ! Interpolate um observed to the HOC grid
    ! Time
    um_hoc_grid = factor_interp( time_frac, um_given(:,i2), um_given(:,i1) )

    ! Interpolate vm observed to the HOC grid
    vm_hoc_grid = factor_interp( time_frac, vm_given(:,i2), vm_given(:,i1) )

    ! Added by Joshua Fasching (October 27 2007)
    um_hoc_grid (1) = um_hoc_grid(2)
    vm_hoc_grid (1) = vm_hoc_grid(2)


    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing
    endif
    return
  end subroutine arm_97_tndcy
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

    use constants, only: Cp, Lv, grav ! Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use diag_ustar_mod, only: diag_ustar ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    use interpolation, only: factor_interp ! Procedure(s)

    use surface_flux, only: compute_ubar, compute_momentum_flux ! Procedure(s)
    use error_code, only: clubb_debug ! Procedure(s)

    use time_dependant_input, only: &
      time_select, &
      time_sfc_given, &
      LH_given, &
      SH_given, &
      l_time_dependant

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
    if( l_time_dependant ) then
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
