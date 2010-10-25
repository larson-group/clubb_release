!----------------------------------------------------------------------
! $Id$
module lba

  ! Description:
  !   Contains subroutines for the LBA case.
  !----------------------------------------------------------------------

  implicit none

  private ! Default Scope

  public :: lba_tndcy, lba_sfclyr

  contains

  !----------------------------------------------------------------------
  subroutine lba_tndcy( thlm_forcing, rtm_forcing, & 
                        sclrm_forcing, edsclrm_forcing )
    !       Description:
    !       Subroutine to set theta and water tendencies for LBA case.

    !       References:
    !----------------------------------------------------------------------

    use grid_class, only: gr !  Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    use array_index, only:  & 
        iisclr_thl, iisclr_rt ! Variable(s)


    implicit none

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) :: & 
      thlm_forcing, & ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing     ! Total water mixing ratio tendency            [kg/kg/s]

    real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar forcing [units vary]

    real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
      edsclrm_forcing ! Passive eddy-scalar forcing [units vary]

    ! ---- Begin Code ----

    ! Large-scale temperature tendency
    thlm_forcing(:) = 0.0

    ! Large-scale advective moisture tendency
    rtm_forcing(:) = 0.0

    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

    return
  end subroutine lba_tndcy

  !----------------------------------------------------------------------
  subroutine lba_sfclyr( time, z, rho0, & 
                         thlm_sfc, ubar,  & 
                         wpthlp_sfc, wprtp_sfc, ustar )

    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS BOMEX specifications

    !       References:
    !       Grabowski, et al. (2005)
    !----------------------------------------------------------------------

    use constants_clubb, only: pi, grav, Lv, Cp ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    implicit none

    intrinsic :: max, sqrt

    ! Constant Parameters
    real, parameter ::  & 
      z0    = 0.035  ! ARM mom. roughness height

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      time      ! Current time        [s]

    real, intent(in) ::  & 
      z,         & ! Height at zt=2      [m] 
      rho0,      & ! Density at zm=1     [kg/m^3] 
      thlm_sfc,  & ! thlm at (2)         [m/s]
      ubar

    ! Output variables
    real, intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variables
    real :: ft, bflx

    ! Compute heat and moisture fluxes
    ! From Table A.1.
    ft = real( max( 0._time_precision,  & 
                 cos( 0.5 * pi * ( (5.25 - ( time/3600.)) / 5.25 ) ) & 
            ) )

    wpthlp_sfc =  ( 270. * ft**1.5 ) / ( rho0 * Cp )
    wprtp_sfc  =  ( 554. * ft**1.3 ) / ( rho0 * Lv )

    bflx = grav/thlm_sfc * wpthlp_sfc

    ! Compute ustar
    ustar = diag_ustar( z, bflx, ubar, z0 )

    return
  end subroutine lba_sfclyr


end module lba

