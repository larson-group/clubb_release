!----------------------------------------------------------------------
! $Id$
module lba

  !       Description:
  !       Contains subroutines for the LBA case.
  !----------------------------------------------------------------------

  implicit none

  private ! Default Scope

  public :: lba_tndcy, lba_sfclyr, lba_init

  private :: zrad, krad, ntimes, nzrad

  ! Constant Parameters
  integer, parameter :: ntimes = 36, nzrad = 33

  real, dimension(nzrad) :: & 
    zrad

  real, dimension(nzrad, ntimes) :: & 
    krad

  integer, parameter :: per_line = 5


  contains

  !----------------------------------------------------------------------
  subroutine lba_tndcy( time, & 
                        thlm_forcing, rtm_forcing, & 
                        sclrm_forcing, edsclrm_forcing )
    !       Description:
    !       Subroutine to set theta and water tendencies for LBA case.

    !       References:
    !----------------------------------------------------------------------

    use grid_class, only: gr !  Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    use interpolation, only: zlinterp_fnc ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    use array_index, only:  & 
        iisclr_thl, iisclr_rt ! Variable(s)

    use interpolation, only: factor_interp ! Procedure(s)

    use parameters_radiation, only: rad_scheme ! Variable(s)

    implicit none

    ! Input
    real(kind=time_precision), intent(in) :: time ! Model time [s]

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) :: & 
      thlm_forcing, & ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing     ! Total water mixing ratio tendency            [kg/kg/s]

    real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar forcing [units vary]

    real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
      edsclrm_forcing ! Passive eddy-scalar forcing [units vary]

    ! Local Variables
    real, dimension(gr%nnzp) :: radht
    real, dimension(nzrad) :: radhtz
    real :: a
    integer :: i1, i2


    if ( trim( rad_scheme ) == "lba" ) then

      ! Calculate radiative heating rate
      if ( time <=  600. ) then
        radhtz = krad(:,1)

      else if ( time >= ntimes * 600. ) then
        radhtz = krad(:,ntimes)

      else
        i1 = 1
        do while ( i1 <= ntimes-1 )
          i2 = i1 + 1
          if ( time >= 600. * i1 .and. time < 600. * i2  ) then
            a  = real(( time - 600. * i1 )/( 600. * i2 - 600. * i1))
            radhtz(:) = factor_interp( a, krad(:,i2), krad(:,i1) )
            i1     = ntimes
          end if
          i1 = i2
        end do
      end if ! time <= times(1)

      radht = zlinterp_fnc( gr%nnzp, nzrad, gr%zt, zrad, radhtz )

      ! Radiative theta-l tendency

      thlm_forcing = radht

    else ! Compute heating rate interactively with BUGSrad

      thlm_forcing = 0.0

    end if ! Simplified radiation

    ! Boundary conditions
    thlm_forcing(1) = 0.0  ! Below surface

    ! Large scale advective moisture tendency
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

  !----------------------------------------------------------------
  subroutine lba_init( iunit, file_path )
    !
    !       Description:
    !       This subroutine initializes the module by reading in forcing
    !       data used in the tndcy subroutine.
    !----------------------------------------------------------------

    use file_functions, only: file_read_1d, file_read_2d ! Procedure(s)

    implicit none

    integer, intent(in) :: iunit ! File unit number

    character(len=*), intent(in) :: &
      file_path ! Path to the forcing files

    call file_read_1d( iunit, & 
      file_path//'lba_heights.dat', & 
      nzrad, per_line, zrad )

    call file_read_2d( iunit, & 
      file_path//'lba_rad.dat', & 
      nzrad, ntimes, per_line, krad )

    return
  end subroutine lba_init

end module lba

