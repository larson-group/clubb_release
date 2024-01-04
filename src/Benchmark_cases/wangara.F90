!----------------------------------------------------------------------
! $Id$
module wangara

! References:
!   See below
!----------------------------------------------------------------------
  implicit none

  public :: wangara_tndcy, wangara_sfclyr

  private ! Default Scope

  contains
!----------------------------------------------------------------------
  subroutine wangara_tndcy( gr, sclr_dim, edsclr_dim, sclr_idx, &
                            wm_zt, wm_zm,  & 
                            thlm_forcing, rtm_forcing, & 
                            sclrm_forcing, edsclrm_forcing )
! Description:
!   Subroutine to set theta and water tendencies for Wangara case
! References:
!   ``A PDF-Based Model for Boundary Layer Clouds. Part II:
!   Model results'' Golaz, et al. (2002)
!   JAS, Vol. 59, pp. 3552--3571.
!----------------------------------------------------------------------

    use array_index, only: &
      sclr_idx_type

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use grid_class, only: &
      grid

    implicit none

    ! Input Variables
    type(grid), target, intent(in) :: &
      gr

    integer, intent(in) :: &
      sclr_dim, & 
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(gr%nz) :: & 
      wm_zt,        & ! w wind on thermodynamic grid                [m/s]
      wm_zm,        & ! w wind on momentum grid                     [m/s]
      thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
      rtm_forcing     ! Total water mixing ratio tendency           [kg/kg/s]

    real( kind = core_rknd ), intent(out), dimension(gr%nz,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar tendency [units/s]

    real( kind = core_rknd ), intent(out), dimension(gr%nz,edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar tendency [units/s]

    ! No large-scale subsidence for now
    wm_zt = 0.0_core_rknd
    wm_zm = 0.0_core_rknd

    ! No large-scale water tendency or cooling

    rtm_forcing  = 0.0_core_rknd
    thlm_forcing = 0.0_core_rknd

    ! Test scalars with thetal and rt if desired
    if ( sclr_idx%iisclr_thl > 0 ) sclrm_forcing(:,sclr_idx%iisclr_thl) = thlm_forcing
    if ( sclr_idx%iisclr_rt  > 0 ) sclrm_forcing(:,sclr_idx%iisclr_rt)  = rtm_forcing

    if ( sclr_idx%iiedsclr_thl > 0 ) edsclrm_forcing(:,sclr_idx%iiedsclr_thl) = thlm_forcing
    if ( sclr_idx%iiedsclr_rt  > 0 ) edsclrm_forcing(:,sclr_idx%iiedsclr_rt)  = rtm_forcing

    return
  end subroutine wangara_tndcy

!----------------------------------------------------------------------
  subroutine wangara_sfclyr( time, & 
                             wpthlp_sfc, wprtp_sfc, ustar )
! Description:
!   This subroutine computes surface fluxes of horizontal momentum,
!   heat and moisture for Wangara day 33
!
! References:
!   ``A PDF-Based Model for Boundary Layer Clouds. Part II:
!   Model results'' Golaz, et al. (2002)
!   JAS, Vol. 59, pp. 3552--3571.
!----------------------------------------------------------------------

    use constants_clubb, only: pi, fstderr, sec_per_day ! Variable(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    implicit none

    intrinsic :: mod, max, cos, sqrt, present

    ! Input variables
    real(kind=time_precision), intent(in) ::  & 
      time    ! Current time  [s]

    ! Output variables
    real( kind = core_rknd ), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variables
    real(kind=time_precision) :: &
      time_utc, time_est, &
      est_offset ! The offset for Australia EST time [s]



!---------------------BEGIN CODE-------------------------


    ! Declare the value of ustar.
    ustar = 0.13_core_rknd

    ! Compute UTC time of the day in seconds

    time_utc = mod( time, real( sec_per_day, kind=time_precision ) )

    ! Now convert UTC time to Australia EST (local time)
    est_offset = 36000._time_precision
    time_est = mod( time_utc + est_offset, real( sec_per_day, kind=time_precision ) )

    if ( time_est < 27000._time_precision .or. time_est > 63000._time_precision ) then
      write(fstderr,*) "wangara_sfclyr: error local time must" & 
        //" be between 730 and 1730._core_rknd"
      write(fstderr,*) 'time_est = ',time_est
      error stop
    end if

    ! Compute heat and moisture fluxes

    wpthlp_sfc = real(0.18_core_rknd * cos( (real( time_est, kind = core_rknd )-&
         45000.0_core_rknd)/36000.0_core_rknd * pi ), kind = core_rknd) ! Known magic number
    wprtp_sfc  = 1.3e-4_core_rknd * wpthlp_sfc ! Known magic number

    return
  end subroutine wangara_sfclyr

!----------------------------------------------------------------------

end module wangara
