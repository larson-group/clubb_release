!----------------------------------------------------------------------
!$Id$
module atex_long

  !       Description:
  !       Contains subroutines for the GCSS ATEX case.
  !----------------------------------------------------------------------

  use clubb_precision, only: core_rknd

  use grid_class, only: &
    grid ! Type
  
  implicit none

  public :: atex_long_tndcy, atex_long_sfclyr

  private ! Default Scope

  contains

  subroutine calc_forcings( ngrdcol, gr, thlm_forcing, rtm_forcing )

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol
    type (grid), intent(in) :: gr

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt) :: & 
      thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
      rtm_forcing     ! Total water mixing ratio tendency           [kg/kg/s]

    !--------------------- Local Variables ---------------------
    integer :: i, k

    !--------------------- Begin Code ---------------------
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
      do k = 1, gr%nzt

        ! Theta-l tendency ! Hing - known magic numbers
        if      ( gr%zt(i,k) >=    0._core_rknd .and. gr%zt(i,k) < 1400._core_rknd ) then
          thlm_forcing(i,k) = -3.5805e-5_core_rknd
        else if ( gr%zt(i,k) >= 1400._core_rknd .and. gr%zt(i,k) < 1650._core_rknd ) then
          thlm_forcing(i,k) = -3.5805e-5_core_rknd &
                              +1.1935e-5_core_rknd*(gr%zt(i,k)-1400._core_rknd)*0.004_core_rknd 
        else if ( gr%zt(i,k) >= 1650._core_rknd .and. gr%zt(i,k) < 2990._core_rknd ) then
          thlm_forcing(i,k) = -2.3870e-5_core_rknd &
                              -0.1155e-5_core_rknd*(gr%zt(i,k)-1650._core_rknd)/1350._core_rknd
        else
          thlm_forcing(i,k) =         0._core_rknd
        end if

        ! Moisture tendency ! Hing - known magic numbers
        if ( gr%zt(i,k) >= 0._core_rknd .and. gr%zt(i,k) < 1050._core_rknd ) then
          ! Brian - known magic number
          rtm_forcing(i,k) = -1.58e-8_core_rknd * ( 1._core_rknd - gr%zt(i,k)/1050._core_rknd )  
        else
          ! Brian
          rtm_forcing(i,k) = 0.0_core_rknd       
        end if
      end do
    end do

  end subroutine calc_forcings

  !======================================================================
  subroutine atex_long_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &
                              gr, time, &
                              wm_zt, wm_zm, & 
                              thlm_forcing, rtm_forcing, &
                              sclrm_forcing, edsclrm_forcing )
  ! Description:
  !   Subroutine to set theta-l and water tendencies and subsidence for the long ATEX case

  ! References:
  !   Ong, H. The nontraditional Coriolis terms and trade-wind cumuli

  !----------------------------------------------------------------------

  use grid_class, only: &
    grid ! Type

  use grid_class, only: &
    zt2zm_api ! Procedure(s)

  use clubb_precision, only: &
    time_precision, & ! Variable(s)
    core_rknd

  use array_index, only: &
    sclr_idx_type

  implicit none

  !--------------------- Input Variables ---------------------
  integer, intent(in) :: &
    ngrdcol, &
    sclr_dim, &
    edsclr_dim

  type (sclr_idx_type), intent(in) :: &
    sclr_idx

  type (grid), intent(in) :: gr

  real(kind=time_precision), intent(in) :: &
    time ! Current time [s]

  !--------------------- Output Variables ---------------------
  real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt) :: &
    wm_zt,        & ! w wind on thermodynamic grid                [m/s]
    thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
    rtm_forcing     ! Total water mixing ratio tendency           [kg/kg/s]

  real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzm) :: &
    wm_zm           ! w wind on momentum grid                     [m/s]

  real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt, sclr_dim) :: &
    sclrm_forcing   ! Passive scalar tendency         [units/s]

  real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt, edsclr_dim) :: &
    edsclrm_forcing ! Eddy-passive scalar tendency    [units/s]

  !--------------------- Local Variables ---------------------
  integer :: i, k

  !--------------------- Begin Code ---------------------
  ! Large scale subsidence ! Hing - known magic numbers
  !$acc parallel loop gang vector collapse(2) default(present)
  do i = 1, ngrdcol
    do k = 1, gr%nzt

      if      ( gr%zt(i,k) >=    0._core_rknd .and. gr%zt(i,k) < 1050._core_rknd ) then
        wm_zt(i,k) = -0.00636_core_rknd *  gr%zt(i,k)                  / 1050._core_rknd
      else if ( gr%zt(i,k) >= 1050._core_rknd .and. gr%zt(i,k) < 1650._core_rknd ) then
        wm_zt(i,k) = -0.00636_core_rknd &
                     -0.00079_core_rknd * (gr%zt(i,k)-1050._core_rknd) /  600._core_rknd
      else
        wm_zt(i,k) = -0.00715_core_rknd
      end if

    end do
  end do

  ! Spin up period ! Hing - known magic number
!  if ( time < 43200.0_time_precision ) then
!    wm_zt(:,:) = wm_zt(:,:) * time / 43200.0_time_precision
!  end if

  wm_zm = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, wm_zt )

  ! Boundary conditions.
  !$acc parallel loop gang vector default(present)
  do i = 1, ngrdcol
    wm_zm(i,1) = 0.0_core_rknd       ! At surface
    wm_zm(i,gr%nzm) = 0.0_core_rknd  ! Model top
  end do

  call calc_forcings( ngrdcol, gr,              &                ! intent(in)
                      thlm_forcing, rtm_forcing )                ! intent(out)

  ! Spin up period ! Hing - known magic number
!  if ( time < 43200.0_time_precision ) then
!    thlm_forcing(:,:) = thlm_forcing(:,:) * time / 43200.0_time_precision
!     rtm_forcing(:,:) =  rtm_forcing(:,:) * time / 43200.0_time_precision
!  end if

  if ( sclr_dim > 0 ) then
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzt
      do i = 1, ngrdcol
        ! Test scalars with thetal and rt if desired
        if ( sclr_idx%iisclr_thl > 0 ) sclrm_forcing(i,k,sclr_idx%iisclr_thl) = thlm_forcing(i,k)
        if ( sclr_idx%iisclr_rt  > 0 ) sclrm_forcing(i,k,sclr_idx%iisclr_rt)  = rtm_forcing(i,k)
      end do
    end do
  end if

  if ( edsclr_dim > 0 ) then
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzt
      do i = 1, ngrdcol
        if ( sclr_idx%iiedsclr_thl > 0 ) &
                edsclrm_forcing(i,k,sclr_idx%iiedsclr_thl) = thlm_forcing(i,k)
        if ( sclr_idx%iiedsclr_rt  > 0 ) &
                edsclrm_forcing(i,k,sclr_idx%iiedsclr_rt)  = rtm_forcing(i,k)
      end do
    end do
  end if

  return

  end subroutine atex_long_tndcy

  !======================================================================
  subroutine atex_long_sfclyr( ngrdcol, time, ubar, & 
                               thlm_sfc, rtm_sfc, exner_sfc, & 
                               wpthlp_sfc, wprtp_sfc, ustar, T_sfc )
  ! Description:
  !   This subroutine computes surface fluxes of
  !   heat and moisture according to GCSS ATEX specifications

  ! References:   
  !   B. Stevens et al., 2000: Simulations of trade-wind cumuli 
  !   under a strong inversion, J. Atmos. Sci, 58,1870-1891.
  !   http://www.atmos.washington.edu/~breth/GCSS/Stevens_etal_
  !   ATEX_JAS2001.pdf
  !----------------------------------------------------------------------

  use sfc_flux, only: compute_wpthlp_sfc, compute_wprtp_sfc

  use interpolation, only: linear_interp_factor ! Procedure(s)

  use time_dependent_input, only: time_sfc_given, T_sfc_given, & ! Variable(s)
                                  time_select                    ! Procedure(s)

  use clubb_precision, only: time_precision, core_rknd ! Variable(s)

  implicit none

  ! Input variables
  integer, intent(in) :: &
    ngrdcol

  real(time_precision), intent(in) :: &
    time       ! the current time [s]

  real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
    ubar,    & ! mean sfc wind speed                           [m/s]
    thlm_sfc,& ! theta_l at first model layer                  [K]
    rtm_sfc, & ! Total water mixing ratio at first model layer [kg/kg]
    exner_sfc  ! Exner function                                [-]

  ! Output variables
  real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
    wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
    wprtp_sfc,   &    ! w'rt' surface flux        [(m kg)/(kg s)]
    ustar,       &
    T_sfc          ! Surface temperature                           [K]
    
  ! Local Variable
  real( kind = core_rknd ) :: & 
    time_frac, & ! Time fraction used for interpolation
    T_sfc_interp

  real( kind = core_rknd ), dimension(ngrdcol) :: & 
    C_10, &      ! Coefficient
    adjustment   ! The adjustment for compute_wprtp_sfc

  integer :: &
    before_time, after_time, i ! The 

  !-----------------BEGIN CODE-----------------------

  !$acc enter data create( C_10, adjustment )

  ! Interpolate T_sfc from time_dependent_input

  call time_select( time, size(time_sfc_given), time_sfc_given, &
                    before_time, after_time, time_frac )

  T_sfc_interp = linear_interp_factor( time_frac, T_sfc_given(after_time), &
                                       T_sfc_given(before_time) )
                
  !$acc parallel loop gang vector default(present)
  do i = 1, ngrdcol
    C_10(i)       = 0.0013_core_rknd
    adjustment(i) = 0.0194664_core_rknd ! 0.981*sat_mixrat_liq_k( 101600.,298.,saturation_flatau )
                                        ! Hing
    ustar(i)      = 0.3_core_rknd
    T_sfc(i)      = T_sfc_interp
  end do

  ! Compute wpthlp_sfc and wprtp_sfc
  call compute_wpthlp_sfc( ngrdcol, C_10, ubar, thlm_sfc, T_sfc, exner_sfc, &
                           wpthlp_sfc ) 

  call compute_wprtp_sfc( ngrdcol, C_10, ubar, rtm_sfc, adjustment, &
                          wprtp_sfc )

  !$acc exit data delete( C_10, adjustment )

  return

  end subroutine atex_long_sfclyr

!----------------------------------------------------------------------
end module atex_long
