!----------------------------------------------------------------------
!$Id$
module atex

  !       Description:
  !       Contains subroutines for the GCSS ATEX case.
  !----------------------------------------------------------------------

  implicit none

  public :: atex_tndcy, atex_sfclyr

  private ! Default Scope

  contains

  !======================================================================
  subroutine atex_tndcy( gr, time, time_initial, &
                         rtm, &
                         wm_zt, wm_zm, & 
                         thlm_forcing, rtm_forcing, & 
                         sclrm_forcing, edsclrm_forcing )
  ! Description:
  !   Subroutine to set theta-l and water tendencies for ATEX case

  ! References:
  !   B. Stevens et al., 2000: Simulations of trade-wind cumuli 
  !   under a strong inversion, J. Atmos. Sci, 58,1870-1891.
  !   http://www.atmos.washington.edu/~breth/GCSS/Stevens_etal_
  !   ATEX_JAS2001.pdf

  !----------------------------------------------------------------------

  use constants_clubb, only: fstderr ! Constant(s)

  use grid_class, only: grid ! Type

  use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

  use grid_class, only: zt2zm ! Procedure(s)

  use clubb_precision, only: time_precision, core_rknd ! Variable(s)

  use error_code, only: &
        clubb_fatal_error, &            ! Constant
        err_code                        ! Error indicator

  use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)
   
  implicit none

  type (grid), target, intent(in) :: gr

  ! Input Variables
  real(kind=time_precision), intent(in) ::  & 
    time,         & ! Current time     [s]
    time_initial ! Initial time     [s]

  real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
    rtm      ! Total water mixing ratio        [kg/kg]

  ! Output Variables
  real( kind = core_rknd ), intent(out), dimension(gr%nz) :: & 
    wm_zt,        & ! w wind on thermodynamic grid                [m/s]
    wm_zm,        & ! w wind on momentum grid                     [m/s]
    thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
    rtm_forcing     ! Total water mixing ratio tendency           [kg/kg/s]


  real( kind = core_rknd ), intent(out), dimension(gr%nz, sclr_dim) :: & 
    sclrm_forcing   ! Passive scalar tendency         [units/s]

  real( kind = core_rknd ), intent(out), dimension(gr%nz, edsclr_dim) :: & 
    edsclrm_forcing ! Eddy-passive scalar tendency    [units/s]

  ! Internal variables
  integer :: i
  real( kind = core_rknd ) :: z_inversion

  ! Forcings are applied only after t = 5400 s
  wm_zt = 0._core_rknd
  wm_zm = 0._core_rknd


  thlm_forcing = 0._core_rknd

  rtm_forcing  = 0._core_rknd

  if ( time >= time_initial + 5400.0_time_precision ) then

  !  Identify height of 6.5 g/kg moisture level

     i = 2
     do while ( i <= gr%nz .and. rtm(i) > 6.5e-3_core_rknd )
        i = i + 1
     end do
     if ( i == gr%nz+1 .or. i == 2 ) then
       write(fstderr,*) "Identification of 6.5 g/kg level failed"
       write(fstderr,*) "Subroutine: atex_tndcy. File: atex.F"
       write(fstderr,*) "i = ", i
       write(fstderr,*) "rtm(i) = ",rtm(i)
       err_code = clubb_fatal_error
       return
     end if
     z_inversion = gr%zt(i-1)

  !          Large scale subsidence

     do i = 2, gr%nz

        if ( gr%zt(i) > 0._core_rknd .and. gr%zt(i) <= z_inversion ) then
           wm_zt(i)  & 
             = -0.0065_core_rknd * gr%zt(i)/z_inversion ! Known magic number
        else if ( gr%zt(i) > z_inversion .and. gr%zt(i) <= z_inversion+300._core_rknd ) then
           wm_zt(i) & 
             = - 0.0065_core_rknd * ( 1._core_rknd - (gr%zt(i)-z_inversion)/&
                 300._core_rknd ) ! Known magic number
        else
           wm_zt(i) = 0._core_rknd
        end if

     end do

     wm_zm = zt2zm( gr, wm_zt )

     ! Boundary conditions.
     wm_zt(1) = 0.0_core_rknd        ! Below surface
     wm_zm(1) = 0.0_core_rknd        ! At surface
     wm_zm(gr%nz) = 0.0_core_rknd  ! Model top

     ! Theta-l tendency

     do i = 2, gr%nz

        if ( gr%zt(i) > 0._core_rknd .and. gr%zt(i) < z_inversion ) then
           thlm_forcing(i) = -1.1575e-5_core_rknd * ( 3._core_rknd - &
             gr%zt(i)/z_inversion ) ! Known magic number
        else if ( gr%zt(i) > z_inversion .and. gr%zt(i) <= z_inversion+300._core_rknd ) then
           thlm_forcing(i) = -2.315e-5_core_rknd * ( 1._core_rknd - &
             (gr%zt(i)-z_inversion)/300._core_rknd ) ! Known magic number
        else
           thlm_forcing(i) = 0.0_core_rknd
        end if

     end do

     ! Moisture tendency
     do i = 2, gr%nz

        if ( gr%zt(i) > 0._core_rknd .and. gr%zt(i) < z_inversion ) then
           rtm_forcing(i) = -1.58e-8_core_rknd * ( 1._core_rknd - &
             gr%zt(i)/z_inversion )  ! Brian - known magic number
        else
           rtm_forcing(i) = 0.0_core_rknd       ! Brian
        end if

     end do

     ! Boundary conditions
     thlm_forcing(1) = 0.0_core_rknd  ! Below surface
     rtm_forcing(1)  = 0.0_core_rknd  ! Below surface

  end if ! time >= time_initial + 5400.0_core_rknd

  ! Test scalars with thetal and rt if desired
  if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
  if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

  if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
  if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

  return
  end subroutine atex_tndcy

  !======================================================================
  subroutine atex_sfclyr( time, ubar, & 
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
  real(time_precision), intent(in) :: &
    time       ! the current time [s]

  real( kind = core_rknd ), intent(in) ::  &
    ubar,    & ! mean sfc wind speed                           [m/s]
    thlm_sfc,& ! theta_l at first model layer                  [K]
    rtm_sfc, & ! Total water mixing ratio at first model layer [kg/kg]
    exner_sfc  ! Exner function                                [-]

  ! Output variables
  real( kind = core_rknd ), intent(out) ::  & 
    wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
    wprtp_sfc,   &    ! w'rt' surface flux        [(m kg)/(kg s)]
    ustar,       &
    T_sfc          ! Surface temperature                           [K]
    
  ! Local Variable
  real( kind = core_rknd ) :: & 
    C_10, &      ! Coefficient
    time_frac, & ! Time fraction used for interpolation
    adjustment   ! The adjustment for compute_wprtp_sfc

  integer :: &
    before_time, after_time ! The 

  !-----------------BEGIN CODE-----------------------

  ! Interpolate T_sfc from time_dependent_input

  call time_select( time, size(time_sfc_given), time_sfc_given, &
                    before_time, after_time, time_frac )

  T_sfc = linear_interp_factor( time_frac, T_sfc_given(after_time), &
                                    T_sfc_given(before_time) )

  ! Compute wpthlp_sfc and wprtp_sfc

  C_10 = 0.0013_core_rknd
  ustar = 0.3_core_rknd
  adjustment = 0.0198293_core_rknd

  wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, T_sfc, exner_sfc )
  wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, adjustment )

  return
  end subroutine atex_sfclyr

!----------------------------------------------------------------------
end module atex
