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
  subroutine atex_tndcy( time, time_initial, &
                         rtm, &
                         err_code, &
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

  use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

  use grid_class, only: gr ! Variable(s)

  use grid_class, only: zt2zm ! Procedure(s)

  use clubb_precision, only: time_precision ! Variable(s)

  use error_code, only: clubb_rtm_level_not_found ! Variable(s)

  use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)
   
  implicit none

  ! Input Variables
  real(kind=time_precision), intent(in) ::  & 
    time,         & ! Current time     [s]
    time_initial ! Initial time     [s]

  real, intent(in), dimension(gr%nzmax) :: & 
    rtm      ! Total water mixing ratio        [kg/kg]

  ! Input/output
  integer, intent(inout) :: err_code ! Diagnostic 

  ! Output Variables
  real, intent(out), dimension(gr%nzmax) :: & 
    wm_zt,        & ! w wind on thermodynamic grid                [m/s]
    wm_zm,        & ! w wind on momentum grid                     [m/s]
    thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
    rtm_forcing     ! Total water mixing ratio tendency           [kg/kg/s]


  real, intent(out), dimension(gr%nzmax, sclr_dim) :: & 
    sclrm_forcing   ! Passive scalar tendency         [units/s]

  real, intent(out), dimension(gr%nzmax, edsclr_dim) :: & 
    edsclrm_forcing ! Eddy-passive scalar tendency    [units/s]

  ! Internal variables
  integer :: i
  real :: z_inversion

  ! Forcings are applied only after t = 5400 s
  wm_zt = 0.
  wm_zm = 0.

  thlm_forcing = 0.
  rtm_forcing  = 0.

  if ( time >= time_initial + 5400.0_time_precision ) then

  !  Identify height of 6.5 g/kg moisture level

     i = 2
     do while ( i <= gr%nzmax .and. rtm(i) > 6.5e-3 )
        i = i + 1
     end do
     if ( i == gr%nzmax+1 .or. i == 2 ) then
       write(fstderr,*) "Identification of 6.5 g/kg level failed"
       write(fstderr,*) "Subroutine: atex_tndcy. File: atex.F"
       write(fstderr,*) "i = ", i
       write(fstderr,*) "rtm(i) = ",rtm(i)
       err_code = clubb_rtm_level_not_found
       return
     end if
     z_inversion = gr%zt(i-1)

  !          Large scale subsidence

     do i = 2, gr%nzmax

        if ( gr%zt(i) > 0. .and. gr%zt(i) <= z_inversion ) then
           wm_zt(i)  & 
             = -0.0065 * gr%zt(i)/z_inversion ! Known magic number
        else if ( gr%zt(i) > z_inversion .and. gr%zt(i) <= z_inversion+300. ) then
           wm_zt(i) & 
             = - 0.0065 * ( 1. - (gr%zt(i)-z_inversion)/300. ) ! Known magic number
        else
           wm_zt(i) = 0.
        end if

     end do

     wm_zm = zt2zm( wm_zt )

     ! Boundary conditions.
     wm_zt(1) = 0.0        ! Below surface
     wm_zm(1) = 0.0        ! At surface
     wm_zm(gr%nzmax) = 0.0  ! Model top

     ! Theta-l tendency

     do i = 2, gr%nzmax

        if ( gr%zt(i) > 0. .and. gr%zt(i) < z_inversion ) then
           thlm_forcing(i) = -1.1575e-5 * ( 3. - gr%zt(i)/z_inversion ) ! Known magic number
        else if ( gr%zt(i) > z_inversion .and. gr%zt(i) <= z_inversion+300. ) then
           thlm_forcing(i) = -2.315e-5 * ( 1. - (gr%zt(i)-z_inversion)/300. ) ! Known magic number
        else
           thlm_forcing(i) = 0.0
        end if

     end do

     ! Moisture tendency
     do i = 2, gr%nzmax

        if ( gr%zt(i) > 0. .and. gr%zt(i) < z_inversion ) then
           rtm_forcing(i) = -1.58e-8 * ( 1. - gr%zt(i)/z_inversion )  ! Brian - known magic number
        else
           rtm_forcing(i) = 0.0       ! Brian
        end if

     end do

     ! Boundary conditions
     thlm_forcing(1) = 0.0  ! Below surface
     rtm_forcing(1)  = 0.0  ! Below surface

  end if ! time >= time_initial + 5400.0

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

  use surface_flux, only: compute_wpthlp_sfc, compute_wprtp_sfc

  use interpolation, only: linear_interp_factor ! Procedure(s)

  use time_dependent_input, only: time_sfc_given, T_sfc_given, & ! Variable(s)
                                  time_select                    ! Procedure(s)

  use clubb_precision, only: time_precision ! Variable(s)

  implicit none

  ! Input variables
  real(time_precision), intent(in) :: &
    time       ! the current time [s]

  real, intent(in) ::  &
    ubar,    & ! mean sfc wind speed                           [m/s]
    thlm_sfc,& ! theta_l at first model layer                  [K]
    rtm_sfc, & ! Total water mixing ratio at first model layer [kg/kg]
    exner_sfc  ! Exner function                                [-]

  ! Output variables
  real, intent(out) ::  & 
    wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
    wprtp_sfc,   &    ! w'rt' surface flux        [(m kg)/(kg s)]
    ustar,       &
    T_sfc          ! Surface temperature                           [K]
    
  ! Local Variable
  real :: & 
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

  C_10 = 0.0013
  ustar = 0.3
  adjustment = 0.0198293

  wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, T_sfc, exner_sfc )
  wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, adjustment )

  return
  end subroutine atex_sfclyr

!----------------------------------------------------------------------
end module atex
