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
                         rtm, rho, rcm, exner, &
                         err_code, &
                         wm_zt, wm_zm, Frad, radht, & 
                         thlm_forcing, rtm_forcing, & 
                         sclrm_forcing, edsclrm_forcing )
  !       Description:
  !       Subroutine to set theta-l and water tendencies for ATEX case

  !       References:

  !----------------------------------------------------------------------

  use constants_clubb, only: fstderr ! Constant(s)

  use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

  use parameters_radiation, only: rad_scheme ! Variable(s)

  use grid_class, only: gr ! Variable(s)

  use grid_class, only: zt2zm ! Procedure(s)

  use cloud_rad_module, only: cloud_rad ! Procedure(s)

  use stats_precision, only: time_precision ! Variable(s)

  use error_code, only: clubb_rtm_level_not_found ! Variable(s)

  use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)
   
  use stats_type, only: stat_update_var ! Procedure(s)

  use stats_variables, only: iradht_LW, zt, l_stats_samp ! Variable(s)

  implicit none

  ! Input Variables
  real(kind=time_precision), intent(in) ::  & 
    time,         & ! Current time     [s]
    time_initial ! Initial time     [s]

  real, intent(in), dimension(gr%nnzp) :: & 
    rtm,   & ! Total water mixing ratio        [kg/kg]
    rho,   & ! Density                         [kg/m^3]
    rcm,   & ! Liquid water mixing ratio       [kg/kg]
    exner    ! Exner function                  [-]

  ! Input/output
  integer, intent(inout) :: err_code ! Diagnostic 

  ! Output Variables
  real, intent(out), dimension(gr%nnzp) :: & 
    wm_zt,        & ! w wind on thermodynamic grid                [m/s]
    wm_zm,        & ! w wind on momentum grid                     [m/s]
    Frad,         & ! Radiative flux                              [W/m^2]
    radht,        & ! Radiative heating rate                      [K/s]
    thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
    rtm_forcing     ! Total water mixing ratio tendency           [kg/kg/s]


  real, intent(out), dimension(gr%nnzp, sclr_dim) :: & 
    sclrm_forcing   ! Passive scalar tendency         [units/s]

  real, intent(out), dimension(gr%nnzp, edsclr_dim) :: & 
    edsclrm_forcing ! Eddy-passive scalar tendency    [units/s]

  ! Internal variables
  integer :: i
  real :: zi

  ! Forcings are applied only after t = 5400 s
  wm_zt = 0.
  wm_zm = 0.

  thlm_forcing = 0.
  rtm_forcing  = 0.

  if ( time >= time_initial + 5400.0 ) then

  !  Identify height of 6.5 g/kg moisture level

     i = 2
     do while ( i <= gr%nnzp .and. rtm(i) > 6.5e-3 )
        i = i + 1
     end do
     if ( i == gr%nnzp+1 .or. i == 2 ) then
       write(fstderr,*) "Identification of 6.5 g/kg level failed"
       write(fstderr,*) "Subroutine: atex_tndcy. File: atex.F"
       write(fstderr,*) "i = ", i
       write(fstderr,*) "rtm(i) = ",rtm(i)
       err_code = clubb_rtm_level_not_found
       return
     end if
     zi = gr%zt(i-1)

  !          Large scale subsidence

     do i = 2, gr%nnzp

        if ( gr%zt(i) > 0. .and. gr%zt(i) <= zi ) then
           wm_zt(i)  & 
             = -0.0065 * gr%zt(i)/zi
        else if ( gr%zt(i) > zi .and. gr%zt(i) <= zi+300. ) then
           wm_zt(i) & 
             = - 0.0065 * ( 1. - (gr%zt(i)-zi)/300. )
        else
           wm_zt(i) = 0.
        end if

     end do

     wm_zm = zt2zm( wm_zt )

     ! Boundary conditions.
     wm_zt(1) = 0.0        ! Below surface
     wm_zm(1) = 0.0        ! At surface
     wm_zm(gr%nnzp) = 0.0  ! Model top

     ! Theta-l tendency

     do i = 2, gr%nnzp

        if ( gr%zt(i) > 0. .and. gr%zt(i) < zi ) then
           thlm_forcing(i) = -1.1575e-5 * ( 3. - gr%zt(i)/zi )
        else if ( gr%zt(i) > zi .and. gr%zt(i) <= zi+300. ) then
           thlm_forcing(i) = -2.315e-5 * ( 1. - (gr%zt(i)-zi)/300. )
        else
           thlm_forcing(i) = 0.0
        end if

     end do

     ! Moisture tendency
     do i = 2, gr%nnzp

        if ( gr%zt(i) > 0. .and. gr%zt(i) < zi ) then
           rtm_forcing(i) = -1.58e-8 * ( 1. - gr%zt(i)/zi )  ! Brian
        else
           rtm_forcing(i) = 0.0       ! Brian
        end if

     end do

     ! Boundary conditions
     thlm_forcing(1) = 0.0  ! Below surface
     rtm_forcing(1)  = 0.0  ! Below surface

  end if ! time >= time_initial + 5400.0

  ! Use cloud_rad() to compute radiation
  if ( trim( rad_scheme ) == "simplified" ) then

    call cloud_rad( rho, rcm, exner, Frad, radht, thlm_forcing )

    if ( l_stats_samp ) then
      call stat_update_var( iradht_LW, radht, zt )
    end if

  end if


  ! Test scalars with thetal and rt if desired
  if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
  if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

  if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
  if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

  return
  end subroutine atex_tndcy

  !======================================================================
  subroutine atex_sfclyr( ubar, T_sfc, & 
                          thlm_sfc, rtm_sfc, exner_sfc, & 
                          wpthlp_sfc, wprtp_sfc, ustar )
  ! Description:
  !   This subroutine computes surface fluxes of
  !   heat and moisture according to GCSS ATEX specifications

  ! References:
  !   None
  !----------------------------------------------------------------------

  use surface_flux, only: compute_wpthlp_sfc, compute_wprtp_sfc

  implicit none

  ! Input variables
  real, intent(in) ::  &
    ubar,    & ! mean sfc wind speed                           [m/s]
    T_sfc,    & ! Surface temperature                           [K]
    thlm_sfc,& ! theta_l at first model layer                  [K]
    rtm_sfc, & ! Total water mixing ratio at first model layer [kg/kg]
    exner_sfc  ! Exner function                                [-]

  ! Output variables
  real, intent(out) ::  & 
    wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
    wprtp_sfc, &      ! w'rt' surface flux        [(m kg)/(kg s)]
    ustar
    
  ! Local Variable
  real :: & 
    C_10  ! Coefficient

  !-----------------BEGIN CODE-----------------------

  C_10 = 0.0013
  ustar = 0.3

  wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, T_sfc, exner_sfc )
  wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, 0.0198293 )

  return
  end subroutine atex_sfclyr

!----------------------------------------------------------------------
end module atex
