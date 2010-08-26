!----------------------------------------------------------------------
! $Id$
module dycoms2_rf01

!       Description:
!       Contains subroutines for the DYCOMS II RF01 case.
!----------------------------------------------------------------------
  implicit none

  public :: dycoms2_rf01_tndcy, dycoms2_rf01_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine dycoms2_rf01_tndcy( rho, rho_zm, rtm, rcm, exner, & 
                                 err_code, & 
                                 Frad, radht, & 
                                 thlm_forcing, rtm_forcing, &
                                 sclrm_forcing, edsclrm_forcing )
    !       Description:
    !       Subroutine to set theta and water tendencies for DYCOMS RF01 case.

    !       References:
    !----------------------------------------------------------------------

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm, ddzm ! Procedure(s)

    use constants_clubb, only: fstderr, Cp ! Constant(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use parameters_radiation, only: rad_scheme ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use error_code, only: clubb_rtm_level_not_found ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variables(s)

    use interpolation, only: lin_int ! Procedure(s)

    use stats_type, only: stat_update_var, stat_update_var_pt ! Procedure(s)

    use stats_variables, only:  & 
        iz_inversion, iradht_LW, zt, sfc, l_stats_samp ! Variable(s)

    use parameters_radiation, only: &
      F0,  & ! Variable(s)
      F1,  &
      kappa

    implicit none

    ! External
    intrinsic :: exp, sqrt, present

    ! Constant Parameters
    real, parameter ::  & 
      ls_div =  3.75e-6
!     F0 = 70.0, F1 = 22.0,  &
!     kay = 85.0

    ! Input Variables
    real, dimension(gr%nnzp), intent(in) ::  & 
      rho_zm,  & ! Density on moment. grid         [kg/m^3]
      rho,     & ! Density on thermo. grid         [kg/m^3] 
      rtm,     & ! Total water mixing ratio        [kg/kg]
      rcm,     & ! Cloud water mixing ratio        [kg/kg]
      exner      ! Exner function                  [-]

    ! Input/Output Variables
    integer, intent(inout) :: err_code

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) ::  & 
      radht,         & ! Radiative heating rate                       [K/s]
      Frad,          & ! Radiative flux                               [W/m^2]
      thlm_forcing,  & ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing      ! Total water mixing ratio tendency            [kg/kg/s]

    real, intent(out), dimension(gr%nnzp, sclr_dim) :: & 
      sclrm_forcing   ! Passive scalar tendency         [units/s]

    real, intent(out), dimension(gr%nnzp, edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar tendency    [units/s]

    ! Local variables
    real, dimension(gr%nnzp) :: lwp

    integer :: i
    real :: z_inversion

    thlm_forcing = 0.
    rtm_forcing  = 0.

    ! Identify height of 8.0 g/kg moisture level

    i = 2
    do while ( i <= gr%nnzp .and. rtm(i) > 8.0e-3 )
      i = i + 1
    end do
    if ( i == gr%nnzp+1 .or. i == 2 ) then
      write(fstderr,*) "Identification of 8.0 g/kg level failed"
      write(fstderr,*) "Subroutine: dycoms2_rf01_tndcy.  " & 
        //" File: dycoms2_rf01.F"
      write(fstderr,*) "i = ",i
      write(fstderr,*) "rtm(i) = ",rtm(i)
      err_code = clubb_rtm_level_not_found
      return
    end if
!z_inversion = (gr%zt(i)-gr%zt(i-1))/(rtm(i)-rtm(i-1))*(8.0e-3-rtm(i-1)) &
!   + gr%zt(i-1)
!        x_sfc(1,iz_inversion) = z_inversion
    z_inversion = lin_int( 8.0e-3, rtm(i), rtm(i-1), gr%zt(i), gr%zt(i-1) )

    if ( l_stats_samp ) then
      call stat_update_var_pt( iz_inversion, 1, z_inversion, sfc )
    end if

    ! Theta-l radiative tendency

    if ( trim( rad_scheme ) == "simplified" ) then

      ! Compute liquid water path from top of the model
      ! We define liquid water path on momentum levels

      lwp(gr%nnzp) = 0.0
      do i = gr%nnzp-1, 1, -1
        lwp(i) = lwp(i+1) + rho(i+1) * rcm(i+1) / gr%invrs_dzt(i+1)
      end do
!         x_sfc(1,ilwp) = lwp(1)

      ! Compute IR radiative flux

      do i = 1, gr%nnzp, 1
        Frad(i) = F0 * EXP( -kappa * lwp(i) ) & 
                + F1 * EXP( -kappa * (lwp(1)-lwp(i)) )
        if ( z_inversion > 0 .and. gr%zm(i) > z_inversion ) then
          Frad(i) = Frad(i) & 
                  + rho_zm(i) * cp * ls_div & 
                    * ( 0.25*(gr%zm(i)-z_inversion)**(4.0/3.0) & 
                        + z_inversion*(gr%zm(i)-z_inversion)**(1.0/3.0) )
        end if
      end do

      ! Compute IR heating rate

      radht          = ( -1.0/(Cp*rho) ) * ddzm( Frad ) & 
                     * ( 1.0 / exner )
      radht(1)       = 0.
      radht(gr%nnzp) = 0.

      if ( l_stats_samp ) then
        call stat_update_var( iradht_LW, radht, zt )
      end if
    else

      radht = 0.0 ! Computed elsewhere

    end if ! simplified

    ! Add heating rate to theta-l forcing

    thlm_forcing = thlm_forcing + radht

    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

    return
  end subroutine dycoms2_rf01_tndcy
  
  !======================================================================
  subroutine dycoms2_rf01_sfclyr( sfctype, T_sfc, p_sfc,  & 
                                    exner_sfc, ubar, & 
                                    thlm_sfc, rtm_sfc, rho_zm_sfc, &
                                    wpthlp_sfc, wprtp_sfc, ustar )
  ! Description:
  !   This subroutine computes surface fluxes of
  !   heat and moisture according to GCSS DYCOMS II RF 01 specifications

  ! References:
  !   None
  !----------------------------------------------------------------------
  use constants_clubb, only: Cp, fstderr, Lv ! Variable(s)

  use saturation, only: sat_mixrat_liq ! Variable(s)

  use surface_flux, only: compute_wpthlp_sfc, compute_wprtp_sfc ! Procedure(s)

  implicit none

  ! Input variables
  integer, intent(in) :: &
    sfctype
  real, intent(in) ::  &
    T_sfc,      & ! Surface temperature                           [K]
    p_sfc,      & ! Surface pressure                              [Pa]
    exner_sfc, & ! Exner function                                [-]
    ubar,      & ! mean sfc wind speed                           [m/s]
    thlm_sfc,  & ! theta_l at first model layer                  [K]
    rtm_sfc,   & ! Total water mixing ratio at first model layer [kg/kg]
    rho_zm_sfc   ! Density at the surface                        [kg/m^3]

  ! Output variables
  real, intent(out) ::  & 
    wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
    wprtp_sfc, &      ! w'rt' surface flux        [(m kg)/(kg s)]
    ustar
    
  ! Local Variable
  real :: & 
    Cd  ! Coefficient

  !-----------------BEGIN CODE-----------------------

  Cd = 0.0011

  ustar = 0.25

  ! Compute heat and moisture fluxes
  if ( sfctype == 0 ) then

    wpthlp_sfc =  15.0 / ( rho_zm_sfc * Cp )
    wprtp_sfc  = 115.0 / ( rho_zm_sfc * Lv )

  else if ( sfctype == 1 ) then

    wpthlp_sfc = compute_wpthlp_sfc( Cd, ubar, thlm_sfc, T_sfc, exner_sfc )
    wprtp_sfc = compute_wprtp_sfc( Cd, ubar, rtm_sfc, sat_mixrat_liq( p_sfc, T_sfc ) )

  else  ! Undefined value for sfctype

    write(fstderr,*) "Invalid sfctype value = ", sfctype
    stop

  end if ! sfctype
  return
  end subroutine dycoms2_rf01_sfclyr

!----------------------------------------------------------------------
end module dycoms2_rf01
