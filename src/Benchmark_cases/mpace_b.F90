!----------------------------------------------------------------------
! $Id$
module mpace_b

!       Description:
!       Contains subroutines for the mpace_b intercomparison.
!----------------------------------------------------------------------

  implicit none

  public :: mpace_b_tndcy, mpace_b_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine mpace_b_tndcy( xi_abs, & 
                            rho, p_in_Pa, thvm, rcm, & 
                            wm_zt, wm_zm, thlm_forcing, rtm_forcing, & 
                            Frad, radht, &
                            sclrm_forcing, edsclrm_forcing )

  !        Description:
  !          Subroutine to large-scale subsidence for mpace_b case (Michael
  !          Falk, 21 July 2006).  Added ice and radiation based on Adam
  !          Smith Nov 11 case, 27 July 2006.  Comments and documentation
  !          added 31 July 2006.
  !
  !        References:
  !          Liou, Wallace and Hobbs, Shettle and Weinman
  !-----------------------------------------------------------------------

    use constants_clubb, only: Rd, Cp, Lv, p0, rc_tol, zero_threshold ! Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    use rad_lwsw_module, only: rad_lwsw ! Variable(s)

    use array_index, only: iiedsclr_rt, iiedsclr_thl, iisclr_rt, iisclr_thl ! Variable(s)

    use stats_type, only: stat_update_var ! Procedure(s)

    use stats_variables, only: iFrad_LW, iFrad_SW, iradht_SW,  & ! Variable(s)
                   iradht_LW, zt, zm, l_stats_samp

    use parameters_radiation, only: &
      F0, F1, Fs_values, alvdr, kappa, omega, gc, eff_drop_radius, rad_scheme
     

    implicit none

    ! Local constants, subsidence
    real, parameter :: & 
      grav0 = 9.8,     & ! m/s
      D     = 5.8e-6,  & ! 1/s
      p_sfc  = 101000., & ! Pa
      pinv  = 85000.     ! Pa; ditto

    ! Local constants, LW radiation (from DYCOMS II-RF01)
    !real, parameter :: & 
    !  F0  = 70.0, & 
    !  F1  = 22.0, & 
    !  kap = 85.0

    ! Local constants, SW radiation (Shettle and Weinman)
    !real, parameter :: & 
    !  Fs0    = 1212.75, & 
    !  radius = 1.0e-5, & 
    !  A      = 0.1, & 
    !  gc     = 0.86, & 
    !  omega  = 0.9965

    ! Input Variables
    real, intent(in) ::  & 
      xi_abs ! Cosine of the solar zenith angle [-]

    real, dimension(gr%nnzp), intent(in) :: & 
      rho,     & ! Density of air                         [kg/m^3]
      p_in_Pa, & ! Pressure                               [Pa]
      thvm,    & ! Virtual potential temperature          [K]
      rcm        ! Cloud water mixing ratio               [kg/kg]

    ! Output Variables
    real, dimension(gr%nnzp), intent(out) ::  & 
      wm_zt,        & ! Large-scale vertical motion on t grid   [m/s]
      wm_zm,        & ! Large-scale vertical motion on m grid   [m/s]
      thlm_forcing, & ! Large-scale thlm tendency               [K/s]
      rtm_forcing,  & ! Large-scale rtm tendency                [kg/kg/s]
      Frad,         & ! Total radiative flux                    [W/m^2]
      radht           ! dT/dt, then d Theta/dt, due to rad.     [K/s]

    real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar LS tendency            [units/s]

    real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar LS tendency     [units/s]


    ! Local Variables, radiation scheme
    real, dimension(gr%nnzp) ::  & 
      radht_LW, & ! dT/dt, then d Theta/dt, due to LW rad.  [K/s]
      radht_SW, & ! dT/dt, then d Theta/dt, due to SW rad.  [K/s]
      Frad_LW,  & ! Longwave radiative flux                 [W/m^2]
      Frad_SW     ! Shortwave radiative flux                [W/m^2]


    ! Local Variables, general
    integer :: i, k ! Loop indices


    ! Local Variables, subsidence scheme
    real :: & 
      velocity_omega


    ! Local Variables
    real :: t_tendency

    real, dimension(gr%nnzp) :: & 
      radht_theta, & 
      radht_LW_theta, & 
      radht_SW_theta, & 
    !  LWP,            & ! Liquid water path                             [kg/m^2]
      rcm_rad,        & ! Flipped array of liq. water mixing ratio       [kg/kg]
      rho_rad,        & ! Flipped array of air density                   [kg/m^3]
      dsigm,          & ! Flipped array of grid spacing                  [m]
      coamps_zm,      & ! Flipped array of momentum level altitudes      [m]
      coamps_zt,      & ! Flipped array of thermodynamic level altitudes [m]
      frad_out,       & ! Flipped array of radiaive flux                 [W/m^2]
      frad_lw_out,    & ! Flipped array of LW radiative flux             [W/m^2]
      frad_sw_out,    & ! Flipped array of SW radiative flux             [W/m^2]
      radhtk,         & ! Flipped array of radiative heating             [K/s]
      radht_lw_out,   & ! Flipped array of LW radiative heating          [K/s]
      radht_sw_out      ! Flipped array of SW radiative heating          [K/s]

    ! Local variables, on/off switches for individual schemes
    logical ::  & 
      l_lw_on, & 
      l_sw_on, & 
      !l_subs_on, &
      l_center

    !-----------------------------------------------------------------------

    ! Set which schemes to use
    l_lw_on           = .TRUE.
    l_sw_on           = .TRUE.
    !l_subs_on         = .TRUE.
    l_center          = .TRUE.

    ! Compute vertical motion
    do i=2,gr%nnzp
      velocity_omega = min( D*(p_sfc-p_in_Pa(i)), D*(p_sfc-pinv) )
      wm_zt(i) = -velocity_omega * Rd * thvm(i) / p_in_Pa(i) / grav0
    end do



    ! Boundary condition
    wm_zt(1) = 0.0        ! Below surface

    ! Interpolate
    wm_zm = zt2zm( wm_zt )

    ! Boundary conditions
    wm_zm(1) = 0.0        ! At surface
    wm_zm(gr%nnzp) = 0.0  ! Model top


    ! Compute large-scale tendencies
    do i=1,gr%nnzp
     t_tendency = min( -4.,-15.*(1.-((p_sfc-p_in_Pa(i))/21818.)) ) ! K/day
     thlm_forcing(i) = (t_tendency * ((p_sfc/p_in_Pa(i)) ** (Rd/Cp)))  & 
                      / 86400. ! K/s
     rtm_forcing(i)  = min( 0.164,-3*(1-((p_sfc-p_in_Pa(i))/15171.)) ) /  & 
                   1000. / 86400. ! g/kg/day -> kg/kg/s
    end do

    if ( trim( rad_scheme ) == "simplified" ) then

      do k = 1, gr%nnzp
        rcm_rad(k)  = rcm(gr%nnzp-k+1)
        rho_rad(k) = rho(gr%nnzp-k+1)
        dsigm(k)    = 1.0 / gr%invrs_dzt(gr%nnzp-k+1)
        coamps_zm(k) = gr%zm(gr%nnzp-k+1)
        coamps_zt(k) = gr%zt(gr%nnzp-k+1)
      end do

      if ( xi_abs == 0. ) then
        l_sw_on = .false.
      end if

      call rad_lwsw(rcm_rad, rho_rad, dsigm, & 
                  coamps_zm, coamps_zt, & 
                  Frad_out, Frad_LW_out, Frad_SW_out, & 
                  radhtk, radht_LW_out, radht_SW_out, & 
                  gr%nnzp-1, l_center, & 
                  xi_abs, F0, F1, kappa, eff_drop_radius, real( alvdr ), gc, Fs_values(1), omega, &
                  l_sw_on, l_lw_on)

      do k = 2, gr%nnzp-1
        Frad(k)     = Frad_out(gr%nnzp-k+1)
        Frad_LW(k)  = Frad_LW_out(gr%nnzp-k+1)
        Frad_SW(k)  = Frad_SW_out(gr%nnzp-k+1)

        radht(k)    = radhtk(gr%nnzp-k+1)
        radht_LW(k) = radht_LW_out(gr%nnzp-k+1)
        radht_SW(k) = radht_SW_out(gr%nnzp-k+1)

        radht_theta(k)    = radht(k) * ((p0/p_in_Pa(k))**(Rd/Cp))
        radht_LW_theta(k) = radht_LW(k) * ((p0/p_in_Pa(k))**(Rd/Cp))
        radht_SW_theta(k) = radht_SW(k) * ((p0/p_in_Pa(k))**(Rd/Cp))
      end do ! k = 2..gr%nnzp

      Frad(1)    = Frad(2)
      Frad_LW(1) = Frad_LW(2)
      Frad_SW(1) = Frad_SW(2)
      radht_theta(1)    = radht_theta(2)
      radht_LW_theta(1) = radht_LW_theta(2)
      radht_SW_theta(1) = radht_SW_theta(2)

      Frad(gr%nnzp)    = Frad(gr%nnzp-1)
      Frad_LW(gr%nnzp) = Frad_LW(gr%nnzp-1)
      Frad_SW(gr%nnzp) = Frad_SW(gr%nnzp-1)
      radht_theta(gr%nnzp)    = radht_theta(gr%nnzp-1)
      radht_LW_theta(gr%nnzp) = radht_LW_theta(gr%nnzp-1)
      radht_SW_theta(gr%nnzp) = radht_SW_theta(gr%nnzp-1)

      radht(1:gr%nnzp)    = radht_theta(1:gr%nnzp)
      radht_LW(1:gr%nnzp) = radht_LW_theta(1:gr%nnzp)
      radht_SW(1:gr%nnzp) = radht_SW_theta(1:gr%nnzp)

      do k = 1, gr%nnzp
        thlm_forcing(k) = thlm_forcing(k) + radht_theta(k)
      end do

      if ( l_stats_samp ) then
     
        call stat_update_var( iradht_LW, radht_LW, zt )

        call stat_update_var( iradht_SW, radht_SW, zt )

        call stat_update_var( iFrad_SW, Frad_SW, zm )

        call stat_update_var( iFrad_LW, Frad_LW, zm )

      end if

    end if ! ~ l_bugsrad


    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

    return
  end subroutine mpace_b_tndcy

!----------------------------------------------------------------------
  subroutine mpace_b_sfclyr( rho0, & 
                           wpthlp_sfc, wprtp_sfc, ustar )

  !        Description:
  !          Surface forcing subroutine for mpace_b case.  Written July-
  !          November 2006 by Michael Falk.
  !
  !        References:
  !          mpace_b specification, arm.gov
  !-----------------------------------------------------------------------

    use constants_clubb, only: Cp, Lv ! Variable(s)

    implicit none

    ! External
    intrinsic :: max, sqrt

    ! Parameter Constants
    real, parameter :: &
    ! The values of these are from the mpace_b specification.
      sensible_heat_flx  = 136.5,  & ! Sensible Heat Flux     [W m^-2] 
      latent_heat_flx    = 107.7     ! Latent Heat Flux       [W m^-2] 
    ! eMFc

    ! Input Variables
    real, intent(in)  :: & 
      rho0    ! Air density at surface       [kg/m^3

    ! Output Variables
    real, intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    !-----------------------------------------------------------------------

    ! Declare the value of ustar.
    ustar = 0.25

    ! Compute heat and moisture fluxes
    wpthlp_sfc = sensible_heat_flx/(rho0*Cp)
    wprtp_sfc  = latent_heat_flx/(rho0*Lv)

    return
  end subroutine mpace_b_sfclyr

end module mpace_b
