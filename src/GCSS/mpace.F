#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: mpace.F,v 1.3 2007-02-10 16:55:58 griffinb Exp $
        module mpace

!       Description:
!       Contains subroutines for the MPACE intercomparison.
!----------------------------------------------------------------------

        public :: mpace_tndcy, mpace_sfclyr

        private

        contains

!----------------------------------------------------------------------
        subroutine mpace_tndcy
     .  ( time, time_initial, dt, rlat, thlm,
     .    exner, rhot, rtm, p, thvm, rcm,
     .    wmt, wmm, theta_tndcy, rt_tndcy,
     .    Ncnm, Ncm, Frad, radht,
     .    sclrm_forcing )

!        Description:
!          Subroutine to large-scale subsidence for MPACE case (Michael
!          Falk, 21 July 2006).  Added ice and radiation based on Adam
!          Smith Nov 11 case, 27 July 2006.  Comments and documentation
!          added 31 July 2006.
!
!        References:
!          Liou, Wallace and Hobbs, Shettle and Weinman
!-----------------------------------------------------------------------

        use constants
        use model_flags
        use grid_class
        use icedfs
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none

        ! Local constants, subsidence
        real, parameter ::
     .  grav0 = 9.8,     ! m/s
     .  D     = 5.8e-6,  ! 1/s
     .  psfc  = 101000., ! Pa
     .  p00   = 100000., ! Pa; for theta calcs
     .  pinv  = 85000.   ! Pa; ditto

        ! Local constants, LW radiation (from DYCOMS II-RF01)
        real, parameter ::
     .  F0  = 70.0,
     .  F1  = 22.0,
     .  kap = 85.0

        ! Local constants, SW radiation (Shettle and Weinman)
        real, parameter ::
     .  Fs0    = 1212.75,
     .  radius = 1.0e-5,
     .  A      = 0.1,
     .  gc     = 0.86,
     .  omega  = 0.9965
!    .  rlat = 71.75

        ! Local constants, SW radiation (Liou solar angle scheme)
        real, parameter ::
     .  c_0 = 0.006918,
     .  c_1 = -0.399912,
     .  c_2 = -0.006758,
     .  c_3 = -0.002697,
     .  d_1 = 0.070257,
     .  d_2 = 0.000907,
     .  d_3 = 0.000148

        ! Input Variables
        real, intent(in) :: 
     .  time,         ! Current time of simulation      [s]
     .  time_initial, ! Initial time of simulation      [s]
     .  dt,           ! Current length of timestep      [s]
     .  rlat          ! Latitude                        [Degrees North]

        real, dimension(gr%nnzp), intent(in) ::
     .  thlm,  ! Liquid water potential temperature     [K]
     .  exner, ! Exner function                         [-]
     .  rhot,  ! Density of air                         [kg/m^3]
     .  rtm,   ! Total water mixing ratio               [kg/kg]
     .  p,     ! Pressure                               [Pa]
     .  thvm,  ! Virtual potential temperature          [K]
     .  rcm    ! Cloud water mixing ratio               [kg/kg]


        ! Output Variables
        real, dimension(gr%nnzp), intent(out) :: 
     .  wmt,         ! Large-scale vertical motion on t grid   [m/s]
     .  wmm,         ! Large-scale vertical motion on m grid   [m/s]
     .  theta_tndcy, ! Large-scale thlm tendency               [K/s]
     .  rt_tndcy,    ! Large-scale rtm tendency                [kg/kg/s]
     .  Frad,        ! Total radiative flux                    [W/m^2]
     .  Ncm,         ! Cloud droplet number concentration      [count/m^3]
     .  Ncnm,        ! Cloud nuclei number concentration       [count/m^3]
     .  radht        ! dT/dt, then d Theta/dt, due to rad.     [K/s]


        ! Output Variables (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclrm_forcing ! Passive scalar LS tendency            [units/s]


        ! Local Variables, radiation scheme
        real, dimension(gr%nnzp) :: 
     .  radht_LW,! dT/dt, then d Theta/dt, due to LW rad.  [K/s]
     .  radht_SW,! dT/dt, then d Theta/dt, due to SW rad.  [K/s]
     .  Frad_LW, ! Longwave radiative flux                 [W/m^2]
     .  Frad_SW  ! Shortwave radiative flux                [W/m^2]


        ! Local Variables, general
        integer :: i, k ! Loop indices


        ! Local Variables, subsidence scheme
        real ::
     .  velocity_omega


        ! Local Variables, radiation scheme
        real ::
     .  xi_abs,
     .  sda_t,
     .  sda_delta,
     .  sda_h,
     .  t_since_noon,
     .  julday,
     .  start_time_until_noon,
     .  t_tendency

        real, dimension(gr%nnzp) ::
     .  radht_theta,
     .  radht_LW_theta,
     .  radht_SW_theta,
     .  LWP,
     .  thm_rad,
     .  rcm_rad,
     .  rhot_rad,
     .  dsigm,
     .  frad_out,
     .  frad_lw_out,
     .  frad_sw_out,
     .  radhtk,
     .  radht_lw_out,
     .  radht_sw_out

        ! Local variables, simplified ice (icedfs) microphysics
        real, dimension(gr%nnzp) :: dqc_dt_icedfs

        ! Local variables, on/off switches for individual schemes
        logical :: 
     .  subs_on,
     .  lw_on,
     .  sw_on,
     .  center

!-----------------------------------------------------------------------

        ! Set which schemes to use
        subs_on         = .TRUE.
        lw_on           = .TRUE.
        sw_on           = .TRUE.
        center          = .TRUE.

        ! Compute vertical motion
        do i=2,gr%nnzp
          velocity_omega = min( D*(psfc-p(i)), D*(psfc-pinv) )
          wmt(i) = -velocity_omega * Rd * thvm(i) / p(i) / grav0
        end do
        wmm = zt2zm( wmt )


        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        

        ! Compute large-scale tendencies
        do i=1,gr%nnzp
         t_tendency = min( -4.,-15.*(1.-((psfc-p(i))/21818.)) ) ! K/day
         theta_tndcy(i) = (t_tendency * ((psfc/p(i)) ** (Rd/Cp))) 
     .                    / 86400. ! K/s
         rt_tndcy(i)  = min( 0.164,-3*(1-((psfc-p(i))/15171.)) ) / 
     .                 1000. / 86400. ! g/kg/day -> kg/kg/s
        end do

        ! Compute radiation
        julday = 282
        start_time_until_noon = 18840 + 61200
        t_since_noon   = time - start_time_until_noon
        sda_t = 2*3.14*(julday-1)/365

        sda_delta = c_0 + c_1*cos(sda_t) + d_1*sin(sda_t) +
     .              c_2*cos(2*sda_t) + d_2*sin(2*sda_t) +
     .              c_3*cos(3*sda_t) + d_3*sin(3*sda_t)

        sda_h = 2*3.14*t_since_noon/86400

        xi_abs = sin(rlat*3.14/180) * sin(sda_delta) +
     .           cos(rlat*3.14/180) * cos(sda_delta) * cos(sda_h)

        xi_abs = max(xi_abs,0.)

        if (xi_abs == 0.) then
          sw_on = .FALSE.
        end if

        if (.not. sw_on) then
          xi_abs = 0.
        end if

        if ( .not. lbugsrad ) then
          do k = 1, gr%nnzp
            thm_rad(k)  = thlm(gr%nnzp-k+1) + Lv/Cp*rcm(gr%nnzp-k+1)
            rcm_rad(k)  = rcm(gr%nnzp-k+1)
            rhot_rad(k) = rhot(gr%nnzp-k+1)
            dsigm(k)    = 1.0 / gr%dzt(gr%nnzp-k+1)
          end do

          call rad_lwsw(thm_rad, rcm_rad, rhot_rad, dsigm,
     .                  Frad_out, Frad_LW_out, Frad_SW_out,
     .                  radhtk, radht_LW_out, radht_SW_out,
     .                  gr%nnzp-1, center,
     .                  xi_abs, F0, F1, kap, radius, A, gc, Fs0, omega,
     .                  sw_on, lw_on)

          do k = 1, gr%nnzp-1
            Frad(k)     = Frad_out(gr%nnzp-k+1)
            Frad_LW(k)  = Frad_LW_out(gr%nnzp-k+1)
            Frad_SW(k)  = Frad_SW_out(gr%nnzp-k+1)

            radht(k)    = radhtk(gr%nnzp-k+1)
            radht_LW(k) = radht_LW_out(gr%nnzp-k+1)
            radht_SW(k) = radht_SW_out(gr%nnzp-k+1)

            radht_theta(k)    = radht(k) * ((p00/p(k))**(Rd/Cp))
            radht_LW_theta(k) = radht_LW(k) * ((p00/p(k))**(Rd/Cp))
            radht_SW_theta(k) = radht_SW(k) * ((p00/p(k))**(Rd/Cp))
          end do ! k

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
            theta_tndcy(k) = theta_tndcy(k) + radht_theta(k)
          end do

        end if ! ~ lbugsrad

#ifdef STATS
        if ( .not.lbugsrad .and. lstats_samp ) then

          if ( iradht_LW > 0 ) then
            zt%x(:,iradht_LW) = zt%x(:,iradht_LW) + radht_LW
            zt%n(:,iradht_LW) = zt%n(:,iradht_LW) + 1
          end if
          if ( iradht_SW > 0 ) then
            zt%x(:,iradht_SW) = zt%x(:,iradht_SW) + radht_SW
            zt%n(:,iradht_SW) = zt%n(:,iradht_SW) + 1
          end if

          if ( iFrad_SW > 0 ) then
            zm%x(:,iFrad_SW) = zm%x(:,iFrad_SW) + Frad_SW
            zm%n(:,iFrad_SW) = zm%n(:,iFrad_SW) + 1
          end if
          if ( iFrad_LW > 0 ) then
            zm%x(:,iFrad_LW) = zm%x(:,iFrad_LW) + Frad_LW
            zm%n(:,iFrad_LW) = zm%n(:,iFrad_LW) + 1
          end if

        end if
#endif /*STATS*/


        ! Compute icedfs microphysics
        if ( licedfs ) then
          call altocu_icedf( dt, thlm, rcm, p, exner, rhot,
     .                       dqc_dt_icedfs )
          do k=1,gr%nnzp
            theta_tndcy(k) = theta_tndcy(k)
     .                       - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
            rt_tndcy(k) = rt_tndcy(k) + dqc_dt_icedfs(k)
          end do
        end if

        ! Initialize Ncnm on first timestep
        if ( lcoamps_micro .and. time == time_initial ) then
          Ncnm(1:gr%nnzp) 
     .    = 30.0 * (1.0 + exp(-gr%zt(1:gr%nnzp)/2000.0)) * 1.e6

        else if ( kk_rain ) then
          do k=1, gr%nnzp, 1
            if ( rcm(k) >= rc_tol ) then
              Ncm(k) = 30.0 * (1.0 + exp(-gr%zt(k)/2000.0)) * 1.e6
     .                 * rhot(k) 
            end if
          end do
        end if


        ! Set scalar tendencies
        if ( present( sclrm_forcing ) ) then
          sclrm_forcing(:,SCLR_THETA) = theta_tndcy
          sclrm_forcing(:,SCLR_RT)    = rt_tndcy
        end if

        return
        end subroutine mpace_tndcy

!----------------------------------------------------------------------
        subroutine mpace_sfclyr( rho0, um_sfc, vm_sfc,
     .                           upwp_sfc, vpwp_sfc,
     .                           wpthlp_sfc, wprtp_sfc,
     .                           wpsclrp_sfc, wpedsclrp_sfc )

!        Description:
!          Surface forcing subroutine for MPACE case.  Written July-
!          November 2006 by Michael Falk.
!
!        References:
!          MPACE specification, arm.gov
!-----------------------------------------------------------------------

        use constants

        implicit none

        ! External
        intrinsic :: max, sqrt, present

        ! Parameter Constants
        real, parameter ::
     .  ubmin = 0.25,
     .  ustar = 0.25,
        ! Values of HF & LF are from the MPACE specification.
     .  HF    = 136.5, ! Sensible Heat Flux     [W m^-2] 
     .  LF    = 107.7  ! Latent Heat Flux       [W m^-2] 

        ! Input Variables
        real, intent(in)  ::
     .  rho0,    ! Air density at surface       [kg/m^3
     .  um_sfc,  ! um at zt(2)                  [m/s]
     .  vm_sfc   ! vm at zt(2)                  [m/s]

        ! Output Variables
        real, intent(out) :: 
     .  upwp_sfc,    ! u'w' at (1)      [m^2/s^2]
     .  vpwp_sfc,    ! v'w'at (1)       [m^2/s^2]
     .  wpthlp_sfc,  ! w'th_l' at (1)   [(m K)/s]  
     .  wprtp_sfc    ! w'r_t'(1) at (1) [(m kg)/(s kg)]

        ! Output Variables (optional) 
        real, dimension(sclr_dim), optional, intent(out) ::
     .  wpsclrp_sfc,   ! Passive scalar surface flux      [units m/s]
     .  wpedsclrp_sfc  ! Passive eddy-scalar surface flux [units m/s]

        ! Local Variables
        real :: ubar

!-----------------------------------------------------------------------

        ! Compute heat and moisture fluxes
        wpthlp_sfc = HF/(rho0*Cp)
        wprtp_sfc  = LF/(rho0*Lv)

        ! Compute momentum fluxes
        ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

        upwp_sfc = -um_sfc * ustar*ustar / ubar
        vpwp_sfc = -vm_sfc * ustar*ustar / ubar

        ! Compute passive scalar fluxes (set to rt and theta for now)
        if ( present( wpsclrp_sfc ) ) then
          wpsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpsclrp_sfc(SCLR_RT)    = wprtp_sfc
        end if
        if ( present( wpedsclrp_sfc ) ) then
          wpedsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpedsclrp_sfc(SCLR_RT)    = wprtp_sfc
        end if

        return
        end subroutine mpace_sfclyr

        end module mpace
