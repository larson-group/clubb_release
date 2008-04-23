#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: mpace_a.F,v 1.13 2008-04-23 18:09:57 faschinj Exp $
        module mpace_a

!       Description:
!       Contains subroutines for the mpace_a intercomparison.
!----------------------------------------------------------------------

        public :: mpace_a_tndcy, mpace_a_sfclyr, mpace_a_init

        private

    
        ! These variables were moved here so that they could be 
        ! accessible to all subroutines in mpace_a
        ! Joshua Fasching December 2007

        integer, parameter :: file_ntimes = 139
        integer, parameter :: file_nlevels = 38
        integer, parameter :: per_line = 5
        
        ! File path for forcing the forcing files 
        character(*), parameter :: 
     .     file_path = '../model/mpace_a_forcings/'

        real, dimension(file_nlevels) :: file_pressure
        real, dimension(file_nlevels) :: file_heights
        real, dimension(file_ntimes) :: file_times


! Michael Falk is, on 28 September 2007, removing omega.  We are going to try
! to force the model without specifying it, so we can do the temperature and
! moisture forcings the way Steve Klein wants us to.

!       real, dimension(file_nlevels,file_ntimes) :: omega_forcing ! mb/s
        real, dimension(file_nlevels,file_ntimes) :: dTdt_forcing  ! K/hr
        real, dimension(file_nlevels,file_ntimes) :: dqdt_forcing  ! g/kg/hr
        real, dimension(file_nlevels,file_ntimes) :: vertT_forcing  ! K/hr
        real, dimension(file_nlevels,file_ntimes) :: vertq_forcing  ! g/kg/hr
        real, dimension(file_nlevels,file_ntimes) :: um_obs  ! m/s
        real, dimension(file_nlevels,file_ntimes) :: vm_obs  ! m/s
        real, dimension(file_ntimes) :: file_LH
        real, dimension(file_ntimes) :: file_SH

        contains

!----------------------------------------------------------------------
        subroutine mpace_a_tndcy
     .  ( time, time_initial, dt, rlat, thlm,
     .    exner, rhot, rtm, p, thvm, rcm,
     .    wmt, wmm, theta_tndcy, rt_tndcy,
     .    Ncnm, Ncm, Frad, radht,
     .    um_hoc_grid, vm_hoc_grid,
     .    sclrm_forcing )

!        Description:
!
!        References:
!          Liou, Wallace and Hobbs, Shettle and Weinman
!-----------------------------------------------------------------------

        use constants
        use model_flags
        use grid_class
        use interpolation, only: zlinterp_fnc
        use file_functions, only: file_read_1d, file_read_2d
        use stats_prec, only: time_prec
#ifdef STATS
        use stats_hoc
#endif /*STATS*/

        implicit none

        ! Local constants, subsidence
        real, parameter ::
!     .  grav0 = 9.8,     ! m/s
     .  D     = 5.8e-6,  ! 1/s
     .  psfc  = 101000., ! Pa
     .  pinv  = 85000.   ! Pa; ditto

        ! Local constants, LW radiation (from DYCOMS II-RF01)
        real, parameter ::
     .  F0  = 100.0,
     .  F1  = 20.0,
     .  kap = 119.0

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
        real(kind=time_prec), intent(in) :: 
     .  time,         ! Current time of simulation      [s]
     .  time_initial  ! Initial time of simulation      [s]

        real(kind=time_prec), intent(in) :: 
     .  dt            ! Current length of timestep      [s]

        real, intent(in) :: 
     .  rlat          ! Latitude                        [Degrees North]

        real, dimension(gr%nnzp), intent(in) ::
     .  thlm,  ! Liquid water potential temperature     [K]
     .  exner, ! Exner function                         [-]
     .  rhot,  ! Density of air                         [kg/m^3]
     .  rtm,   ! Total water mixing ratio               [kg/kg]
     .  p,     ! Pressure                               [Pa]
     .  thvm,  ! Virtual potential temperature          [K]
     .  rcm    ! Cloud water mixing ratio               [kg/kg]

        ! Input/Output Variables
        real, dimension(gr%nnzp), intent(inout) :: 
     .  Ncm,         ! Cloud droplet number concentration      [count/m^3]
     .  Ncnm         ! Cloud nuclei number concentration       [count/m^3]

        ! Output Variables
        real, dimension(gr%nnzp), intent(out) :: 
     .  wmt,         ! Large-scale vertical motion on t grid   [m/s]
     .  wmm,         ! Large-scale vertical motion on m grid   [m/s]
     .  theta_tndcy, ! Large-scale thlm tendency               [K/s]
     .  rt_tndcy,    ! Large-scale rtm tendency                [kg/kg/s]
     .  Frad,        ! Total radiative flux                    [W/m^2]
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
        integer :: i, j, k ! Loop indices


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

        ! Local variables, on/off switches for individual schemes
        logical :: 
     .  subs_on,
     .  lw_on,
     .  sw_on,
     .  center

! Open external files (21 Aug 2007, Michael Falk)

      integer left_time,right_time
      real :: ratio

!      real, dimension(file_nlevels) :: omega_column
      real, dimension(file_nlevels) :: dTdt_column
      real, dimension(file_nlevels) :: dqdt_column
      real, dimension(file_nlevels) :: vertT_column
      real, dimension(file_nlevels) :: vertq_column
      real, dimension(file_nlevels) :: um_column
      real, dimension(file_nlevels) :: vm_column

!      real, dimension(gr%nnzp) :: omega_hoc_grid
      real, dimension(gr%nnzp) :: dTdt_hoc_grid
      real, dimension(gr%nnzp) :: dqdt_hoc_grid
      real, dimension(gr%nnzp) :: vertT_hoc_grid
      real, dimension(gr%nnzp) :: vertq_hoc_grid

        real, dimension(gr%nnzp), intent(out) :: 
     .  um_hoc_grid,      ! Observed wind, for nudging         [m/s]
     .  vm_hoc_grid       ! Observed wind, for nudging         [m/s]

      ! This code block takes the model time, finds the time before it and the time after it on the list,
      ! and marks them left_time and right_time for interpolation.  If the time is before the first or
      ! after the last time in the file, it just uses the first or last time without interpolation.
      if (time <= file_times(1)) then
        print *,'Time is at or before the first time in the list.'
        left_time = 1
        right_time = 1
      else if (time >= file_times(file_ntimes)) then
        print *,'Time is at or after the last time in the list.'
        left_time = file_ntimes
        right_time = file_ntimes
      else
        do k=1,file_ntimes-1
          if ((time > file_times(k)) .AND.
     .        (time <=file_times(k+1))) then
            left_time = k
            right_time = k+1
          end if
        end do
      end if

      ratio = (time - file_times(left_time)) /                 ! This is the ratio "a" needed for linear interpolation in time.
     .        (file_times(right_time) - file_times(left_time)) ! at the first time a=0; at the second time a=1.

      do k=1,file_nlevels
!        omega_column(k) = ratio *			       ! Do linear interpolation in time
!     .                      (omega_forcing(k,right_time)
!     .                      -omega_forcing(k,left_time))
!     .                     + omega_forcing(k,left_time)        

        dTdt_column(k) = ratio *                               ! Do linear interpolation in time
     .                      (dTdt_forcing(k,right_time)
     .                      -dTdt_forcing(k,left_time))
     .                     + dTdt_forcing(k,left_time)        

        dqdt_column(k) = ratio *                               ! Do linear interpolation in time
     .                      (dqdt_forcing(k,right_time)
     .                      -dqdt_forcing(k,left_time))
     .                     + dqdt_forcing(k,left_time)        

        vertT_column(k) = ratio *                               ! Do linear interpolation in time
     .                      (vertT_forcing(k,right_time)
     .                      -vertT_forcing(k,left_time))
     .                     + vertT_forcing(k,left_time)        

        vertq_column(k) = ratio *                               ! Do linear interpolation in time
     .                      (vertq_forcing(k,right_time)
     .                      -vertq_forcing(k,left_time))
     .                     + vertq_forcing(k,left_time)        

        um_column(k) = ratio *                                  ! Do linear interpolation in time
     .                      (um_obs(k,right_time)
     .                      -um_obs(k,left_time))
     .                     + um_obs(k,left_time)        

        vm_column(k) = ratio *                                  ! Do linear interpolation in time
     .                      (vm_obs(k,right_time)
     .                      -vm_obs(k,left_time))
     .                     + vm_obs(k,left_time)        

      end do

!      omega_hoc_grid = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt,    ! Do linear interpolation in space
!     .                         file_heights,omega_column)            ! using zlinterp_fnc
      dTdt_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt,
     .                         file_heights,dTdt_column)
      dqdt_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt,
     .                         file_heights,dqdt_column)
      vertT_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt,
     .                         file_heights,vertT_column)
      vertq_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt,
     .                         file_heights,vertq_column)
      um_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt,
     .                         file_heights,um_column)
      vm_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt,
     .                         file_heights,vm_column)

      um_hoc_grid (1) = um_hoc_grid(2)
      vm_hoc_grid (1) = vm_hoc_grid(2)

! eMFc

!-----------------------------------------------------------------------

        ! Set which schemes to use
        subs_on         = .TRUE.
        lw_on           = .TRUE.
        sw_on           = .TRUE.
        center          = .TRUE.

        ! Compute vertical motion
        do i=2,gr%nnzp
!          velocity_omega = omega_hoc_grid(i) * 100 / 3600 ! convering mb/hr to Pa/s
!          wmt(i) = -velocity_omega * Rd * thvm(i) / p(i) / grav
           wmt(i) = 0.
! End of Michael Falk's obliteration of omega.
        end do

        ! Boundary condition
        wmt(1) = 0.0        ! Below surface

        ! Interpolation
        wmm = zt2zm( wmt )

        ! Boundary condition
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        

        ! Compute large-scale tendencies
        do i=1,gr%nnzp
         theta_tndcy(i) = ((dTdt_hoc_grid(i) + vertT_hoc_grid(i))
     .                    * ((psfc/p(i)) ** (Rd/Cp)))
     .                    / 3600. ! K/s
         rt_tndcy(i)  = (dqdt_hoc_grid(i)+vertq_hoc_grid(i))
     .    / 1000. / 3600. ! g/kg/hr -> kg/kg/s
        end do

        ! Compute radiation
        julday = 282
        start_time_until_noon = 29051 + 50400
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

          call rad_lwsw( thm_rad, rcm_rad, rhot_rad, dsigm,
     .                   Frad_out, Frad_LW_out, Frad_SW_out,
     .                   radhtk, radht_LW_out, radht_SW_out,
     .                   gr%nnzp-1, center,
     .                   xi_abs, F0, F1, kap, radius, A, gc, Fs0, omega,
     .                   sw_on, lw_on )

          do k = 2, gr%nnzp-1
            Frad(k)     = Frad_out(gr%nnzp-k+1)
            Frad_LW(k)  = Frad_LW_out(gr%nnzp-k+1)
            Frad_SW(k)  = Frad_SW_out(gr%nnzp-k+1)

            radht(k)    = radhtk(gr%nnzp-k+1)
            radht_LW(k) = radht_LW_out(gr%nnzp-k+1)
            radht_SW(k) = radht_SW_out(gr%nnzp-k+1)

            radht_theta(k)    = radht(k) * ((p0/p(k))**(Rd/Cp))
            radht_LW_theta(k) = radht_LW(k) * ((p0/p(k))**(Rd/Cp))
            radht_SW_theta(k) = radht_SW(k) * ((p0/p(k))**(Rd/Cp))
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

        ! Initialize Ncnm on first timestep
        if ( lcoamps_micro .and. time == time_initial ) then
          Ncnm(1:gr%nnzp)
     .    = 30.0 * (1.0 + exp(-gr%zt(1:gr%nnzp)/2000.0)) * 1.e6

        else if ( kk_rain ) then
          ! Note: Khairoutdinov and Kogan microphysics has only been
          ! tested for marine stratocumulous clouds, and does not
          ! account for snow and ice.
          do k=1, gr%nnzp, 1
            if ( rcm(k) >= rc_tol ) then
              Ncm(k) = 30.0 * (1.0 + exp(-gr%zt(k)/2000.0)) * 1.e6
     .                 / rhot(k) 
            end if
          end do
        end if


        ! Set scalar tendencies
        if ( present( sclrm_forcing ) ) then
          sclrm_forcing(:,SCLR_THETA) = theta_tndcy
          sclrm_forcing(:,SCLR_RT)    = rt_tndcy
        end if

        return
        end subroutine mpace_a_tndcy

!----------------------------------------------------------------------
        subroutine mpace_a_sfclyr( time, rho0, um_sfc, vm_sfc,
     .                             upwp_sfc, vpwp_sfc,
     .                             wpthlp_sfc, wprtp_sfc, ustar,
     .                             wpsclrp_sfc, wpedsclrp_sfc )

!        Description:
!          Surface forcing subroutine for mpace_a case.  Written 
!          October 2007 by Michael Falk.
!
!        References:
!          mpace_a specification from arm.gov
!-----------------------------------------------------------------------

        use constants
        use file_functions, only: file_read_1d, file_read_2d
        use stats_prec, only: time_prec

        implicit none

        ! External
        intrinsic :: max, sqrt, present

        ! Parameter Constants
        real, parameter ::
     .  ubmin = 0.25
!     .  ustar = 0.25

        ! Input Variables
        real(kind=time_prec), intent(in) ::
     .  time     ! current model time           [s]

        real, intent(in)  ::
     .  rho0,    ! Air density at surface       [kg/m^3
     .  um_sfc,  ! um at zt(2)                  [m/s]
     .  vm_sfc   ! vm at zt(2)                  [m/s]

        ! Output Variables
        real, intent(out) :: 
     .  upwp_sfc,    ! u'w' at (1)      [m^2/s^2]
     .  vpwp_sfc,    ! v'w'at (1)       [m^2/s^2]
     .  wpthlp_sfc,  ! w'th_l' at (1)   [(m K)/s]  
     .  wprtp_sfc,   ! w'r_t'(1) at (1) [(m kg)/(s kg)]
     .  ustar        ! surface friction velocity [m/s]

        ! Output Variables (optional) 
        real, dimension(sclr_dim), optional, intent(out) ::
     .  wpsclrp_sfc,   ! Passive scalar surface flux      [units m/s]
     .  wpedsclrp_sfc  ! Passive eddy-scalar surface flux [units m/s]

        ! Local Variables
         real ::
     .   ubar,
     .   latent_heat_flx,
     .   sensible_heat_flx

         integer ::
     .   i,j,k


      integer left_time,right_time

!-----------------------------------------------------------------------

      ! choose which times to use
      if (time <= file_times(1)) then
        print *,'Time is at or before the first time in the list.'
        left_time = 1
        right_time = 1
      else if (time >= file_times(file_ntimes)) then
        print *,'Time is at or after the last time in the list.'
        left_time = file_ntimes
        right_time = file_ntimes
      else
        do k=1,file_ntimes-1
          if ((time > file_times(k)) .AND.
     .        (time <=file_times(k+1))) then
            left_time = k
            right_time = k+1
          end if
        end do
      end if

      latent_heat_flx = ((time-file_times(left_time)) / 
     .     (file_times(right_time)-file_times(left_time))) *
     .     (file_LH(right_time)
     .     -file_LH(left_time))
     .     + file_LH(left_time)        

      sensible_heat_flx = ((time-file_times(left_time)) / 
     .     (file_times(right_time)-file_times(left_time))) *
     .     (file_SH(right_time)
     .     -file_SH(left_time))
     .     + file_SH(left_time)        

        ! Compute heat and moisture fluxes
        wpthlp_sfc = sensible_heat_flx/(rho0*Cp)
        wprtp_sfc  = latent_heat_flx/(rho0*Lv)

        ! Compute momentum fluxes
        ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

        ! Declare the value of ustar.
        ustar = 0.25

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
        end subroutine mpace_a_sfclyr
!----------------------------------------------------------------
        subroutine mpace_a_init()
!
!       Description:
!       This subroutine initializes the module by reading in forcing
!       data used in the tndcy and sfclyr subroutines.
!----------------------------------------------------------------
          use file_functions, only: file_read_1d, file_read_2d
          implicit none
          call file_read_1d(10,
     .      file_path//'mpace_a_press.dat',
     .      file_nlevels, per_line, file_pressure)

          call file_read_1d(10,
     .      file_path//'mpace_a_heights.dat',
     .      file_nlevels, per_line, file_heights)

          call file_read_1d(10,
     .      file_path//'mpace_a_times.dat',
     .      file_ntimes, per_line, file_times)

!      call file_read_2d(10,
!     . file_path//'mpace_a_omega.dat',
!     . file_nlevels, file_ntimes, per_line, omega_forcing)

          call file_read_2d(10,
     .      file_path//'mpace_a_dTdt.dat',
     .      file_nlevels, file_ntimes, per_line, dTdt_forcing)

          call file_read_2d(10,
     .      file_path//'mpace_a_dqdt_horiz.dat',
     .      file_nlevels, file_ntimes, per_line, dqdt_forcing)

          call file_read_2d(10,
     .      file_path//'mpace_a_verts.dat',
     .      file_nlevels, file_ntimes, per_line, vertT_forcing)

          call file_read_2d(10,
     .      file_path//'mpace_a_vertq.dat',
     .      file_nlevels, file_ntimes, per_line, vertq_forcing)

          call file_read_2d(10,
     .      file_path//'mpace_a_um_obs.dat',
     .      file_nlevels, file_ntimes, per_line, um_obs)

          call file_read_2d(10,
     .      file_path//'mpace_a_vm_obs.dat',
     .      file_nlevels, file_ntimes, per_line, vm_obs)

          call file_read_1d(10,
     .      file_path//'mpace_a_lh.dat',
     .      file_ntimes, per_line, file_LH)

          call file_read_1d(10,
     .      file_path//'mpace_a_sh.dat',
     .      file_ntimes, per_line, file_SH)

           
        end subroutine mpace_a_init

        end module mpace_a
