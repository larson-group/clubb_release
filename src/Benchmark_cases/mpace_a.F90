!----------------------------------------------------------------------
! $Id$
  module mpace_a

! Description:
!   Contains subroutines for the mpace_a intercomparison.
!----------------------------------------------------------------------

  implicit none

  public :: mpace_a_tndcy, mpace_a_sfclyr, mpace_a_init

  private ! Default Scope

    
  ! These variables were moved here so that they could be 
  ! accessible to all subroutines in mpace_a
  ! Joshua Fasching December 2007

  integer, parameter :: file_ntimes = 139
  integer, parameter :: file_nlevels = 38
  integer, parameter :: per_line = 5
  
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
  subroutine mpace_a_tndcy( time, time_initial, rlat, & 
                            rho, p_in_Pa, rcm, & 
                            Ncnm, Ncm, &
                            wm_zt, wm_zm, thlm_forcing, rtm_forcing, & 
                            Frad, radht, & 
                            um_hoc_grid, vm_hoc_grid, & 
                            sclrm_forcing, edsclrm_forcing )

!        Description:
!
!        References:
!          Liou, Wallace and Hobbs, Shettle and Weinman
!-----------------------------------------------------------------------

  use constants, only: Cp, Rd, Lv, p0, rc_tol, zero_threshold ! Variable(s)

  use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

  use model_flags, only: l_bugsrad, l_coamps_micro, l_kk_rain ! Variable(s)

  use grid_class, only: gr ! Variable(s)

  use grid_class, only: zt2zm ! Procedure(s)

  use interpolation, only: zlinterp_fnc, factor_interp ! Procedure(s)

  use stats_precision, only: time_precision ! Variable(s)

  use rad_lwsw_mod, only: rad_lwsw ! Procedure(s)

  use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

  use error_code, only: clubb_debug ! Procedure(s)
 
  use stats_type, only: stat_update_var ! Procedure(s)

  use stats_variables, only: iradht_LW, iradht_SW, iFrad_LW,  & ! Variable(s)
                 iFrad_SW, zt, zm, l_stats_samp
 
  implicit none

  ! Local constants, subsidence
  real, parameter :: & 
!     .  grav0 = 9.8,     ! m/s
!     .  D     = 5.8e-6,  ! 1/s
  psfc  = 101000.   ! Pa
!     .  pinv  = 85000.   ! Pa; ditto

  ! Local constants, LW radiation (from DYCOMS II-RF01)
  real, parameter :: & 
  F0  = 100.0, & 
  F1  = 20.0, & 
  kap = 119.0

  ! Local constants, SW radiation (Shettle and Weinman)
  real, parameter :: & 
  Fs0    = 1212.75, & 
  radius = 1.0e-5, & 
  A      = 0.1, & 
  gc     = 0.86, & 
  omega  = 0.9965
!    .  rlat = 71.75

  ! Local constants, SW radiation (Liou solar angle scheme)
  real, parameter :: & 
  c_0 = 0.006918, & 
  c_1 = -0.399912, & 
  c_2 = -0.006758, & 
  c_3 = -0.002697, & 
  d_1 = 0.070257, & 
  d_2 = 0.000907, & 
  d_3 = 0.000148

  ! Input Variables
  real(kind=time_precision), intent(in) ::  & 
  time,         & ! Current time of simulation      [s]
  time_initial    ! Initial time of simulation      [s]

  real, intent(in) ::  & 
  rlat          ! Latitude                        [Degrees North]

  real, dimension(gr%nnzp), intent(in) :: & 
  rho,     & ! Density of air                         [kg/m^3]
  p_in_Pa, & ! Pressure                               [Pa]
  rcm        ! Cloud water mixing ratio               [kg/kg]

  ! Input/Output Variables
  real, dimension(gr%nnzp), intent(inout) ::  & 
  Ncm,     & ! Cloud droplet number concentration      [count/m^3]
  Ncnm       ! Cloud nuclei number concentration       [count/m^3]

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
  edsclrm_forcing ! Eddy-passive scalar forcing         [units/s]

  ! Local Variables, radiation scheme
  real, dimension(gr%nnzp) ::  & 
  radht_LW, & ! dT/dt, then d Theta/dt, due to LW rad.  [K/s]
  radht_SW, & ! dT/dt, then d Theta/dt, due to SW rad.  [K/s]
  Frad_LW,  & ! Longwave radiative flux                 [W/m^2]
  Frad_SW     ! Shortwave radiative flux                [W/m^2]


  ! Local Variables, general
  integer :: i, k ! Loop indices


  ! Local Variables, subsidence scheme
!        real ::
!     .  velocity_omega


  ! Local Variables, radiation scheme
  real :: & 
  xi_abs, & 
  sda_t, & 
  sda_delta, & 
  sda_h, & 
  t_since_noon, & 
  julday, & 
  start_time_until_noon

  real, dimension(gr%nnzp) :: & 
  radht_theta, & 
  radht_LW_theta, & 
  radht_SW_theta,  & 
!     .  LWP,            ! Liquid water path                              [kg/m^2]
  rcm_rad,         & ! Flipped array of liq. water mixing ratio       [kg/kg]
  rho_rad,         & ! Flipped array of air density                   [kg/m^3]
  dsigm,           & ! Flipped array of grid spacing                  [m]
  coamps_zm,       & ! Flipped array of momentum level altitudes      [m]
  coamps_zt,       & ! Flipped array of thermodynamic level altitudes [m]
  frad_out,        & ! Flipped array of radiaive flux                 [W/m^2]
  frad_lw_out,     & ! Flipped array of LW radiative flux             [W/m^2]
  frad_sw_out,     & ! Flipped array of SW radiative flux             [W/m^2]
  radhtk,          & ! Flipped array of radiative heating             [K/s]
  radht_lw_out,    & ! Flipped array of LW radiative heating          [K/s]
  radht_sw_out       ! Flipped array of SW radiative heating          [K/s]

  ! Local variables, on/off switches for individual schemes
  logical ::  & 
  l_lw_on, & 
  l_sw_on, & 
!     .  l_subs_on,
  l_center

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

  real, dimension(gr%nnzp), intent(out) ::  & 
  um_hoc_grid,       & ! Observed wind, for nudging         [m/s]
  vm_hoc_grid       ! Observed wind, for nudging         [m/s]

! This code block takes the model time, finds the time before it and the time after it on 
! the list, and marks them left_time and right_time for interpolation.  If the time is 
! before the first or after the last time in the file, it just uses the first or last 
! time without interpolation.

left_time = -1
right_time = -1

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
    if ((time > file_times(k)) .AND. & 
        (time <=file_times(k+1))) then
      left_time = k
      right_time = k+1
    end if
  end do
end if

if( left_time == -1 .or. right_time == -1 ) then
  call clubb_debug(1, "file_times not sorted in mpace_a_tndcy.")
endif

! This is the ratio "a" needed for linear interpolation in time.
ratio = real((time - file_times(left_time)) /  &          ! at the first time a=0;
        (file_times(right_time) - file_times(left_time))) ! at the second time a=1.

do k=1,file_nlevels
!        omega_column(k) = ratio *			       ! Do linear interpolation in time
!     .                      (omega_forcing(k,right_time)
!     .                      -omega_forcing(k,left_time))
!     .                     + omega_forcing(k,left_time)

  dTdt_column(k)  = factor_interp( ratio, dTdt_forcing(k, right_time), dTdt_forcing(k, left_time) )
  dqdt_column(k)  = factor_interp( ratio, dqdt_forcing(k, right_time), dqdt_forcing(k, left_time) )
  vertT_column(k) = factor_interp( ratio, vertT_forcing(k,right_time),vertT_forcing(k,left_time) )
  vertq_column(k) = factor_interp( ratio, vertq_forcing(k,right_time),vertq_forcing(k,left_time) )
  um_column(k)    = factor_interp( ratio, um_obs(k, right_time), um_obs(k, left_time) )
  vm_column(k)    = factor_interp( ratio, vm_obs(k, right_time), vm_obs(k, left_time) )
end do

!     Do linear interpolation in space
!     using zlinterp_fnc
dTdt_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt, & 
                         file_heights,dTdt_column)
dqdt_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt, & 
                         file_heights,dqdt_column)
vertT_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt, & 
                         file_heights,vertT_column)
vertq_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt, & 
                         file_heights,vertq_column)
um_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt, & 
                         file_heights,um_column)
vm_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt, & 
                         file_heights,vm_column)

um_hoc_grid (1) = um_hoc_grid(2)
vm_hoc_grid (1) = vm_hoc_grid(2)

! eMFc

!-----------------------------------------------------------------------

  ! Set which schemes to use
  l_lw_on           = .TRUE.
  l_sw_on           = .TRUE.
!        l_subs_on         = .TRUE.
  l_center          = .TRUE.

  ! Compute vertical motion
  do i=2,gr%nnzp
!          velocity_omega = omega_hoc_grid(i) * 100 / 3600 ! convering mb/hr to Pa/s
!          wm_zt(i) = -velocity_omega * Rd * thvm(i) / p_in_Pa(i) / grav
     wm_zt(i) = 0.
! End of Michael Falk's obliteration of omega.
  end do

  ! Boundary condition
  wm_zt(1) = 0.0        ! Below surface

  ! Interpolation
  wm_zm = zt2zm( wm_zt )

  ! Boundary condition
  wm_zm(1) = 0.0        ! At surface
  wm_zm(gr%nnzp) = 0.0  ! Model top
  

  ! Compute large-scale tendencies
  do i=1,gr%nnzp
   thlm_forcing(i) = ((dTdt_hoc_grid(i) + vertT_hoc_grid(i)) & 
                    * ((psfc/p_in_Pa(i)) ** (Rd/Cp))) & 
                    / 3600. ! K/s
   rtm_forcing(i)  = (dqdt_hoc_grid(i)+vertq_hoc_grid(i)) & 
    / 1000. / 3600. ! g/kg/hr -> kg/kg/s
  end do

  ! Compute radiation
  julday = 282
  start_time_until_noon = 29051 + 50400
  t_since_noon   = real( time - start_time_until_noon )
  sda_t = 2*3.14*(julday-1)/365

  sda_delta = c_0 + c_1*cos(sda_t) + d_1*sin(sda_t) + & 
              c_2*cos(2*sda_t) + d_2*sin(2*sda_t) + & 
              c_3*cos(3*sda_t) + d_3*sin(3*sda_t)

  sda_h = 2*3.14*t_since_noon/86400

  xi_abs = sin(rlat*3.14/180) * sin(sda_delta) + & 
           cos(rlat*3.14/180) * cos(sda_delta) * cos(sda_h)

  xi_abs = max(xi_abs,zero_threshold)

  if (xi_abs == 0.) then
    l_sw_on = .FALSE.
  end if

  if (.not. l_sw_on) then
    xi_abs = 0.
  end if

  if ( .not. l_bugsrad ) then
    do k = 1, gr%nnzp
      rcm_rad(k)  = rcm(gr%nnzp-k+1)
      rho_rad(k) = rho(gr%nnzp-k+1)
      dsigm(k)    = 1.0 / gr%dzt(gr%nnzp-k+1)
      coamps_zm(k) = gr%zm(gr%nnzp-k+1)
      coamps_zt(k) = gr%zt(gr%nnzp-k+1)
    enddo

    call rad_lwsw( rcm_rad, rho_rad, dsigm, & 
                   coamps_zm, coamps_zt, & 
                   Frad_out, Frad_LW_out, Frad_SW_out, & 
                   radhtk, radht_LW_out, radht_SW_out, & 
                   gr%nnzp-1, l_center, & 
                   xi_abs, F0, F1, kap, radius, A, gc, Fs0, omega, & 
                   l_sw_on, l_lw_on )

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
      thlm_forcing(k) = thlm_forcing(k) + radht_theta(k)
    end do

  end if ! ~ l_bugsrad

  if ( .not.l_bugsrad .and. l_stats_samp ) then
 
    call stat_update_var( iradht_LW, radht_LW, zt )

    call stat_update_var( iradht_SW, radht_SW, zt )

    call stat_update_var( iFrad_SW, Frad_SW, zm )

    call stat_update_var( iFrad_LW, Frad_LW, zm )

  end if
 

  ! Initialize Ncnm on first timestep
  if ( l_coamps_micro .and. time == time_initial ) then
    Ncnm(1:gr%nnzp) = 30.0 * (1.0 + exp(-gr%zt(1:gr%nnzp)/2000.0)) * 1.e6

  else if ( l_kk_rain ) then
    ! Note: Khairoutdinov and Kogan microphysics has only been
    ! tested for marine stratocumulous clouds, and does not
    ! account for snow and ice.
    do k=1, gr%nnzp, 1
      if ( rcm(k) >= rc_tol ) then
        Ncm(k) = 30.0 * (1.0 + exp(-gr%zt(k)/2000.0)) * 1.e6 & 
                 / rho(k) 
      end if
    end do
  end if


  ! Test scalars with thetal and rt if desired
  if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
  if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

  if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
  if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

  return
  end subroutine mpace_a_tndcy

!----------------------------------------------------------------------
  subroutine mpace_a_sfclyr( time, rho0, um_sfc, vm_sfc, & 
                             upwp_sfc, vpwp_sfc, & 
                             wpthlp_sfc, wprtp_sfc, ustar, & 
                             wpsclrp_sfc, wpedsclrp_sfc )
!        Description:
!          Surface forcing subroutine for mpace_a case.  Written 
!          October 2007 by Michael Falk.
!
!        References:
!          mpace_a specification from arm.gov
!-----------------------------------------------------------------------

  use constants, only: Cp, Lv ! Variable(s)

  use parameters_model, only: sclr_dim, edsclr_dim  ! Variable(s)

  use stats_precision, only: time_precision ! Variable(s)

  use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

  use error_code, only: clubb_debug ! Procedure(s)

  use interpolation, only: factor_interp ! Procedure(s)

  implicit none

  ! External
  intrinsic :: max, sqrt, present

  ! Parameter Constants
  real, parameter :: & 
  ubmin = 0.25
!     .  ustar = 0.25

  ! Input Variables
  real(kind=time_precision), intent(in) :: & 
  time     ! current model time           [s]

  real, intent(in)  :: & 
  rho0,     & ! Air density at surface       [kg/m^3]
  um_sfc,   & ! um at zt(2)                  [m/s]
  vm_sfc      ! vm at zt(2)                  [m/s]

  ! Output Variables
  real, intent(out) ::  & 
  upwp_sfc,     & ! u'w' at (1)      [m^2/s^2]
  vpwp_sfc,     & ! v'w'at (1)       [m^2/s^2]
  wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
  wprtp_sfc,    & ! w'r_t' at (1)    [(m kg)/(s kg)]
  ustar           ! surface friction velocity [m/s]

  real, dimension(sclr_dim), intent(out) :: & 
  wpsclrp_sfc    ! Passive scalar surface flux      [units m/s]

  real, dimension(edsclr_dim), intent(out) :: & 
  wpedsclrp_sfc  ! Passive eddy-scalar surface flux [units m/s]

  ! Local Variables
   real :: & 
   ubar, & 
   latent_heat_flx, & 
   sensible_heat_flx

   integer :: k

   integer :: & 
   left_time, right_time

   real :: ratio
!-----------------------------------------------------------------------

left_time = -1
right_time = -1

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
    if ((time > file_times(k)) .AND. & 
        (time <= file_times(k+1))) then
      left_time = k
      right_time = k+1
    end if
  end do
end if

! Sanity check to make certain that the values read into
! file_times are sorted. Joshua Fasching June 2008
if ( left_time == -1 .or. right_time == -1 ) then
        call clubb_debug(1, "file_times not sorted in MPACE_A")
endif

ratio = real(((time-file_times(left_time)) /  & 
     (file_times(right_time)-file_times(left_time))))

latent_heat_flx = factor_interp( ratio, file_LH(right_time), file_LH(left_time) )

sensible_heat_flx = factor_interp( ratio, file_SH(right_time), file_SH(left_time) )

 ! Compute heat and moisture fluxes
  wpthlp_sfc = sensible_heat_flx/(rho0*Cp)
  wprtp_sfc  = latent_heat_flx/(rho0*Lv)

  ! Compute momentum fluxes
  ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

  ! Declare the value of ustar.
  ustar = 0.25

  upwp_sfc = -um_sfc * ustar*ustar / ubar
  vpwp_sfc = -vm_sfc * ustar*ustar / ubar

  ! Let passive scalars be equal to rt and theta_l for now
  if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
  if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

  if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
  if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc

  return
  end subroutine mpace_a_sfclyr
!----------------------------------------------------------------
  subroutine mpace_a_init( iunit, file_path )
!
!       Description:
!       This subroutine initializes the module by reading in forcing
!       data used in the tndcy and sfclyr subroutines.
!----------------------------------------------------------------
    use file_functions, only: file_read_1d, file_read_2d ! Procedure(s)

    implicit none

    integer, intent(in) :: iunit ! File unit number

    character(len=*), intent(in) :: &
      file_path ! Path to the forcing files

    ! ---- Begin Code ----

    call file_read_1d( iunit, & 
      file_path//'mpace_a_press.dat', & 
      file_nlevels, per_line, file_pressure )

    call file_read_1d( iunit, & 
      file_path//'mpace_a_heights.dat', & 
      file_nlevels, per_line, file_heights )

    call file_read_1d( iunit, & 
      file_path//'mpace_a_times.dat', & 
      file_ntimes, per_line, file_times )

!      call file_read_2d( iunit,
!     . file_path//'mpace_a_omega.dat',
!     . file_nlevels, file_ntimes, per_line, omega_forcing)

    call file_read_2d( iunit, & 
      file_path//'mpace_a_dTdt.dat', & 
      file_nlevels, file_ntimes, per_line, dTdt_forcing )

    call file_read_2d( iunit, & 
      file_path//'mpace_a_dqdt_horiz.dat', & 
      file_nlevels, file_ntimes, per_line, dqdt_forcing )

    call file_read_2d( iunit, & 
      file_path//'mpace_a_verts.dat', & 
      file_nlevels, file_ntimes, per_line, vertT_forcing )

    call file_read_2d( iunit, & 
      file_path//'mpace_a_vertq.dat', & 
      file_nlevels, file_ntimes, per_line, vertq_forcing )

    call file_read_2d( iunit, & 
      file_path//'mpace_a_um_obs.dat', & 
      file_nlevels, file_ntimes, per_line, um_obs )

    call file_read_2d( iunit, & 
      file_path//'mpace_a_vm_obs.dat', & 
      file_nlevels, file_ntimes, per_line, vm_obs )

    call file_read_1d( iunit, & 
      file_path//'mpace_a_lh.dat', & 
      file_ntimes, per_line, file_LH )

    call file_read_1d( iunit, & 
      file_path//'mpace_a_sh.dat', & 
      file_ntimes, per_line, file_SH )

    return 
  end subroutine mpace_a_init

  end module mpace_a
