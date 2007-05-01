#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: jun25.F,v 1.2 2007-05-01 22:57:30 dschanen Exp $
        module jun25

!       Description:
!       Contains subroutines for the June 11 case.
!----------------------------------------------------------------------

        public :: jun25_altocu_tndcy

        private

        contains

!-----------------------------------------------------------------------
        subroutine jun25_altocu_tndcy
     .             ( time, time_initial, rlat, dt, thlm, 
     .               rcm, p, exner, rhot, rtm, wmt,
     .               wmm, thlm_forcing, rtm_forcing, 
     .               Frad, radht, 
     .               sclr_tndcy ) ! Optional
!       Description:
!       Computes subsidence, radiation, and LS tendencies for the June
!       25th Altocumulous case.

!       References:
!-----------------------------------------------------------------------
        use grid_class
        use constants
        use model_flags
        use icedfs
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none

        ! Constant parameters

        ! Input variables
        real, intent(in) :: 
     .  time,         ! Time of simulation since start        [s]
     .  time_initial, ! Initial time of simulation            [s]
     .  rlat,         ! Reference latitude should be 37.6     [Degrees North]
     .  dt            ! Length of model timestep              [s]

        real, dimension(gr%nnzp), intent(in) :: 
     .  thlm,  ! Liquid potential temperature           [K]
     .  rcm,   ! Liquid water mixing ratio              [kg/kg]
     .  p,     ! Basic state pressure                   [Pa]
     .  exner, ! Exner function                         [-]
     .  rhot   ! Density of reference state on t grid   [kg/m^3]

        ! Input/Output Variables
        real, dimension(gr%nnzp), intent(inout) :: 
     .  rtm    ! Total water mixing ratio               [kg/kg]

        ! Output variables
        real, dimension(gr%nnzp), intent(inout) :: 
     .  wmt,         ! Vertical ascent/descent on therm. grid      [m/s]
     .  wmm,         ! Vertical ascent/descent on moment. grid     [m/s]
     .  thlm_forcing ! Change in liq. water potential temperature 
                     ! due to radiative heating and ice diffusion  [K/s]

        real, dimension(gr%nnzp), intent(inout) :: 
     .  rtm_forcing, ! Change in total water due to ice diffusion  [kg/kg/s]
     .  Frad,        ! Total radiative flux (LW + SW)              [W/m^2]
     .  radht        ! Total radiative heating (LW +SW)            [K/s]

        ! Output variables (optional)
        real, optional, dimension(gr%nnzp,sclr_dim),intent(out) :: 
     .  sclr_tndcy ! Large-scale tendency for passive scalars

!-----------------------------------------------------------------------
! LOCAL VARIABLES
!
!
! TEMPORARY ARRAYS USED FOR COAMPS RADIATIVE SCHEME
! (SEE COMMENTS BELOW FOR SCHEME DESCRIPTION)
! 
! frad_out      : temporary "flipped" array of total radiative flux
!                 from rad_lwsw                                 Unit: W/m^2
! frad_lw_out   : temporary "flipped" array of LW radiative flux
!                 from rad_lwsw                                 Unit: W/m^2
! frad_sw_out   : temporary "flipped" array of SW radiative flux
!                 from rad_lwsw                                 Unit: W/m^2
!
! radhtk        : temporary "flipped" array of total radiative heating
!                 from rad_lwsw                                 Unit: K/s
! radht_lw_out  : temporary "flipped" array of LW radiative heating
!                 from rad_lwsw                                 Unit: K/s
! radht_sw_out  : temporary "flipped" array of SW radiative heating
!                 from rad_lwsw                                 Unit: K/s
!
! INTERPOLATION ARRAYS AND CONSTANTS
! zsubs         : heights at which wmt data is supplied
!                 (used for subsidence interpolation)           Unit: m
! tsubs         : times after initialization at which wmt data is supplied
!                 (NOT USED IN NOV.11 CASE)                     Unit: s
! wtX(Y)        : vertical velocity specified at height Y and time X
!                 (ONLY wt1 IS USED IN NOV.11 CASE)             Unit: m/s
! w1-w2         : vertical velocity before (w1) and after (w2)
!                 the current time at the specified level
!                 (NOT USED IN NOV.11 CASE)                     Unit: m/s
!
! ADDITIONAL PARAMETERS FOR NOV.11 SUBSIDENCE (NOT FOR JUN.25 CASE)
! subs_on       : logical variable tells us whether to turn subsidence on
!                                                               Unit: NONE
! wmax          : defines value of maximum subsidence in profile
!                                                               Unit: cm/s
! zi            : defines approx. height of inversion within cloud
!                 (subsidence is equal to wmax at this height)  Unit: m
! dazi          : defines height above inversion
!                 (above this height, subsidence linearly tapers off to zero)
!                                                               Unit: m
! dbzi          : defines height above inversion
!                 (below this height, subsidence linearly tapers off to zero)
!                                                               Unit: m
! dac           : defines height above cloud
!                 (at / above this height, we have NO subsidence)
!                                                               Unit: m
! dbc           : defines height below cloud
!                 (at / below this height, we have NO subsidence)
!                                                               Unit: m
!
! RADIATION PARAMETERS
! sw_on         : logical variable passed to radiation scheme
!                 - is SW radiation on?                         Unit: NONE
! lw_on         : logical variable passed to radiation scheme
!                 - is LW radiation on?                         Unit: NONE
! center        : use centered differencing (as opposed to a one-sided
!                 forward difference) in radiation code         Unit: NONE
! 
! xi_abs        : cosine of the solar zenith angle              Unit: NONE
! F0            : coefficient for cloud top heating (see Stevens)
!                                                               Unit: W/m^2
! F1            : coefficient for cloud base heating (see Stevens)
!                                                               Unit: W/m^2
! kap           : "a constant" according to Duynkerke eqn. 5, where his
!                 value is 130 m^2/kg.                          Unit: m^2/kg
! radius        : effective droplet radius                      Unit: m
! AA            : albedo -- sea surface, according to Lenderink.
!                                                               Unit: NONE
! gc            : asymmetry parameter, "g" in Duynkerke.        Unit: NONE
! Fs0           : The incident incoming SW insolation at cloud top in the 
!                 direction of the incoming beam (not the vertical).
!                                                               Unit: W/m^2
! omega         : single-scattering albedo                      Unit: NONE
! 
! SOLAR ZENITH ANGLE PARAMETERS (NOT USED IN NOV.11 CASE)
! c0            : coefficient for calculation of declination angle from
!                 Liou Table 2.2 and Eqn. 2.2.10                Unit: NONE
! c1            : same as above                                 Unit: NONE
! c2            : same as above                                 Unit: NONE
! c3            : same as above                                 Unit: NONE
! d1            : Same as above                                 Unit: NONE
! d2            : Same as above                                 Unit: NONE
! d3            : Same as above                                 Unit: NONE
!
! sda_t         : Linear function of day of the year.
!                 sda_t=0 January 1 and sda_t -> 2*pi December 31.
! sda_delta     : Solar declination angle function from Liou 2.2.10 
! sda_h         : Hour angle (Angle through which the earth must to turn
!                 to put sun directly overhead on a point's meridian)
!                 (Angle between current time and solar noon)
! t_since_noon  : Number of seconds since noon (after noon > 0)
! julday        : Julian day of the year (January 1=1; December 31=365)
! 
! start_time_until_noon: number of seconds between start time and solar noon
!                                                               Unit: s
! 
! Fs0 INTERPOLATION PARAMETERS 
! nparam        : Number of Fs0 values in the list.
! xilist        : Values of cosine of solar zenith angle corresponding to
!                 the values in Fslist
! Fslist        : Values of Fs0 corresponding to the values in xilist.
!-----------------------------------------------------------------------

        !-----------------------------
        ! Local ice diffusional arrays
        !-----------------------------
        real, dimension(gr%nnzp) :: 
     .  dqc_dt_icedfs ! Change in liquid water over time, 
                      ! due to ice diffusional growth [1/s]

        !------------------------
        ! Local radiation arrays
        !------------------------
        real, dimension(gr%nnzp) :: 
     .  Frad_LW, ! Long wave radiative flux     [W/m^2]
     .  Frad_SW, ! Short wave radiative flux    [W/m^2]
     .  radht_LW,! Long wave radiative heating  [K/s]
     .  radht_SW ! Short wave radiative heating [K/s]

        real, dimension(gr%nnzp) :: 
     .  LWP,      ! Liquid water path from domain top            [kg/m^2]
     .  thm_rad,  ! "flipped" array of potential temperature     [K]
     .  rcm_rad,  ! "flipped" array of liquid water mixing ratio [kg/kg]
     .  rhot_rad, ! "flipped" array of air density               [kg/m^3]
     .  dsigm     ! "flipped" array of grid spacing              [m]


        real, dimension(gr%nnzp) :: 
     .  frad_out,
     .  frad_lw_out,
     .  frad_sw_out

        real, dimension(gr%nnzp) :: 
     .  radhtk,
     .  radht_lw_out,
     .  radht_sw_out

        !---------------------------------------------------------------
        ! Working arrays for subsidence interpolation
        !---------------------------------------------------------------
        real, dimension(5) :: 
     .  zsubs ! [m]

        real, dimension(6) :: tsubs

        real, dimension(5) ::
     .  wt1, wt2, wt3, wt4, wt5, wt6

        real, dimension(gr%nnzp) ::
     .  w1, w2

        !---------------------------------------------------------------
        ! LW Radiative constants
        !---------------------------------------------------------------
        real, parameter :: 
     .  F0   = 107.0, ! [W/m^2]
     .  F1   = 61.0,  ! [W/m^2]
     .  kap  = 100.0  ! [m^2/kg]

        !---------------------------------------------------------------
        ! Working arrays for SW radiation interpolation
        !---------------------------------------------------------------
        integer, parameter :: nparam = 12

        real, dimension(nparam) :: xilist, Fslist

        !---------------------------------------------------------------
        ! SW Radiative constants
        !---------------------------------------------------------------
        real, parameter :: 
     .  radius = 1.0e-5,
     .  AA     = 0.1,
     .  gc     = 0.85,
     .  omega  = 0.992

        !---------------------------------------------------------------
        ! Additional SW radiative variables
        !---------------------------------------------------------------
        real :: xi_abs, Fs0

        !---------------------------------------------------------------
        ! Variables used for calculating solar zenith angle
        ! (Currently not used in Nov.11 case)
        !
        ! NOTE BY ADAM SMITH, 13 April 2006
        ! In COAMPS, these constants are called c0, c1, c2, c3, d1, d2,
        ! and d3.  Since HOC uses these variables names elsewhere, I have
        ! renamed them here.
        !---------------------------------------------------------------
        real :: c_0, c_1, c_2, c_3, d_1, d_2, d_3,
     .  sda_t, sda_delta, sda_h, t_since_noon, julday,
     .  start_time_until_noon

        !---------------------------------------------------------------
        ! Toggle for implementing differencing method in interpolations
        !---------------------------------------------------------------
        logical :: center

        !---------------------------------------------------------------
        ! Toggles for activating/deactivating forcings
        !---------------------------------------------------------------
        logical :: subs_on, lw_on, sw_on

        !---------------------------------------------------------------
        ! Variable used for working within vertical arrays
        !---------------------------------------------------------------
        integer :: k

        !---------------------------------------------------------------
        ! END OF VARIABLE DECLARATION
        !---------------------------------------------------------------


!-----------------------------------------------------------------------
! Toggles for activating/deactivating forcings
! To turn off a specific forcing, set the corresponding toggle to .FALSE.
!-----------------------------------------------------------------------
         subs_on   = .TRUE.
         lw_on     = .TRUE.
         sw_on     = .TRUE.

!-----------------------------------------------------------------------
! Toggle for centered/forward differencing (in interpolations)
! To use centered differencing, set the toggle to .TRUE.
! To use forward differencing, set the toggle to .FALSE.
!-----------------------------------------------------------------------
         center    = .TRUE.

!-----------------------------------------------------------------------
! Calculation of Solar Declination Angle
!
! This section of code calculates the solar declination angle based
! on the simulaiton's time, date, and latitude.  These parameters are
! applied in the scheme, and a value of xi_abs is calculated.
! 
! If xi_abs <= 0.0 at any time step, we shut off shortwave radiation
! for that time step (by setting sw_on = .FALSE.).  If xi_abs > 0.0,
! we reactivate solar radiation by setting sw_on = .TRUE.
!
! This version is based on solar_declin.m, written by Dr. Vince Larson
! and based on Liou 2.2.10.
!-----------------------------------------------------------------------

        julday = 176
        start_time_until_noon = 30540
        t_since_noon   = time - start_time_until_noon
        sda_t = 2.0*pi*(julday-1)/365

        c_0 = 0.006918
        c_1 = -0.399912
        c_2 = -0.006758
        c_3 = -0.002697
        d_1 = 0.070257
        d_2 = 0.000907
        d_3 = 0.000148

        sda_delta = c_0 + c_1*cos(sda_t)   + d_1*sin(sda_t)   +
     .                    c_2*cos(2*sda_t) + d_2*sin(2*sda_t) +
     .                    c_3*cos(3*sda_t) + d_3*sin(3*sda_t)

        sda_h = 2.0*pi*t_since_noon/86400.0

        xi_abs = sin(rlat*pi/180) * sin(sda_delta) +
     .           cos(rlat*pi/180) * cos(sda_delta) * cos(sda_h)

        xi_abs = max(xi_abs,0.)

!-----------------------------------------------------------------------
! Modification by Adam Smith 26 June 2006
! It is difficult to remember to set xi_abs = 0 when we want to shut off
! solar radiation.  If sw_on = .FALSE. above, we will automatically set
! xi_abs to 0 to avoid confusion or errors.
!-----------------------------------------------------------------------
      if ( .not. sw_on ) then
        xi_abs = 0.
      end if

!-----------------------------------
! End of ajsmith4's Modification
!-----------------------------------

      if (xi_abs == 0.) then
        sw_on = .FALSE.
      else
        sw_on = .TRUE.
      end if


!-----------------------------------------------------------------------
!                                                                      c
! Fs0 Interpolation Parameters-- these also from Kurt Kotenberg's      c
! BUGSrad output.  Fs0 changes somewhat over the range of solar zenith c
! angles, and we obtained these values by solving                      c
! Fs0 = F_vertical / xi_abs .                                          c
!                                                                      c
! The linear_interpolation function returns Fs0.                       c
!-----------------------------------------------------------------------

      xilist(1) = 0.0
      xilist(2) = 0.01
      xilist(3) = 0.1
      xilist(4) = 0.2
      xilist(5) = 0.3
      xilist(6) = 0.4
      xilist(7) = 0.5
      xilist(8) = 0.6
      xilist(9) = 0.7
      xilist(10) = 0.8
      xilist(11) = 0.9
      xilist(12) = 1.0

      Fslist(1)  = 0.0
      Fslist(2)  = 715.86
      Fslist(3)  = 1073.577
      Fslist(4)  = 1165.0905
      Fslist(5)  = 1204.7033
      Fslist(6)  = 1227.6898
      Fslist(7)  = 1243.1772
      Fslist(8)  = 1254.5893
      Fslist(9)  = 1263.5491
      Fslist(10) = 1270.8668
      Fslist(11) = 1277.0474
      Fslist(12) = 1282.3994

      call linear_interpolation( nparam, xilist, Fslist, xi_abs, Fs0 )

!-----------------------------------------------------------------------
! Subsidence Parameters
!-----------------------------------------------------------------------

        ! Modification for setting June 25th gr%zm(1) to the correct
        ! value in meters for the actual altitude, rather than 0.
        ! -dschanen 1 May 2007
!       zsubs(1) = 0
!       zsubs(2) = 360
!       zsubs(3) = 1090
!       zsubs(4) = 1890
!       zsubs(5) = 2500

        zsubs(1) = gr%zm(1)
        zsubs(2) = gr%zm(1) + 360.
        zsubs(3) = gr%zm(1) + 1090.
        zsubs(4) = gr%zm(1) + 1890.
        zsubs(5) = gr%zm(1) + 2500.

        tsubs(1) = 0
        tsubs(2) = 10800
        tsubs(3) = 28800
        tsubs(4) = 36000
        tsubs(5) = 36000
        tsubs(6) = 36000

        wt1(1) = 0.
        wt1(2) = .004
        wt1(3) = .004
        wt1(4) = .004
        wt1(5) = 0.

        wt2(1) = 0.
        wt2(2) = .004
        wt2(3) = .004
        wt2(4) = .004
        wt2(5) = 0.

        wt3(1) = 0.
        wt3(2) = -.003
        wt3(3) = -.003
        wt3(4) = -.003
        wt3(5) = 0.

        wt4(1) = 0.
        wt4(2) = -.003
        wt4(3) = -.003
        wt4(4) = -.003
        wt4(5) = 0.

        wt5(1) = 0.
        wt5(2) = -.003
        wt5(3) = -.003
        wt5(4) = -.003
        wt5(5) = 0.

        wt6(1) = 0.
        wt6(2) = -.003
        wt6(3) = -.003
        wt6(4) = -.003
        wt6(5) = 0.


!-----------------------------------------------------------------------
! SPECIAL METHOD USED TO CALCULATE RADIATION                           c
! Grid descriptions by Adam Smith, 27 June 2006                        c
!                                                                      c
! In order to verify our HOC simulations are working properly, we      c
! have first developed a series of 3D simulations using the COAMPS-LES c
! model.  This large-eddy simulation (LES) simulation uses specific    c
! methods to calculate radiation, subsidence, and other microphysical  c
! processes.  To make the two models simluate clouds as closely as     c
! possible, we use the same radiation scheme in both models.           c
!                                                                      c
! In COAMPS-LES, we use a separate subroutine, rad_lwsw, to implement  c
! all radiation code.  This allows the subroutine to be duplicated     c
! exactly in many different models.  However, the subroutine uses the  c
! COAMPS vertical grid.  Therefore, for HOC to implement this code     c
! correctly, we must modify some of our variable profiles before       c
! calling the radiation subroutine.                                    c
!                                                                      c
! The following diagram describes the differences in model grids:      c
!                                                                      c
!       COAMPS-LES                                    HOC              c
!                                                                      c
!  k= 1   (m) ----------    <MODEL TOP>    ---------- k=kk+1      (m)  c
!                                                                      c
!  k= 1   (t) ----------                   ---------- k=kk+1      (t)  c
!                                                                      c
!  k= 2   (m) ----------                   ---------- k=kk        (m)  c
!                                                                      c
!  k= 2   (t) ----------                   ---------- k=kk        (t)  c
!                                                                      c
!            .                  .                     .                c
!            .                  .                     .                c
!            .                  .                     .                c
!                                                                      c
!  k=kk-1 (m) ----------  m = mass levels  ---------- k=3         (m)  c
!                                                                      c
!  k=kk-1 (t) ----------  t = thermo       ---------- k=3         (t)  c
!                                  levels                              c
!  k=kk   (m) ----------                   ---------- k=2         (m)  c
!                                                                      c
!  k=kk   (t) ----------  kk = number of   ---------- k=2         (t)  c
!                              vertical                                c
!  k=kk+1 (m) ----------       heights     ---------- k=1         (m)  c
!                                                                      c
! //////////////////////// MODEL SURFACE ///////////////////////////// c
!                                          ---------- k=1         (t)  c
!                                                                      c
!                                                                      c
! The major difference in the grids is that HOC uses an additional     c
! thermodynamic level below the model "surface".  This means that all  c
! HOC thermodynamic heights are shifted down one vertical level, and   c
! HOC also has one fewer momentum level than COAMPS.  Therefore, we    c
! use one additional vertical level in HOC, to make sure that the      c
! vertical domain matches in both models.                              c
!                                                                      c
! Also, the COAMPS grid indices are numbered from the top of the model c
! downward, while the HOC grid indices are numbered from the bottom    c
! up.  Therefore, since we are using a COAMPS radiation scheme, we     c
! flip moisture and temperature profiles that are passed into the      c
! rad_lwsw subroutine.  The rad scheme will produce results in using   c
! the COAMPS grid scheme, so all radiation output will be flipped      c
! back to the HOC grid before being applied to the model.              c
!                                                                      c
! Finally, since the COAMPS scheme does not have a gridpoint below     c
! model surface, we add that point to all radiative output files once  c
! they are converted back to HOC setup.  This allows all averages and  c
! calculations to be done correctly.                                   c
!                                                                      c
!                                                                      c
! Computation of radiative fluxes on staggered grid                    c
! Comments by Michael Falk, 16 February 2005.                          c
!                                                                      c
! Frad (and its components Frad_LW and Frad_SW) should be computed on  c
! w points, not on mass points, which is apparent from its formulation c
! and from its location in stats_sw instead of stats_sm.  The grid     c
! looks like this:                                                     c
!                                                                      c
!                                                                      c
! -----Frad----------------------------------    k = 1  (w level)      c
!     /    \            |-dwm                                          c
! -LWP------radht----------------------------    k = 1  (mass level)   c
!     \    /            |-dmw                                          c
! -----Frad----------------------------------    k = 2  (w level)      c
!     /    \                                                           c
! -LWP------radht----------------------------    k = 2  (mass level)   c
!     \    /                                                           c
! -----Frad----------------------------------    k = 3  (w level)      c
!     /    \                                                           c
! -LWP------radht----------------------------    k = 3  (mass level)   c
!                                                                      c
! If you consider Frad to take place on mass levels, then computing    c
! LWP is a forward difference and is only first-order accurate, while  c
! if Frad computed in between LWP levels, it is a centered difference  c
! which is second-order accurate.                                      c
!                                                                      c
! The coding implementation requires that Frad depend on LWP(k) and    c
! LWP(k-1) since the w level for a given k is at a higher altitude     c
! than the mass level.  radht, back on mass levels, depends on Frad(k) c
! and Frad(k+1).                                                       c
!                                                                      c
! Additionally, these computations assume that the distance between    c
! mass levels (dsigma) is constant, and that the w levels (spaced by   c
! dsigmw) always fall exactly halfway in between the mass levels.  If  c
! this is not the case, consider dwm to be the distance between a w    c
! level and the mass level below it, and dmw to be the distance        c
! between a mass level and the w level below it.  Then, the            c
! formulation for Frad_LW, for instance, would use a weighted average: c
!                                                                      c
! (dwm/(dwm+dmw)) * lwp(k) + (dmw/(dwm+dmw)) * lwp(k-1)                c
! which, for dwm always == dmw, reduces to                             c
! (1/2) * (lwp(k)) + (1/2) * (lwp(k-1))                                c
! which is identical to the current formulation.                       c
! ((lwp(k)+lwp(k-1))/2)                                                c
!                                                                      c
!                                                                      c
!-----------------------------------------------------------------------

        !--------------------------------------------------------------- 
        ! We only implement this section if we choose not to use the
        ! BUGSrad interactive radiation scheme.
        !---------------------------------------------------------------

        if ( .not. lbugsrad ) then

        !----------------------------------------------------------------
        ! This code transforms these profiles from HOC grid to COAMPS
        ! grid.  The COAMPS-grid profiles are then passed to rad_lwsw
        ! for implementation.
        !----------------------------------------------------------------
          do k = 1, gr%nnzp
            thm_rad(k) = thlm(gr%nnzp-k+1) 
     .                 + Lv/(Cp*exner(k))*rcm(gr%nnzp-k+1)
            rcm_rad(k) = rcm(gr%nnzp-k+1)
            rhot_rad(k) = rhot(gr%nnzp-k+1)
            dsigm(k) = 1.0 / gr%dzt(gr%nnzp-k+1)
          end do

        !----------------------------------------------------------------
        ! Calling the radiation subroutine, which uses the COAMPS
        ! grid method.  All input and output profiles use the COAMPS
        ! grid setup.
        !----------------------------------------------------------------
          call rad_lwsw(thm_rad, rcm_rad, rhot_rad, dsigm,
     .                  Frad_out, Frad_LW_out, Frad_SW_out,
     .                  radhtk, radht_LW_out, radht_SW_out,
     .                  gr%nnzp-1, center,
     .                  xi_abs, F0, F1, kap, radius, AA, gc, Fs0, omega,
     .                  sw_on, lw_on)


        !-------------------------------------------------------------
        ! This code transforms the radiation results back into HOC
        ! grid setup.  These Frad and radht arrays are actually
        ! applied to the HOC model.
        !
        ! The radht results are initially calculated in terms of
        ! standard temperature (T).  However, HOC calculates
        ! temperature in terms of potential temperature (theta).
        ! Therefore, we multiply all radht results by (1.0/exner)
        ! to convert from T to theta.
        !-------------------------------------------------------------
          do k = 1, gr%nnzp-1
            Frad(k)     = Frad_out(gr%nnzp-k+1)
            Frad_LW(k)  = Frad_LW_out(gr%nnzp-k+1)
            Frad_SW(k)  = Frad_SW_out(gr%nnzp-k+1)

            radht(k)    = ( 1.0/exner(k) ) * radhtk(gr%nnzp-k+1)
            radht_LW(k) = ( 1.0/exner(k) ) * radht_LW_out(gr%nnzp-k+1)
            radht_SW(k) = ( 1.0/exner(k) ) * radht_SW_out(gr%nnzp-k+1)
          end do

          Frad(1) = Frad(2)
          Frad_LW(1) = Frad_LW(2)
          Frad_SW(1) = Frad_SW(2)

          radht(1) = radht(2)
          radht_LW(1) = radht_LW(2)
          radht_SW(1) = radht_SW(2)

        END IF ! ~lbugsrad


        !-------------------------------------------------------------
        ! Compute the loss of total water due to diffusional
        ! growth of ice.  This is defined on thermodynamic levels.
        !-------------------------------------------------------------
        IF (time >= time_initial + 3600.0 .and. licedfs ) THEN

           CALL altocu_icedf( dt, thlm, rcm, p, exner, rhot,
     .                        dqc_dt_icedfs )

        ELSE

           DO k = 1, gr%nnzp, 1
              dqc_dt_icedfs(k) = 0.0
           ENDDO

        END IF


        !---------------------------------------------------------------
        ! Using linear interpolation to calculate subsidence
        ! Original code by Michael Falk
        ! Added for Jun.25 case by Adam Smith, 13 April 2006
        !---------------------------------------------------------------
        if ( time < tsubs(1) ) then
        do k=1,gr%nnzp
          call linear_interpolation(5,zsubs,wt1,gr%zt(k),wmt(k))
        end do
 
        else if (time < tsubs(2)) then
        do k=2,gr%nnzp
          call linear_interpolation(5,zsubs,wt1,gr%zt(k),w1(k))
          call linear_interpolation(5,zsubs,wt2,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(1))/(tsubs(2)-tsubs(1))*(w2(k)-w1(k))+w1(k)
        end do
 
        else if (time < tsubs(3)) then
        do k=2,gr%nnzp
          call linear_interpolation(5,zsubs,wt2,gr%zt(k),w1(k))
          call linear_interpolation(5,zsubs,wt3,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(2))/(tsubs(3)-tsubs(2))*(w2(k)-w1(k))+w1(k)
        end do
 
        else if (time < tsubs(4)) then
        do k=2,gr%nnzp
          call linear_interpolation(5,zsubs,wt3,gr%zt(k),w1(k))
          call linear_interpolation(5,zsubs,wt4,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(3))/(tsubs(4)-tsubs(3))*(w2(k)-w1(k))+w1(k)
        end do
 
        else if (time < tsubs(5)) then
        do k=2,gr%nnzp
          call linear_interpolation(5,zsubs,wt4,gr%zt(k),w1(k))
          call linear_interpolation(5,zsubs,wt5,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(4))/(tsubs(5)-tsubs(4))*(w2(k)-w1(k))+w1(k)
        end do
 
        else if (time < tsubs(6)) then
        do k=2,gr%nnzp
          call linear_interpolation(5,zsubs,wt5,gr%zt(k),w1(k))
          call linear_interpolation(5,zsubs,wt6,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(5))/(tsubs(6)-tsubs(5))*(w2(k)-w1(k))+w1(k)
        end do
 
        else if (time >= tsubs(6)) then
        do k=2,gr%nnzp
          call linear_interpolation(5,zsubs,wt6,gr%zt(k),wmt(k))
        end do
        end if

        wmt(1) = wmt(2)

        wmm = zt2zm(wmt)


        !---------------------------------------------------------------
        ! Enter the final theta-l and rtm tendencies
        !---------------------------------------------------------------
        DO k = 1, gr%nnzp, 1
           IF ( .not. lbugsrad ) THEN
             thlm_forcing(k) = radht(k)
     .                    - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           ELSE
             thlm_forcing(k) = - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           END IF
           rtm_forcing(k) = dqc_dt_icedfs(k)
        ENDDO

        if ( present( sclr_tndcy ) ) then
          sclr_tndcy(:,SCLR_THETA) = thlm_forcing
          sclr_tndcy(:,SCLR_RT)    = rtm_forcing
        end if 

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

        return
        end subroutine jun25_altocu_tndcy

        end module jun25
