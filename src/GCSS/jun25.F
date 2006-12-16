#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: jun25.F,v 1.1 2006-12-16 00:32:45 dschanen Exp $
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c LOCAL VARIABLES                                                             c
c dqc_dt_icedfs : change in liquid water over time, due to ice diffusional    c
c                 growth.                                       Unit: 1/s     c
c                                                                             c
c LWP           : liquid water path from domain top.            Unit: kg/m^2  c
c                                                                             c
c TEMPORARY ARRAYS USED FOR COAMPS RADIATIVE SCHEME                           c
c (SEE COMMENTS BELOW FOR SCHEME DESCRIPTION)                                 c
c thm_rad       : temporary "flipped" array of potential temperature          c
c                                                               Unit: K       c
c rcm_rad       : temporary "flipped" array of liquid water mixing ratio      c
c                                                               Unit: kg/kg   c
c rhot_rad      : temporary "flipped" array of air density      Unit: kg/m^3  c
c dsigm         : temporary "flipped" array of grid spacing     Unit: m       c
c                                                                             c
c frad_out      : temporary "flipped" array of total radiative flux           c
c                 from rad_lwsw                                 Unit: W/m^2   c
c frad_lw_out   : temporary "flipped" array of LW radiative flux              c
c                 from rad_lwsw                                 Unit: W/m^2   c
c frad_sw_out   : temporary "flipped" array of SW radiative flux              c
c                 from rad_lwsw                                 Unit: W/m^2   c
c                                                                             c
c radhtk        : temporary "flipped" array of total radiative heating        c
c                 from rad_lwsw                                 Unit: K/s     c
c radht_lw_out  : temporary "flipped" array of LW radiative heating           c
c                 from rad_lwsw                                 Unit: K/s     c
c radht_sw_out  : temporary "flipped" array of SW radiative heating           c
c                 from rad_lwsw                                 Unit: K/s     c
c                                                                             c
c INTERPOLATION ARRAYS AND CONSTANTS                                          c
c zsubs         : heights at which wmt data is supplied                       c
c                 (used for subsidence interpolation)           Unit: m       c
c tsubs         : times after initialization at which wmt data is supplied    c
c                 (NOT USED IN NOV.11 CASE)                     Unit: s       c
c wtX(Y)        : vertical velocity specified at height Y and time X          c
c                 (ONLY wt1 IS USED IN NOV.11 CASE)             Unit: m/s     c
c w1-w2         : vertical velocity before (w1) and after (w2)                c
c                 the current time at the specified level                     c
c                 (NOT USED IN NOV.11 CASE)                     Unit: m/s     c
c                                                                             c
c ADDITIONAL PARAMETERS FOR NOV.11 SUBSIDENCE (NOT FOR JUN.25 CASE)           c
c subs_on       : logical variable tells us whether to turn subsidence on     c
c                                                               Unit: NONE    c
c wmax          : defines value of maximum subsidence in profile              c
c                                                               Unit: cm/s    c
c zi            : defines approx. height of inversion within cloud            c
c                 (subsidence is equal to wmax at this height)  Unit: m       c
c dazi          : defines height above inversion                              c
c                 (above this height, subsidence linearly tapers off to zero) c
c                                                               Unit: m       c
c dbzi          : defines height above inversion                              c
c                 (below this height, subsidence linearly tapers off to zero) c
c                                                               Unit: m       c
c dac           : defines height above cloud                                  c
c                 (at / above this height, we have NO subsidence)             c
c                                                               Unit: m       c
c dbc           : defines height below cloud                                  c
c                 (at / below this height, we have NO subsidence)             c
c                                                               Unit: m       c
c                                                                             c
c RADIATION PARAMETERS                                                        c
c sw_on         : logical variable passed to radiation scheme                 c
c                 - is SW radiation on?                         Unit: NONE    c
c lw_on         : logical variable passed to radiation scheme                 c
c                 - is LW radiation on?                         Unit: NONE    c
c center        : use centered differencing (as opposed to a one-sided        c
c                 forward difference) in radiation code         Unit: NONE    c
c                                                                             c
c xi_abs        : cosine of the solar zenith angle              Unit: NONE    c
c F0            : coefficient for cloud top heating (see Stevens)             c
c                                                               Unit: W/m^2   c
c F1            : coefficient for cloud base heating (see Stevens)            c
c                                                               Unit: W/m^2   c
c kap           : "a constant" according to Duynkerke eqn. 5, where his       c
c                 value is 130 m^2/kg.                          Unit: m^2/kg  c
c radius        : effective droplet radius                      Unit: m       c
c AA            : albedo -- sea surface, according to Lenderink.              c
c                                                               Unit: NONE    c
c gc            : asymmetry parameter, "g" in Duynkerke.        Unit: NONE    c
c Fs0           : The incident incoming SW insolation at cloud top in the     c
c                 direction of the incoming beam (not the vertical).          c
c                                                               Unit: W/m^2   c
c omega         : single-scattering albedo                      Unit: NONE    c
c                                                                             c
c SOLAR ZENITH ANGLE PARAMETERS (NOT USED IN NOV.11 CASE)                     c
c c0            : coefficient for calculation of declination angle from       c
c                 Liou Table 2.2 and Eqn. 2.2.10                Unit: NONE    c
c c1            : same as above                                 Unit: NONE    c
c c2            : same as above                                 Unit: NONE    c
c c3            : same as above                                 Unit: NONE    c
c d1            : Same as above                                 Unit: NONE    c
c d2            : Same as above                                 Unit: NONE    c
c d3            : Same as above                                 Unit: NONE    c
c                                                                             c
c sda_t         : Linear function of day of the year.                         c
c                 sda_t=0 January 1 and sda_t -> 2*pi December 31.            c
c sda_delta     : Solar declination angle function from Liou 2.2.10           c
c sda_h         : Hour angle (Angle through which the earth must to turn      c
c                 to put sun directly overhead on a point's meridian)         c
c                 (Angle between current time and solar noon)                 c
c t_since_noon  : Number of seconds since noon (after noon > 0)               c
c julday        : Julian day of the year (January 1=1; December 31=365)       c
c                                                                             c
c start_time_until_noon: number of seconds between start time and solar noon  c
c                                                               Unit: s       c
c                                                                             c
c Fs0 INTERPOLATION PARAMETERS                                                c
c nparam        : Number of Fs0 values in the list.                           c
c xilist        : Values of cosine of solar zenith angle corresponding to     c
c                 the values in Fslist                                        c
c Fslist        : Values of Fs0 corresponding to the values in xilist.        c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Local ice diffusional arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, DIMENSION(gr%nnzp):: dqc_dt_icedfs

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Local radiation arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real, dimension(gr%nnzp) :: 
     .  Frad_LW, ! Long wave radiative flux     [W/m^2]
     .  Frad_SW, ! Short wave radiative flux    [W/m^2]
     .  radht_LW,! Long wave radiative heating  [K/s]
     .  radht_SW ! Short wave radiative heating [K/s]

        REAL, DIMENSION(gr%nnzp):: LWP
        REAL, DIMENSION(gr%nnzp):: thm_rad
        REAL, DIMENSION(gr%nnzp):: rcm_rad
        REAL, DIMENSION(gr%nnzp):: rhot_rad
        REAL, DIMENSION(gr%nnzp):: dsigm

        REAL, DIMENSION(gr%nnzp):: frad_out
        REAL, DIMENSION(gr%nnzp):: frad_lw_out
        REAL, DIMENSION(gr%nnzp):: frad_sw_out

        REAL, DIMENSION(gr%nnzp):: radhtk
        REAL, DIMENSION(gr%nnzp):: radht_lw_out
        REAL, DIMENSION(gr%nnzp):: radht_sw_out

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Working arrays for subsidence interpolation !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL zsubs(5)
        REAL tsubs(6)
        REAL wt1(5), wt2(5), wt3(5), wt4(5), wt5(5), wt6(5)
        REAL w1(gr%nnzp), w2(gr%nnzp)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! LW Radiative constants !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, PARAMETER :: F0   = 107.0  ! W/m^2
        REAL, PARAMETER :: F1   = 61.0   ! W/m^2
        REAL, PARAMETER :: kap  = 100.0  ! m^2/kg

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Working arrays for SW radiation interpolation !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER, PARAMETER :: nparam = 12
        REAL, DIMENSION(nparam) :: xilist, Fslist

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! SW Radiative constants !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, PARAMETER :: radius = 1.0e-5
        REAL, PARAMETER :: AA     = 0.1
        REAL, PARAMETER :: gc     = 0.85
        REAL, PARAMETER :: omega  = 0.992

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Additional SW radiative variables !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL xi_abs, Fs0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variables used for calculating solar zenith angle                 !!
        !! (Currently not used in Nov.11 case)                               !!
        !!                                                                   !!
        !! NOTE BY ADAM SMITH, 13 April 2006                                 !!
        !! In COAMPS, these constants are called c0, c1, c2, c3, d1, d2,     !!
        !! and d3.  Since HOC uses these variables names elsewhere, I have   !!
        !! renamed them here.                                                !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL c_0, c_1, c_2, c_3, d_1, d_2, d_3
        REAL sda_t, sda_delta, sda_h, t_since_noon, julday
        REAL start_time_until_noon

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Toggle for implementing differencing method in interpolations !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        LOGICAL :: center

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Toggles for activating/deactivating forcings                  !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        LOGICAL :: icedfs_on, subs_on, lw_on, sw_on

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variable used for working within vertical arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER:: k

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! END OF VARIABLE DECLARATION !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Toggles for activating/deactivating forcings                            cc
cc To turn off a specific forcing, set the corresponding toggle to .FALSE. cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         icedfs_on = .FALSE.
         subs_on   = .TRUE.
         lw_on     = .TRUE.
         sw_on     = .TRUE.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Toggle for centered/forward differencing (in interpolations) cc
cc To use centered differencing, set the toggle to .TRUE.       cc
cc To use forward differencing, set the toggle to .FALSE.       cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         center    = .TRUE.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Calculation of Solar Declination Angle                               c
c                                                                      c
c This section of code calculates the solar declination angle based    c
c on the simulaiton's time, date, and latitude.  These parameters are  c
c applied in the scheme, and a value of xi_abs is calculated.          c
c                                                                      c
c If xi_abs <= 0.0 at any time step, we shut off shortwave radiation   c
c for that time step (by setting sw_on = .FALSE.).  If xi_abs > 0.0,   c
c we reactivate solar radiation by setting sw_on = .TRUE.              c
c                                                                      c
c This version is based on solar_declin.m, written by Dr. Vince Larson c
c and based on Liou 2.2.10.                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Modification by Adam Smith 26 June 2006                                   !!
!! It is difficult to remember to set xi_abs = 0 when we want to shut off    !!
!! solar radiation.  If sw_on = .FALSE. above, we will automatically set     !!
!! xi_abs to 0 to avoid confusion or errors.                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( .not. sw_on ) then
        xi_abs = 0.
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End of ajsmith4's Modification !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (xi_abs == 0.) then
        sw_on = .FALSE.
      else
        sw_on = .TRUE.
      end if


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Fs0 Interpolation Parameters-- these also from Kurt Kotenberg's      c
c BUGSrad output.  Fs0 changes somewhat over the range of solar zenith c
c angles, and we obtained these values by solving                      c
c Fs0 = F_vertical / xi_abs .                                          c
c                                                                      c
c The linear_interpolation function returns Fs0.                       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subsidence Parameters !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        zsubs(1) = 0
        zsubs(2) = 360
        zsubs(3) = 1090
        zsubs(4) = 1890
        zsubs(5) = 2500

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


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c SPECIAL METHOD USED TO CALCULATE RADIATION                           c
c Grid descriptions by Adam Smith, 27 June 2006                        c
c                                                                      c
c In order to verify our HOC simulations are working properly, we      c
c have first developed a series of 3D simulations using the COAMPS-LES c
c model.  This large-eddy simulation (LES) simulation uses specific    c
c methods to calculate radiation, subsidence, and other microphysical  c
c processes.  To make the two models simluate clouds as closely as     c
c possible, we use the same radiation scheme in both models.           c
c                                                                      c
c In COAMPS-LES, we use a separate subroutine, rad_lwsw, to implement  c
c all radiation code.  This allows the subroutine to be duplicated     c
c exactly in many different models.  However, the subroutine uses the  c
c COAMPS vertical grid.  Therefore, for HOC to implement this code     c
c correctly, we must modify some of our variable profiles before       c
c calling the radiation subroutine.                                    c
c                                                                      c
c The following diagram describes the differences in model grids:      c
c                                                                      c
c       COAMPS-LES                                    HOC              c
c                                                                      c
c  k= 1   (m) ----------    <MODEL TOP>    ---------- k=kk+1      (m)  c
c                                                                      c
c  k= 1   (t) ----------                   ---------- k=kk+1      (t)  c
c                                                                      c
c  k= 2   (m) ----------                   ---------- k=kk        (m)  c
c                                                                      c
c  k= 2   (t) ----------                   ---------- k=kk        (t)  c
c                                                                      c
c            .                  .                     .                c
c            .                  .                     .                c
c            .                  .                     .                c
c                                                                      c
c  k=kk-1 (m) ----------  m = mass levels  ---------- k=3         (m)  c
c                                                                      c
c  k=kk-1 (t) ----------  t = thermo       ---------- k=3         (t)  c
c                                  levels                              c
c  k=kk   (m) ----------                   ---------- k=2         (m)  c
c                                                                      c
c  k=kk   (t) ----------  kk = number of   ---------- k=2         (t)  c
c                              vertical                                c
c  k=kk+1 (m) ----------       heights     ---------- k=1         (m)  c
c                                                                      c
c //////////////////////// MODEL SURFACE ///////////////////////////// c
c                                          ---------- k=1         (t)  c
c                                                                      c
c                                                                      c
c The major difference in the grids is that HOC uses an additional     c
c thermodynamic level below the model "surface".  This means that all  c
c HOC thermodynamic heights are shifted down one vertical level, and   c
c HOC also has one fewer momentum level than COAMPS.  Therefore, we    c
c use one additional vertical level in HOC, to make sure that the      c
c vertical domain matches in both models.                              c
c                                                                      c
c Also, the COAMPS grid indices are numbered from the top of the model c
c downward, while the HOC grid indices are numbered from the bottom    c
c up.  Therefore, since we are using a COAMPS radiation scheme, we     c
c flip moisture and temperature profiles that are passed into the      c
c rad_lwsw subroutine.  The rad scheme will produce results in using   c
c the COAMPS grid scheme, so all radiation output will be flipped      c
c back to the HOC grid before being applied to the model.              c
c                                                                      c
c Finally, since the COAMPS scheme does not have a gridpoint below     c
c model surface, we add that point to all radiative output files once  c
c they are converted back to HOC setup.  This allows all averages and  c
c calculations to be done correctly.                                   c
c                                                                      c
c                                                                      c
c Computation of radiative fluxes on staggered grid                    c
c Comments by Michael Falk, 16 February 2005.                          c
c                                                                      c
c Frad (and its components Frad_LW and Frad_SW) should be computed on  c
c w points, not on mass points, which is apparent from its formulation c
c and from its location in stats_sw instead of stats_sm.  The grid     c
c looks like this:                                                     c
c                                                                      c
c                                                                      c
c -----Frad----------------------------------    k = 1  (w level)      c
c     /    \            |-dwm                                          c
c -LWP------radht----------------------------    k = 1  (mass level)   c
c     \    /            |-dmw                                          c
c -----Frad----------------------------------    k = 2  (w level)      c
c     /    \                                                           c
c -LWP------radht----------------------------    k = 2  (mass level)   c
c     \    /                                                           c
c -----Frad----------------------------------    k = 3  (w level)      c
c     /    \                                                           c
c -LWP------radht----------------------------    k = 3  (mass level)   c
c                                                                      c
c If you consider Frad to take place on mass levels, then computing    c
c LWP is a forward difference and is only first-order accurate, while  c
c if Frad computed in between LWP levels, it is a centered difference  c
c which is second-order accurate.                                      c
c                                                                      c
c The coding implementation requires that Frad depend on LWP(k) and    c
c LWP(k-1) since the w level for a given k is at a higher altitude     c
c than the mass level.  radht, back on mass levels, depends on Frad(k) c
c and Frad(k+1).                                                       c
c                                                                      c
c Additionally, these computations assume that the distance between    c
c mass levels (dsigma) is constant, and that the w levels (spaced by   c
c dsigmw) always fall exactly halfway in between the mass levels.  If  c
c this is not the case, consider dwm to be the distance between a w    c
c level and the mass level below it, and dmw to be the distance        c
c between a mass level and the w level below it.  Then, the            c
c formulation for Frad_LW, for instance, would use a weighted average: c
c                                                                      c
c (dwm/(dwm+dmw)) * lwp(k) + (dmw/(dwm+dmw)) * lwp(k-1)                c
c which, for dwm always == dmw, reduces to                             c
c (1/2) * (lwp(k)) + (1/2) * (lwp(k-1))                                c
c which is identical to the current formulation.                       c
c ((lwp(k)+lwp(k-1))/2)                                                c
c                                                                      c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        !! We only implement this section if we choose not to use the !!
        !! BUGSRAD radiation scheme.                                  !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IF ( .not. lbugsrad ) THEN

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! This code transforms these profiles from HOC grid to COAMPS !!
        !! grid.  The COAMPS-grid profiles are then passed to rad_lwsw !!
        !! for implementation.                                         !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do k = 1, gr%nnzp
            thm_rad(k) = thlm(gr%nnzp-k+1) 
     .                 + Lv/(Cp*exner(k))*rcm(gr%nnzp-k+1)
            rcm_rad(k) = rcm(gr%nnzp-k+1)
            rhot_rad(k) = rhot(gr%nnzp-k+1)
            dsigm(k) = 1.0 / gr%dzt(gr%nnzp-k+1)
          end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Calling the radiation subroutine, which uses the COAMPS    !!
        !! grid method.  All input and output profiles use the COAMPS !!
        !! grid setup.                                                !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call rad_lwsw(thm_rad, rcm_rad, rhot_rad, dsigm,
     .                  Frad_out, Frad_LW_out, Frad_SW_out,
     .                  radhtk, radht_LW_out, radht_SW_out,
     .                  gr%nnzp-1, center,
     .                  xi_abs, F0, F1, kap, radius, AA, gc, Fs0, omega,
     .                  sw_on, lw_on)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! This code transforms the radiation results back into HOC !!
        !! grid setup.  These Frad and radht arrays are actually    !!
        !! applied to the HOC model.                                !!
        !!                                                          !!
        !! The radht results are initially calculated in terms of   !!
        !! standard temperature (T).  However, HOC calculates       !!
        !! temperature in terms of potential temperature (theta).   !!
        !! Therefore, we multiply all radht results by (1.0/exner)  !!
        !! to convert from T to theta.                              !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

        END IF ! ~ lbugsrad


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Compute the loss of total water due to diffusional       !!
        !! growth of ice.  This is defined on thermodynamic levels. !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (time >= time_initial + 3600.0 .and. icedfs_on) THEN

           CALL altocu_icedf( dt, thlm, rcm, p, exner, rhot,
     .                        dqc_dt_icedfs )

        ELSE

           DO k = 1, gr%nnzp, 1
              dqc_dt_icedfs(k) = 0.0
           ENDDO

        END IF


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Using linear interpolation to calculate subsidence !!
        !! Original code by Michael Falk                      !!
        !! Added for Jun.25 case by Adam Smith, 13 April 2006 !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (time.lt.tsubs(1)) then
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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Enter the final theta-l and rtm tendencies !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

        RETURN
        END SUBROUTINE jun25_altocu_tndcy

        end module jun25
