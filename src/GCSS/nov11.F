#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: nov11.F,v 1.1 2006-12-16 00:32:45 dschanen Exp $
        module nov11

!       Description:
!       Contains subroutines for the Nov. 11 case.
!----------------------------------------------------------------------

        public :: nov11_altocu_tndcy

        private

        contains

!-----------------------------------------------------------------------
        subroutine nov11_altocu_tndcy
     .             ( time, time_initial, dt, thlm,
     .               rcm, p, exner, rhot, rtm, wmt,
     .               wmm, thl_tndcy, rt_tndcy,
     .               Frad, radht, sclr_tndcy )

        use grid_class
        use constants
        use model_flags
        use icedfs
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c INPUT VARIABLES                                                             c
c time         :  time of simulation since start                 Unit: s      c
c time_initial :  initial time of simulation                     Unit: s      c
c dt           :  length of model time step                      Unit: s      c
c thlm         :  profile of liquid water potential temperature  Unit: K      c
c rcm          :  profile of liquid (cloud) water mixing ratio   Unit: kg/kg  c
c p            :  basic state pressure                           Unit: hPa    c
c exner        :  basic state Exner function                     Unit: NONE   c
c rhot         :  density of reference state at thermodynamic levels          c
c                                                                Unit: kg/m^3 c
c                                                                             c
c OUTPUT VARIABLES                                                            c
c rtm          :  profile of total water mixing ratio            Unit: kg/kg  c
c wmt          :  profile of vertical ascent/descent at thermodynamic levels  c
c                                                                Unit: m/s    c
c wmm          :  profile of vertical ascent/descent at momentum levels       c
c                                                                Unit: m/s    c
c                                                                             c
c thl_tndcy    :  change in liq. water potential temperature due to           c
c                 radiative heating and ice diffusion            Unit: K      c
c rt_tndcy     :  change in total water due to ice diffusion     Unit: kg/kg  c
c                                                                             c
c Frad         :  profile of total radiative flux (lw + sw)      Unit: W/m^2  c
c                                                                             c
c radht        :  profile of total radiative heating (lw + sw)   Unit: K/s    c
c                                                                             c
c sclr_tndcy   :                                                 Unit: ??     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !!!!!!!!!!!!!!!!!!!!!
        !! Input variables !!
        !!!!!!!!!!!!!!!!!!!!!
        REAL, INTENT(IN):: time
        REAL, INTENT(IN):: time_initial
        REAL, INTENT(IN):: dt
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: thlm
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rcm
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: p
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: exner
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rhot

        !!!!!!!!!!!!!!!!!!!!!!
        !! Output variables !!
        !!!!!!!!!!!!!!!!!!!!!!

        REAL, DIMENSION(1:gr%nnzp), INTENT(INOUT):: rtm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: wmt
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: wmm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: thl_tndcy
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: rt_tndcy

        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Frad

        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: radht

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Output variables (optional) !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, OPTIONAL, DIMENSION(gr%nnzp, sclr_dim), 
     .                  INTENT(OUT) :: sclr_tndcy

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
c ADDITIONAL PARAMETERS FOR NOV.11 SUBSIDENCE                                 c
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
        REAL, DIMENSION(1:gr%nnzp):: dqc_dt_icedfs

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Local radiation arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, DIMENSION(1:gr%nnzp) :: 
     .  Frad_LW, ! Long wave radiative flux     [W/m^2]
     .  Frad_SW, ! Short wave radiative flux    [W/m^2]
     .  radht_LW,! Long wave radiative heating  [K/s]
     .  radht_SW ! Short wave radiative heating [K/s]
        REAL, DIMENSION(1:gr%nnzp):: LWP
        REAL, DIMENSION(1:gr%nnzp):: thm_rad
        REAL, DIMENSION(1:gr%nnzp):: rcm_rad
        REAL, DIMENSION(1:gr%nnzp):: rhot_rad
        REAL, DIMENSION(1:gr%nnzp):: dsigm

        REAL, DIMENSION(1:gr%nnzp):: frad_out
        REAL, DIMENSION(1:gr%nnzp):: frad_lw_out
        REAL, DIMENSION(1:gr%nnzp):: frad_sw_out

        REAL, DIMENSION(1:gr%nnzp):: radhtk
        REAL, DIMENSION(1:gr%nnzp):: radht_lw_out
        REAL, DIMENSION(1:gr%nnzp):: radht_sw_out

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Working arrays for subsidence interpolation !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real zsubs(7)
        real wt1(7)       ! ONLY wt1 IS NEEDED FOR NOV.11 CASE

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Subsidence constant and variables (for Nov.11 case only) !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL:: wmax, zi, dazi, dbzi, dbc, dac

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! LW Radiative constants !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, PARAMETER:: F0   = 104.0  ! W/m^2
        REAL, PARAMETER:: F1   = 62.0   ! W/m^2
        REAL, PARAMETER:: kap  = 94.2 ! m^2/kg

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Working arrays for SW radiation interpolation !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER, PARAMETER :: nparam = 2
        REAL, DIMENSION(1:nparam):: xilist, Fslist

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! SW Radiative constants !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, PARAMETER:: radius = 1.0e-5
        REAL, PARAMETER:: A      = 0.1
        REAL, PARAMETER:: gc     = 0.86
        REAL, PARAMETER:: omega  = .9965

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Additional SW radiative variables !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL:: xi_abs, Fs0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variables used for calculating solar zenith angle                 !!
        !! (Currently not used in Nov.11 case)                               !!
        !!                                                                   !!
        !! NOTE BY ADAM SMITH, 13 April 2006                                 !!
        !! In COAMPS, these constants are called c0, c1, c2, c3, d1, d2,     !!
        !! and d3.  Since HOC uses these variables elsewhere, I have renamed !!
        !! them here.                                                        !!
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
         icedfs_on = .TRUE.
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
c                                                                      c
c FOR NOV.11 CASE                                                      c
c ---------------                                                      c
c The Nov.11 case uses a constant value for the cosine of the solar    c
c zenith angle.  Because the simulation is only 4 hours long, and it   c
c starts at 18Z (11am local time), we have rather constant sunlight    c
c through the entire simulation period.  As a result, we will use a    c
c constant value of xi_abs = 0.4329.  For now, the following           c
c calculation of xi_abs will be commented out, then followed by the    c
c manual declaration of xi_abs as a constant.  We may later decide to  c
c use the code later to calculate a more accurate xi_abs.              c
c                                                                      c
c Comment by Adam Smith (ajsmith4) on 26 June 2006                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      julday = 315
c      start_time_until_noon = 7916
c      t_since_noon   = time - start_time_until_noon
c      sda_t = 2*pi*(julday-1)/365
c
c      c0 = 0.006918
c      c1 = -0.399912
c      c2 = -0.006758
c      c3 = -0.002697
c      d1 = 0.070257
c      d2 = 0.000907
c      d3 = 0.000148
c
c      sda_delta = c0 + c1*cos(sda_t) + d1*sin(sda_t) +
c     .            c2*cos(2*sda_t) + d2*sin(2*sda_t) +
c     .            c3*cos(3*sda_t) + d3*sin(3*sda_t)
c
c      sda_h = 2*pi*t_since_noon/86400
c
c      xi_abs = sin(rlat*pi/180) * sin(sda_delta) +
c     .         cos(rlat*pi/180) * cos(sda_delta) * cos(sda_h)
c
c      xi_abs = max(xi_abs,0.)

       xi_abs = .4329

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Modification by Adam Smith 26 June 2006                                   !!
!! It is difficult to remember to set xi_abs = 0 when we want to shut off    !!
!! solar radiation.  If sw_on = .FALSE. above, we will automatically set     !!
!! xi_abs to 0 to avoid confusion or errors.                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (.not. sw_on) then
        xi_abs = 0.
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End of ajsmith4's Modification !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (xi_abs .eq. 0.) then
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
c                                                                      c
c FOR NOV.11 CASE                                                      c
c ---------------                                                      c
c As explained above, the solar declination angle is presumed to be    c
c constant for Nov.11 cases.  Because of this, we also use a constant  c
c value of Fs0 for solar radiation.                                    c
c                                                                      c
c Because the Jun.25 case uses a linear_interpolation to calculate     c
c Fs0, we will duplicate that scheme here to keep all of our group's   c
c code consistent.  This means that Fslist will be a 2D array with the c
c same Fs0 at all heights.  Therefore, Fs0 will be constant over the   c
c entire model run.                                                    c
c                                                                      c
c Comment by Adam Smith on 26 June 2006                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      xilist(1) = 0.
      xilist(2) = 1.

      Fslist(1) = 1212.75
      Fslist(2) = 1212.75

      call linear_interpolation( nparam, xilist, Fslist, xi_abs, Fs0 )


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Subsidence Parameters                                                c
c                                                                      c
c FOR NOV.11 CASE                                                      c
c ---------------                                                      c
c The Nov.11 case uses a constant subsidence profile, initiated after  c
c 1 hour of model runtime.  This initial hour is used to "spinup" the  c
c simulation and produce a realistic cloud.                            c
c                                                                      c
c In jun25.F (in mjfalk's /coamps/mod/consolidated5 directory on       c
c condella), subsidence varies over time.  As a result, he uses a      c
c number of arrays defining subsidence profiles for different times.   c
c Since Nov.11's subsidence does not vary with time, we only need one  c
c of these arrays, and the rest have been removed.  The array listing  c
c different times for subsidence has been removed as well.             c
c                                                                      c
c Comment by Adam Smith on 26 June 2006                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subsidence Parameters !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        wmax =  -0.03
          zi = 1400.0
        dazi =  400.0
        dbzi = 1000.0
         dbc =  300.0
         dac =  200.0

        zsubs(1) = 0.
        zsubs(2) = zi-dbzi-dbc
        zsubs(3) = zi-dbzi
        zsubs(4) = zi
        zsubs(5) = zi+dazi
        zsubs(6) = zi+dazi+dac
        zsubs(7) = 2400.

          wt1(1) = 0.
          wt1(2) = 0.
          wt1(3) = wmax
          wt1(4) = wmax
          wt1(5) = wmax
          wt1(6) = 0.
          wt1(7) = 0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!! SPECIAL NOV.11 CONDITION FOR TOTAL WATER ABOVE CLOUD          !!
!! One hour after the initial time, the total water above cloud  !!
!! is adjusted to be 0.89 of what it previously was.             !!
!!                                                               !!
!! The conditional statement here is set so that if the timestep !!
!! is such that there is no timestep at exactly 3600.0 seconds,  !!
!! then the operation still happnens at the first timestep and   !!
!! only the first timestep after 3600.0 seconds.                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF ( time >= time_initial + 3600.0    .AND.
     .       time <  time_initial + 3600.0 + dt ) THEN

           DO k = 1, gr%nnzp, 1
              IF ( gr%zt(k) > 1555.0 ) THEN
                 rtm(k) = 0.89 * rtm(k)
              END IF
           ENDDO

        END IF


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
c     \    /            |                                        c
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
     .                  xi_abs, F0, F1, kap, radius, A, gc, Fs0, omega,
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
        !!                                                          !!
        !! NOTE FOR JUN.25 CASE                                     !!
        !! --------------------                                     !!
        !! Ice diffusion is generally not used in the Jun.25 case.  !!
        !! The code is still available here if needed.              !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (time >= time_initial + 3600.0 .and. icedfs_on) THEN

           CALL altocu_icedf( dt, thlm, rcm, p, exner, rhot,
     .                              dqc_dt_icedfs )

        ELSE

           DO k = 1, gr%nnzp, 1
              dqc_dt_icedfs(k) = 0.0
           ENDDO

        END IF


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
c Using linear interpolation scheme to interpolate subsidence        c
c                                                                    c
c FOR NOV.11 CASE                                                    c
c ---------------                                                    c
c As mentioned above, we want to implement a constant subsidence     c
c profile throughout the entire simulation (except for the 1 hour    c
c initial "spinup" period).  Because we do not have variations in    c
c subsidence over time, all tsubs sections (used in Jun.25 case)     c
c have been removed below.  Only one loop remains below, which will  c
c implement the interpolation subroutine, then return the w_ls       c
c value.  Unlike jun25, we use this w_ls value directly instead of   c
c using interpolation to calculate a value between time steps.       c
c                                                                    c
c DIAGRAM OF NOV. 11 SUBSIDENCE PROFILE                              c
c -------------------------------------                              c
c       |      <- No subsidence in this region                       c
c       |                                                            c
c       ------------ Height = (zi + dazi + dac)                      c
c        \                                                           c
c         \    <- Subsidence tapers linearly in this region          c
c          \                                                         c
c           -------- Height = (zi + dazi)                            c
c           |                                                        c
c           |  <- Subsidence equals wmax in this region              c
c           |                                                        c
c           -------- Height = zi                                     c
c           |                                                        c
c           -------- Height = (zi - dbzi)                            c
c          /                                                         c
c         /    <- Subsidence tapers linearly in this region          c
c        /                                                           c
c       ------------ Height = (zi - dbzi - dbc)                      c
c       |                                                            c
c       |      <- No subsidence in this region                       c
c                                                                    c
c                                                                    c
c Comment by Adam Smith on 26 June 2006                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do k=2,gr%nnzp
          if ( (time >= time_initial + 3600.0) .and. subs_on ) then
            call linear_interpolation(7,zsubs,wt1,gr%zt(k),wmt(k))
          else
c           If time is not yet one hour, we have no subsidence
            wmt(k) = 0.0
          end if

          wmt(1) = wmt(2)
        end do

        wmm = zt2zm(wmt)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Enter the final theta-l and rtm tendencies !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO k = 1, gr%nnzp, 1
           IF ( .not. lbugsrad ) THEN
             thl_tndcy(k) = radht(k)
     .                    - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           ELSE
             thl_tndcy(k) = - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           END IF
           rt_tndcy(k) = dqc_dt_icedfs(k)
        ENDDO

        IF ( present( sclr_tndcy ) ) THEN
          sclr_tndcy(:,SCLR_THETA) = thl_tndcy
          sclr_tndcy(:,SCLR_RT)    = rt_tndcy
        END IF

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
        END SUBROUTINE nov11_altocu_tndcy

        end module nov11
