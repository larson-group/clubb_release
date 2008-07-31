!----------------------------------------------------------------------
! $Id: clex9_nov02.F90,v 1.8 2008-07-31 16:10:43 faschinj Exp $
  module clex9_nov02

!       Description:
!       Contains subroutines for the Nov. 11 case.

!       References:
!----------------------------------------------------------------------

  implicit none

  public :: clex9_nov02_tndcy

  ! Note: bottom point not at surface, so there is no sfc
  ! subroutine

  ! Used to start the microphysics after predetermined amount of time
  logical, private ::  & 
  tdelay_lcoamps_micro, tdelay_licedfs 

!$omp   threadprivate(tdelay_lcoamps_micro, tdelay_licedfs)

  private ! Default Scope

  contains

!-----------------------------------------------------------------------
  subroutine clex9_nov02_tndcy & 
             ( time, time_initial, rlat, rlon, & 
               rcm, exner, rho, wmt, & 
               wmm, thlm_forcing, rtm_forcing, & 
               Frad, radht, Ncnm, sclrm_forcing )

!       Description:
!       Compute subsidence, radiation, and large-scale tendencies.

!       References:
!----------------------------------------------------------------------

  use grid_class, only: gr ! Variable(s)

  use grid_class, only: zt2zm ! Procedure(s)

  use constants, only: Lv, Cp ! Variable(s)

  use parameters, only: sclr_dim ! Variable(s)

  use model_flags, only: l_bugsrad, l_coamps_micro, l_licedfs ! Variable(s)

  use stats_precision, only: time_precision ! Variable(s)

  use cos_solar_zen_mod, only: cos_solar_zen ! Procedure(s)

  use interpolation, only: linear_interpolation ! Procedure(s)

  use rad_lwsw_mod, only: rad_lwsw ! Procedure(s)

  use array_index, only: iisclr_thl, iisclr_rt

 
  use stats_type, only: stat_update_var ! Procedure(s)

  use stats_variables, only:  & 
      iradht_LW, iradht_SW, zt, zm, l_stats_samp,  & ! Variable(s)
                 iFrad_SW, iFrad_LW
 

  implicit none

  ! Local constants

  integer, parameter ::  & 
  nparam = 12 ! Number of Fs0 values in the SW radiation lists

  ! LW Radiative constants
  real, parameter ::  & 
  F0   = 104.0,  & ! Coefficient for cloud top heating (see Stevens) [W/m^2]
  F1   = 62.0,   & ! Coefficient for cloud base heating (see Stevens)[W/m^2]
  kap  = 94.2      ! A "constant" according to Duynkerke eqn. 5, 
                   ! where his value is 130 m^2/kg [m^2/kg]

  ! SW Radiative constants
  real, parameter ::  & 
  radius = 1.0e-5, & ! Effective droplet radius                      [m]
  A      = 0.1,    & ! Albedo -- sea surface, according to Lenderink [-]
  gc     = 0.86,   & ! Asymmetry parameter, "g" in Duynkerke.        [-]
  omega  = 0.9965    ! Single-scattering albedo                      [-]


  ! Toggles for activating/deactivating forcings
  logical, parameter ::  & 
  subs_on   = .true., & 
  lw_on     = .true.

  ! Toggle for centered/forward differencing (in interpolations)
  ! To use centered differencing, set the toggle to .true.
  ! To use forward differencing, set the toggle to  .false.
  logical, parameter :: & 
  center = .false.

  ! Input variables
  real(kind=time_precision), intent(in) :: & 
  time,            & ! Current time          [s]
  time_initial       ! Initial time          [s]

  real, intent(in) :: & 
  rlat,            & ! Latitude              [degrees_N]
  rlon               ! Longitude             [degrees_E]

  real, intent(in), dimension(gr%nnzp) :: & 
  rcm,     & ! Cloud water mixing ratio      [kg/kg]
  exner,   & ! Exner function                [-]
  rho       ! Density                       [kg/m^3]

  ! Output variables
  real, intent(out), dimension(gr%nnzp) :: & 
  wmt,             & ! Mean vertical wind on the thermo. grid  [m/s]
  wmm,             & ! Mean vertical wind on the moment. grid  [m/s]
  thlm_forcing,    & ! Theta_l forcing                         [K/s]
  rtm_forcing,     & ! Total water forcing                     [kg/kg/s]
  Frad,            & ! Radiative flux                          [W/m^2]
  radht,           & ! Radiative heating                       [K/s]
  Ncnm               ! Cloud nuclei number concentration       [num/m^3]

  ! Output variables (optional)
  real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
  sclrm_forcing   ! Passive scalar forcing                  [units/s]


  ! Local variables

  ! Local radiation arrays
  real, dimension(gr%nnzp) ::  & 
  Frad_LW,  & ! Long wave radiative flux     [W/m^2]
  Frad_SW,  & ! Short wave radiative flux    [W/m^2]
  radht_LW, & ! Long wave radiative heating  [K/s]
  radht_SW    ! Short wave radiative heating [K/s]

  real, dimension(gr%nnzp) ::  & 
!     .  LWP,       ! Liquid water path                              [kg/m^2]
  rcm_rad,    & ! Flipped array of liq. water mixing ratio       [kg/kg]
  rhot_rad,   & ! Flipped array of air density                   [kg/m^3]
  dsigm,      & ! Flipped array of grid spacing                  [m]
  coamps_zm,  & ! Flipped array of momentum level altitudes      [m]
  coamps_zt     ! Flipped array of thermodynamic level altitudes [m]

  real, dimension(gr%nnzp) ::  & 
  frad_out,    & ! Flipped array of radiaive flux            [W/m^2]
  frad_lw_out, & ! Flipped array of LW radiative flux        [W/m^2]
  frad_sw_out    ! Flipped array of SW radiative flux        [W/m^2] 

  real, dimension(gr%nnzp) ::  & 
  radhtk,       & ! Flipped array of radiative heating       [K/s]
  radht_lw_out, & ! Flipped array of LW radiative heating    [K/s]
  radht_sw_out    ! Flipped array of SW radiative heating    [K/s]

  ! Working arrays for subsidence interpolation
  real, dimension(7) ::  & 
  zsubs, & ! Heights at which wmt data is supplied (used for subsidence interpolation) [m]
  wt1      ! ONLY wt1 IS NEEDED FOR NOV.11 CASE

  ! Subsidence constant and variables
  real :: & 
  wmax,  & ! Defines value of maximum subsidence in profile  [cm/s]
  zi,    & ! Defines approx. height of inversion within cloud 
           ! (subsidence is equal to wmax at this height) [m]
  dazi,  & ! Defines height above inversion (above this height 
           ! subsidence linearly tapers off to zero)     [m]
  dbzi,  & ! Defines height above inversion (below this height, 
           ! subsidence linearly tapers off to zero)    [m]
  dbc,   & ! Defines height below cloud (at / below this height, we have NO subsidence) [m]
  dac      ! Defines height above cloud (at / above this height, we have NO subsidence) [m]

  ! Working arrays for SW radiation interpolation

  real, dimension(nparam) ::  & 
  xilist, & ! Values of cosine of solar zenith angle corresponding 
            !   to the values in Fslist
  Fslist    ! Values of Fs0 corresponding to the values in xilist.


  ! Additional SW radiative variables

  real :: & 
  xi_abs, & ! Cosine of the solar zenith angle  [-]
  Fs0       ! The incident incoming SW insolation at cloud top in the
            !   direction of the incoming beam (not the vertical) [W/m^2]

  logical :: sw_on
 
  ! Variable used for working within vertical arrays

  integer :: k

  sw_on = .true. ! This is necessay to use the xi_abs value below
                 ! Joshua Fasching June 2008

!-----------------------------------------------------------------------
! FOR NOV.11 CASE
! ---------------
! The Nov.11 case uses a constant value for the cosine of the solar
! zenith angle.  Because the simulation is only 4 hours long, and it
! starts at 18Z (11am local time), we have rather constant sunlight
! through the entire simulation period.  As a result, we will use a
! constant value of xi_abs = 0.4329.  For now, the following
! calculation of xi_abs will be commented out, then followed by the
! manual declaration of xi_abs as a constant.  We may later decide to
! use the code later to calculate a more accurate xi_abs.
!
! Comment by Adam Smith (ajsmith4) on 26 June 2006
!-----------------------------------------------------------------------

!      Replaced the calculation based on time since solar noon 
!       etc., with a generalized function based on time and lat/lon.
!       -dschanen 5 Jan 2007

!      NOTE: The results from this function match the xi_abs results
!            obtained using the COAMPS xi_abs method.  Therefore, we
!            use the "cos_solar_zen" function for this case.
 xi_abs = real( cos_solar_zen( 2, 11, 2001, time, rlat, rlon ) )

!**********************************************************************
! Addition by Adam Smith, 02 April 2008
! Adding a correction to xi_abs so that it cannot be less than zero.
! "xi_abs" must be non-zero to allow linear interpolation to calculate
! the correct amount of solar radiative flux.
!**********************************************************************
 xi_abs = max( xi_abs, 0. )
!****************************
! End of ajsmith4's addition
!****************************

!-----------------------------------------------------------------------
! Modification by Adam Smith 26 June 2006
! It is difficult to remember to set xi_abs = 0 when we want to shut off
! solar radiation.  If sw_on = .FALSE. above, we will automatically set
! xi_abs to 0 to avoid confusion or errors.
!-----------------------------------------------------------------------

if (.not. sw_on) then
  xi_abs = 0.0
end if

!-----------------------------------
! End of ajsmith4's Modification 
!-----------------------------------

if (xi_abs == 0.0 ) then
  sw_on = .false.
else
  sw_on = .true.
end if



!-----------------------------------------------------------------------
! Fs0 Interpolation Parameters-- these also from Kurt Kotenberg's
! BUGSrad output.  Fs0 changes somewhat over the range of solar zenith
! angles, and we obtained these values by solving
! Fs0 = F_vertical / xi_abs . 
!
! The linear_interpolation function returns Fs0.
!
! FOR NOV.11 CASE
! ---------------
! As explained above, the solar declination angle is presumed to be
! constant for Nov.11 cases.  Because of this, we also use a constant
! value of Fs0 for solar radiation.
!
! Because the Jun.25 case uses a linear_interpolation to calculate
! Fs0, we will duplicate that scheme here to keep all of our group's
! code consistent.  This means that Fslist will be a 2D array with the
! same Fs0 at all heights.  Therefore, Fs0 will be constant over the
! entire model run.
!
! Comment by Adam Smith on 26 June 2006
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

Fslist(1) = 0.0
Fslist(2) = 715.86
Fslist(3) = 1073.577
Fslist(4) = 1165.0905
Fslist(5) = 1204.7033
Fslist(6) = 1227.6898
Fslist(7) = 1243.1772
Fslist(8) = 1254.5893
Fslist(9) = 1263.5491
Fslist(10) = 1270.8668
Fslist(11) = 1277.0474
Fslist(12) = 1282.3994

call linear_interpolation( nparam, xilist, Fslist, xi_abs, Fs0 )

!-----------------------------------------------------------------------
!
! Subsidence Parameters
!
! FOR NOV.11 CASE
! ---------------
! The Nov.11 case uses a constant subsidence profile, initiated after
! 1 hour of model runtime.  This initial hour is used to "spinup" the
! simulation and produce a realistic cloud.
!
! In jun25.F (in mjfalk's /coamps/mod/consolidated5 directory on
! condella), subsidence varies over time.  As a result, he uses a
! number of arrays defining subsidence profiles for different times.
! Since Nov.11's subsidence does not vary with time, we only need one
! of these arrays, and the rest have been removed.  The array listing
! different times for subsidence has been removed as well.
! 
! Comment by Adam Smith on 26 June 2006
!
!-----------------------------------------------------------------------

  !----------------------
  ! Subsidence Parameters
  !----------------------
  wmax =  0.007
    zi = 1600.0 + gr%zm(1)
  dazi = 1000.0
  dbzi = 1000.0
   dbc =  300.0
   dac =  200.0

  zsubs(1) = 0. + gr%zm(1)
  zsubs(2) = zi-dbzi-dbc
  zsubs(3) = zi-dbzi
  zsubs(4) = zi
  zsubs(5) = zi+dazi
  zsubs(6) = zi+dazi+dac
  zsubs(7) = 3000. + gr%zm(1)

    wt1(1) = 0.
    wt1(2) = 0.
    wt1(3) = wmax
    wt1(4) = wmax
    wt1(5) = wmax
    wt1(6) = 0.
    wt1(7) = 0.

!-----------------------------------------------------------------------       
! SPECIAL NOV.11 CONDITION FOR TOTAL WATER ABOVE CLOUD
! One hour after the initial time, the total water above cloud
! is adjusted to be 0.89 of what it previously was.
!
! The conditional statement here is set so that if the timestep
! is such that there is no timestep at exactly 3600.0 seconds,
! then the operation still happnens at the first timestep and
! only the first timestep after 3600.0 seconds.
!-----------------------------------------------------------------------       
!        if ( time >= time_initial + 3600.0  .and.
!     .       time <  time_initial + 3600.0 + dt ) then
!
!          do k = 1, gr%nnzp, 1
!            if ( gr%zt(k) > ( 2587.5 + gr%zm(1) ) ) then
!               rtm(k) = 0.89 * rtm(k)
!            end if
!          end do
!
!        end if 


!-----------------------------------------------------------------------
! SPECIAL METHOD USED TO CALCULATE RADIATION
! Grid descriptions by Adam Smith, 27 June 2006
!
! In order to verify our CLUBB simulations are working properly, we
! have first developed a series of 3D simulations using the COAMPS-LES
! model.  This large-eddy simulation (LES) simulation uses specific
! methods to calculate radiation, subsidence, and other microphysical
! processes.  To make the two models simluate clouds as closely as
! possible, we use the same radiation scheme in both models.
! 
! In COAMPS-LES, we use a separate subroutine, rad_lwsw, to implement
! all radiation code.  This allows the subroutine to be duplicated
! exactly in many different models.  However, the subroutine uses the
! COAMPS vertical grid.  Therefore, for CLUBB to implement this code
! correctly, we must modify some of our variable profiles before
! calling the radiation subroutine.
!
! The following diagram describes the differences in model grids
! (see "ADDITIONAL NOTE" for important grid information):
!
!       COAMPS-LES                                   CLUBB
!
!  k= 1   (m) ----------    <MODEL TOP>    ---------- k=kk+1      (m)
!
!  k= 1   (t) ----------                   ---------- k=kk+1      (t)
!
!  k= 2   (m) ----------                   ---------- k=kk        (m)
!
!  k= 2   (t) ----------                   ---------- k=kk        (t)
!
!            .                  .                     .
!            .                  .                     .
!            .                  .                     .
!
!  k=kk-1 (m) ----------  m = momentum     ---------- k=3         (m)
!                                  levels
!  k=kk-1 (t) ----------  t = thermo       ---------- k=3         (t)
!                                  levels
!  k=kk   (m) ----------                   ---------- k=2         (m)
!
!  k=kk   (t) ----------  kk = number of   ---------- k=2         (t)
!                              vertical
!  k=kk+1 (m) ----------       heights     ---------- k=1         (m)
!
! //////////////////////// MODEL SURFACE /////////////////////////////
!                                          ---------- k=1         (t)
!
!
! The major difference in the grids is that CLUBB uses an additional
! thermodynamic level below the model "surface".  This means that all
! CLUBB thermodynamic heights are shifted down one vertical level, and
! CLUBB also has one fewer momentum level than COAMPS.  Therefore, we
! use one additional vertical level in CLUBB, to make sure that the
! vertical domain matches in both models.
!
! ADDITIONAL NOTE: In order to reconcile the COAMPS and CLUBB grids,
!                  the uppermost level of the COAMPS grid is ignored,
!                  eliminating the need to add level kk+1 to the CLUBB
!                  grid.  For the purposes of this code, the COAMPS
!                  grid levels are re-indexed to start from level 1
!                  at the uppermost useful COAMPS grid level (which
!                  was previously referred to as level 2 in the above
!                  diagram and is at the same altitude as CLUBB
!                  level kk).  Likewise, the lowermost COAMPS grid
!                  level is now indexed at kk (rather than kk+1 in the
!                  above diagram).  Brian Griffin; May 10, 2008.
!
! Also, the COAMPS grid indices are numbered from the top of the model
! downward, while the CLUBB grid indices are numbered from the bottom
! up.  Therefore, since we are using a COAMPS radiation scheme, we
! flip moisture and temperature profiles that are passed into the
! rad_lwsw subroutine.  The rad scheme will produce results in using
! the COAMPS grid scheme, so all radiation output will be flipped
! back to the CLUBB grid before being applied to the model.
!
! Finally, since the COAMPS scheme does not have a gridpoint below
! model surface, we add that point to all radiative output files once
! they are converted back to CLUBB setup.  This allows all averages and
! calculations to be done correctly.
!
!
! Computation of radiative fluxes on staggered grid
! Comments by Michael Falk, 16 February 2005.
!
! Frad (and its components Frad_LW and Frad_SW) should be computed on
! w points, not on mass points, which is apparent from its formulation
! and from its location in stats_sw instead of stats_sm.  The grid
! looks like this:
! 
!
! -----Frad----------------------------------    k = 1  (w level)
!     /    \            |-dwm
! -LWP------radht----------------------------    k = 1  (mass level)
!     \    /            |
! -----Frad----------------------------------    k = 2  (w level)
!     /    \
! -LWP------radht----------------------------    k = 2  (mass level)
!     \    /
! -----Frad----------------------------------    k = 3  (w level)
!     /    \
! -LWP------radht----------------------------    k = 3  (mass level)
!
! If you consider Frad to take place on mass levels, then computing
! LWP is a forward difference and is only first-order accurate, while
! if Frad computed in between LWP levels, it is a centered difference
! which is second-order accurate.
!
! The coding implementation requires that Frad depend on LWP(k) and
! LWP(k-1) since the w level for a given k is at a higher altitude
! than the mass level.  radht, back on mass levels, depends on Frad(k)
! and Frad(k+1).
!
! ADDITIONAL NOTE: For clarification of terminology, a w level on the
!                  COAMPS grid is equivalent to a momentum (m) level
!                  on the CLUBB grid, and a mass level on the COAMPS
!                  grid is equivalent to a thermodynamic (t) level on
!                  the CLUBB grid.  Brian Griffin; May 10, 2008.
!
! Additionally, these computations assume that the distance between
! mass levels (dsigma) is constant, and that the w levels (spaced by
! dsigmw) always fall exactly halfway in between the mass levels.  If
! this is not the case, consider dwm to be the distance between a w
! level and the mass level below it, and dmw to be the distance
! between a mass level and the w level below it.  Then, the
! formulation for Frad_LW, for instance, would use a weighted average:
!
! (dwm/(dwm+dmw)) * lwp(k) + (dmw/(dwm+dmw)) * lwp(k-1)
! which, for dwm always == dmw, reduces to
! (1/2) * (lwp(k)) + (1/2) * (lwp(k-1))
! which is identical to the current formulation.
! ((lwp(k)+lwp(k-1))/2)
!
! ADDITIONAL NOTE: The CLUBB parameterization is now set up to be
!                  compatible with the use of a stretched
!                  (or unevenly-spaced) grid, as well as with the use
!                  of an evenly-spaced grid.  Interpolation functions
!                  are used to compute any weighted averages, rather
!                  than using general numbers such as (1/2), which is
!                  compatible only with an evenly-spaced grid.
!                  Brian Griffin; May 10, 2008.
!
!
!-----------------------------------------------------------------------

  !---------------------------------------------------------------
  ! We only implement this section if we choose not to use the
  ! BUGSRAD radiation scheme
  !---------------------------------------------------------------
  if ( .not. l_bugsrad ) then

  !---------------------------------------------------------------
  ! This code transforms these profiles from CLUBB grid to COAMPS
  ! grid.  The COAMPS-grid profiles are then passed to rad_lwsw
  ! for implementation.
  !---------------------------------------------------------------
    do k = 1, gr%nnzp
      rcm_rad(k) = rcm(gr%nnzp-k+1)
      rhot_rad(k) = rho(gr%nnzp-k+1)
      dsigm(k) = 1.0 / gr%dzt(gr%nnzp-k+1)
      coamps_zm(k) = gr%zm(gr%nnzp-k+1)
      coamps_zt(k) = gr%zt(gr%nnzp-k+1)
    enddo

  !---------------------------------------------------------------
  ! Calling the radiation subroutine, which uses the COAMPS
  ! grid method.  All input and output profiles use the COAMPS
  ! grid setup.
  !---------------------------------------------------------------
    call rad_lwsw( rcm_rad, rhot_rad, dsigm, & 
                   coamps_zm, coamps_zt, & 
                   Frad_out, Frad_LW_out, Frad_SW_out, & 
                   radhtk, radht_LW_out, radht_SW_out, & 
                   gr%nnzp-1, center, & 
                   xi_abs, F0, F1, kap, radius, A, gc, Fs0, omega, & 
                   sw_on, lw_on )

  !---------------------------------------------------------------
  ! This code transforms the radiation results back into CLUBB
  ! grid setup.  These Frad and radht arrays are actually
  ! applied to the CLUBB model.
  !
  ! The radht results are initially calculated in terms of
  ! standard temperature (T).  However, CLUBB calculates
  ! temperature in terms of potential temperature (theta).
  ! Therefore, we multiply all radht results by (1.0/exner)
  ! to convert from T to theta.
  !---------------------------------------------------------------
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

  end if ! ~ l_bugsrad

  if ( time == time_initial ) then

    if ( l_coamps_micro ) then
      ! Turn off microphysics for now, re-enable at
      ! time = 3600.
      l_coamps_micro        = .false.
      tdelay_lcoamps_micro = .true.

    else if ( l_licedfs ) then
      l_licedfs        = .false.
      tdelay_licedfs = .true.

    else
      tdelay_lcoamps_micro = .false.
      tdelay_licedfs       = .false.
    end if

  end if

  if ( time >= ( time_initial + 3600.0 )  & 
            .and. tdelay_licedfs ) then
  !---------------------------------------------------------------
  ! Compute the loss of total water due to diffusional
  ! growth of ice.  This is defined on thermodynamic levels.
  !---------------------------------------------------------------
    l_licedfs = .true.
    

  else if ( time == ( time_initial + 3600.0 )  & 
            .and. tdelay_lcoamps_micro ) then

  !---------------------------------------------------------------
  ! Start COAMPS micro after predefined time delay
  !---------------------------------------------------------------

    l_coamps_micro = .true.

    Ncnm(1:gr%nnzp) = 0.0

  end if


!---------------------------------------------------------------------
!
! Using linear interpolation scheme to interpolate subsidence
!
! FOR NOV.11 CASE
! ---------------
! As mentioned above, we want to implement a constant subsidence
! profile throughout the entire simulation (except for the 1 hour 
! initial "spinup" period).  Because we do not have variations in
! subsidence over time, all tsubs sections (used in Jun.25 case)
! have been removed below.  Only one loop remains below, which will
! implement the interpolation subroutine, then return the w_ls
! value.  Unlike jun25, we use this w_ls value directly instead of
! using interpolation to calculate a value between time steps.
!
! DIAGRAM OF NOV. 11 SUBSIDENCE PROFILE
! -------------------------------------
!       |      <- No subsidence in this region
!       |
!       ------------ Height = (zi + dazi + dac)
!        \
!         \    <- Subsidence tapers linearly in this region
!          \
!           -------- Height = (zi + dazi)
!           |
!           |  <- Subsidence equals wmax in this region
!           |
!           -------- Height = zi
!           |
!           -------- Height = (zi - dbzi)
!          /
!         /    <- Subsidence tapers linearly in this region
!        /
!       ------------ Height = (zi - dbzi - dbc)
!       |
!       |      <- No subsidence in this region
!
!
! Comment by Adam Smith on 26 June 2006
!-----------------------------------------------------------------------

  do k=2,gr%nnzp
    if ( (time >= time_initial + 3600.0) .and. subs_on ) then
      call linear_interpolation( 7, zsubs, wt1, gr%zt(k), wmt(k) )
    else
!           If time is not yet one hour, we have no subsidence
      wmt(k) = 0.0
    end if

    wmt(1) = wmt(2)
  end do

  wmm = zt2zm(wmt)

  ! Enter the final theta-l and rtm tendencies
  do k = 1, gr%nnzp, 1

    if ( .not. l_bugsrad ) then
      thlm_forcing(k) = radht(k)
    else
      thlm_forcing(k) = 0.
    end if

    rtm_forcing(k) = 0.
  end do

  ! Test scalars with thetal and rt if desired
  if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
  if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

 
  ! Save LW and SW components of radiative heating and
  ! radiative flux based on simplified radiation.
  if ( .not.l_bugsrad .and. l_stats_samp ) then
    call stat_update_var( iradht_LW, radht_LW, zt )

    call stat_update_var( iradht_SW, radht_SW, zt )

    call stat_update_var( iFrad_SW, Frad_SW, zm )

    call stat_update_var( iFrad_LW, Frad_LW, zm )

  end if
 

  return
  end subroutine clex9_nov02_tndcy

  end module clex9_nov02
