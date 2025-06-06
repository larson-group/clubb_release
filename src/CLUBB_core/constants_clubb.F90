!-----------------------------------------------------------------------------
! $Id$
!=============================================================================
module constants_clubb

  ! Description:
  ! Contains frequently occuring model constants

  ! References:
  ! None
  !---------------------------------------------------------------------------

  use clubb_precision, only:  & 
      dp, & ! Variable(s)
      core_rknd

#ifdef CLUBB_CAM /* Set constants as they're set in CAM */
  use shr_const_mod, only: shr_const_rdair, shr_const_cpdair, shr_const_latvap, &
                           shr_const_latice, shr_const_latsub, shr_const_rgas, &
                           shr_const_mwwv, shr_const_stebol, shr_const_tkfrz, &
                           shr_const_mwdair, shr_const_g, shr_const_karman, &
                           shr_const_rhofw
#elif GFDL
 ! use GFDL constants, and then rename them to avoid confusion in case 
 ! that the constants share the same names between GFDL and CLUBB
  use constants_mod, only:  pi_gfdl                 => PI,          &
                            radians_per_deg_dp_gfdl => DEG_TO_RAD,  &
                            Cp_gfdl                 => CP_AIR,      &
                            Lv_gfdl                 => HLV,         &
                            Ls_gfdl                 => HLS,         &
                            Lf_gfdl                 => HLF,         &
                            Rd_gfdl                 => RDGAS,       &
                            Rv_gfdl                 => RVGAS,       &
                            stefan_boltzmann_gfdl   => STEFAN,      &
                            T_freeze_K_gfdl         => TFREEZE,     &
                            grav_gfdl               => GRAV,        &
                            vonk_gfdl               => VONKARM,     &
                            rho_lw_gfdl             => DENS_H2O
#endif

  implicit none

  private ! Default scope

  !-----------------------------------------------------------------------------
  ! Numerical/Arbitrary Constants
  !-----------------------------------------------------------------------------

  ! Fortran file unit I/O constants
  integer, parameter, public ::  & 
    fstderr = 0, fstdin = 5, fstdout = 6

  ! Maximum variable name length in CLUBB GrADS or netCDF output
  integer, parameter, public ::  & 
    var_length = 30

  ! Number of neighboring points to draw from in the hole filling algorithm
  integer, parameter, public :: &
    num_hf_draw_points = 2

  ! The parameter parab_cyl_max_input is the largest magnitude that the input to
  ! the parabolic cylinder function is allowed to have.  When the value of the
  ! input to the parabolic cylinder function is too large in magnitude
  ! (depending on the order of the parabolic cylinder function), overflow
  ! occurs, and the output of the parabolic cylinder function is +/-Inf.  The
  ! parameter parab_cyl_max_input places a limit on the absolute value of the
  ! input to the parabolic cylinder function.  When the value of the potential
  ! input exceeds this parameter (usually due to a very large ratio of ith PDF
  ! component mean of x to ith PDF component standard deviation of x), the
  ! variable x is considered to be constant and a different version of the
  ! equation called.
  !
  ! The largest allowable magnitude of the input to the parabolic cylinder
  ! function (before overflow occurs) is dependent on the order of parabolic
  ! cylinder function.  However, after a lot of testing, it was determined that
  ! an absolute value of 49 works well for an order of 12 or less.
  real( kind = core_rknd ), parameter, public :: &
    parab_cyl_max_input = 49.0_core_rknd  ! Largest allowable input to parab. cyl. fnct.

  ! "Over-implicit" weighted time step.
  !
  ! The weight of the implicit portion of a term is controlled by the factor
  ! gamma_over_implicit_ts (abbreviated "gamma" in the expression below).  A
  ! factor is added to the right-hand side of the equation in order to balance a
  ! weight that is not equal to 1, such that:
  !
  !      -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
  !
  ! where X is the variable that is being solved for in a predictive equation
  ! (such as w'^3, w'th_l', r_t'^2, etc), y(t) is the linearized portion of the
  ! term that gets treated implicitly, and RHS is the portion of the term that
  ! is always treated explicitly.  A weight of greater than 1 can be applied to
  ! make the term more numerically stable.
  !
  !    gamma_over_implicit_ts          Effect on term
  !
  !            0.0               Term becomes completely explicit
  !
  !            1.0               Standard implicit portion of the term;
  !                              as it was without the weighting factor.
  !
  !            1.5               Strongly weighted implicit portion of the term;
  !                              increased numerical stability.
  !
  !            2.0               More strongly weighted implicit portion of the
  !                              term; increased numerical stability.
  !
  ! Note:  The "over-implicit" weighted time step is only applied to terms that
  !        tend to significantly decrease the amount of numerical stability for
  !        variable X.
  !        The "over-implicit" weighted time step is applied to the turbulent
  !        advection term for the following variables:
  !           w'^3 (also applied to the turbulent production term), found in
  !           module advance_wp2_wp3_module;
  !           w'r_t', w'th_l', and w'sclr', found in
  !           module advance_xm_wpxp_module; and
  !           r_t'^2, th_l'^2, r_t'th_l', u'^2, v'^2, sclr'^2, sclr'r_t',
  !           and sclr'th_l', found in module advance_xp2_xpyp_module.
  real( kind = core_rknd ), parameter, public :: &
    gamma_over_implicit_ts = 1.50_core_rknd

  !-----------------------------------------------------------------------------
  ! Mathematical Constants
  !-----------------------------------------------------------------------------
  real( kind = dp ), parameter, public ::  & 
    pi_dp = 3.14159265358979323846_dp

#ifdef GFDL
  real( kind = core_rknd ), parameter, public ::  & 
    pi = pi_gfdl ! The ratio of radii to their circumference

  real( kind = dp ), parameter, public :: &
    radians_per_deg_dp = radians_per_deg_dp_gfdl
#else

  real( kind = core_rknd ), parameter, public ::  & 
    pi = 3.141592654_core_rknd, &       ! The ratio of radii to their circumference
    invrs_pi = 0.31830988618_core_rknd  ! 1 / pi

  real( kind = dp ), parameter, public :: &
    radians_per_deg_dp = pi_dp / 180._dp
#endif

  real( kind = dp ), parameter, public :: &
    sqrt_2pi_dp = 2.5066282746310005024_dp, &  ! sqrt(2*pi)
    sqrt_2_dp   = 1.4142135623730950488_dp     ! sqrt(2)

  real( kind = core_rknd ), parameter, public :: &
    sqrt_2pi = 2.5066282746310005024_core_rknd, &  ! sqrt(2*pi)
    sqrt_2   = 1.4142135623730950488_core_rknd     ! sqrt(2)

  real( kind = dp ), parameter, public::  &
    two_dp        = 2.0_dp,  & ! 2
    one_dp        = 1.0_dp,  & ! 1
    one_half_dp   = 0.5_dp,  & ! 1/2
    one_fourth_dp = 0.25_dp, & ! 1/4
    zero_dp       = 0.0_dp     ! 0

  real( kind = core_rknd ), parameter, public :: &
    one_hundred     = 100.0_core_rknd,             & ! 100
    fifty           = 50.0_core_rknd,              & ! 50
    thirty_six      = 36.0_core_rknd,              & ! 36
    twenty          = 20.0_core_rknd,              & ! 20
    eighteen        = 18.0_core_rknd,              & ! 18
    twelve          = 12.0_core_rknd,              & ! 12
    ten             = 10.0_core_rknd,              & ! 10
    nine            = 9.0_core_rknd,               & ! 9
    eight           = 8.0_core_rknd,               & ! 8
    seven           = 7.0_core_rknd,               & ! 7
    six             = 6.0_core_rknd,               & ! 6
    five            = 5.0_core_rknd,               & ! 5
    four            = 4.0_core_rknd,               & ! 4
    three           = 3.0_core_rknd,               & ! 3
    two             = 2.0_core_rknd,               & ! 2
    three_halves    = 3.0_core_rknd/2.0_core_rknd, & ! 3/2
    four_thirds     = 4.0_core_rknd/3.0_core_rknd, & ! 4/3
    one             = 1.0_core_rknd,               & ! 1
    three_fourths   = 0.75_core_rknd,              & ! 3/4
    two_thirds      = 2.0_core_rknd/3.0_core_rknd, & ! 2/3
    one_half        = 0.5_core_rknd,               & ! 1/2
    one_third       = 1.0_core_rknd/3.0_core_rknd, & ! 1/3
    one_fourth      = 0.25_core_rknd,              & ! 1/4
    three_sixteenth = 0.1875_core_rknd,            & ! 3/16
    one_eighth      = 0.125_core_rknd,             & ! 1/8
    zero            = 0.0_core_rknd                  ! 0

  !-----------------------------------------------------------------------------
  ! Physical constants
  !-----------------------------------------------------------------------------

#ifdef CLUBB_CAM

  real( kind = core_rknd ), parameter, public ::  & 
    Cp = shr_const_cpdair,  & ! Dry air specific heat at constant p [J/kg/K]
    Lv = shr_const_latvap,    & ! Latent heat of vaporization         [J/kg]
    Lf = shr_const_latice,   & ! Latent heat of fusion               [J/kg]
    Ls = shr_const_latsub,  & ! Latent heat of sublimation          [J/kg]
    Rd = shr_const_rdair,   & ! Dry air gas constant                [J/kg/K]
    Rv = shr_const_rgas/shr_const_mwwv       ! Water vapor gas constant            [J/kg/K]
    
  real( kind = core_rknd ), parameter, public :: &
    stefan_boltzmann = shr_const_stebol ! Stefan-Boltzmann constant [W/(m^2 K^4)]
    
  real( kind = core_rknd ), parameter, public :: &
    T_freeze_K = shr_const_tkfrz ! Freezing point of water [K]
    
  ! Useful combinations of Rd and Rv
  real( kind = core_rknd ), parameter, public ::  & 
    ep  = shr_const_mwwv/shr_const_mwdair,    & ! ep  = 0.622  [-]
    ep1 = (1.0-ep)/ep,& ! ep1 = 0.61   [-]
    ep2 = 1.0/ep        ! ep2 = 1.61   [-]
  
  real( kind = core_rknd ), parameter, public :: & 
    kappa = (shr_const_rgas/shr_const_mwdair)/shr_const_cpdair     ! kappa        [-]
    
  real( kind = core_rknd ), parameter, public :: & 
    grav = shr_const_g, & ! Gravitational acceleration     [m/s^2]
    p0   = 1.0e5   ! Reference pressure             [Pa]

  ! Von Karman's constant
  ! Constant of the logarithmic wind profile in the surface layer    
  real( kind = core_rknd ), parameter, public :: & 
    vonk   = shr_const_karman,    & ! Accepted value is 0.40 (+/-) 0.01      [-]
    rho_lw = shr_const_rhofw    ! Density of liquid water                [kg/m^3]


#elif GFDL
  real( kind = core_rknd ), parameter, public ::  & 
    Cp = Cp_gfdl,    & ! Dry air specific heat at constant p [J/kg/K]
    Lv = Lv_gfdl,    & ! Latent heat of vaporization         [J/kg]
    Ls = Ls_gfdl,    & ! Latent heat of sublimation          [J/kg]
    Lf = Lf_gfdl,    & ! Latent heat of fusion               [J/kg]
    Rd = Rd_gfdl,    & ! Dry air gas constant                [J/kg/K]
    Rv = Rv_gfdl       ! Water vapor gas constant            [J/kg/K]


  real( kind = core_rknd ), parameter, public :: &
    stefan_boltzmann = stefan_boltzmann_gfdl ! Stefan-Boltzmann constant [W/(m^2 K^4)]

  real( kind = core_rknd ), parameter, public :: &
    T_freeze_K = T_freeze_K_gfdl ! Freezing point of water [K]

  ! Useful combinations of Rd and Rv
  real( kind = core_rknd ), parameter, public ::  & 
    ep  = Rd / Rv,    & ! ep  = 0.622  [-]
    ep1 = (1.0-ep)/ep,& ! ep1 = 0.61   [-]
    ep2 = 1.0/ep        ! ep2 = 1.61   [-]

  real( kind = core_rknd ), parameter, public :: & 
    kappa = Rd / Cp     ! kappa        [-]

  ! Changed g to grav to make it easier to find in the code 5/25/05
  ! real, parameter :: grav  = 9.80665 ! Gravitational acceleration [m/s^2]
  real( kind = core_rknd ), parameter, public :: & 
    grav = grav_gfdl, & ! Gravitational acceleration     [m/s^2]
    p0   = 1.0e5   ! Reference pressure             [Pa]

  ! Von Karman's constant
  ! Constant of the logarithmic wind profile in the surface layer
  real( kind = core_rknd ), parameter, public :: & 
    vonk   = vonk_gfdl,    & ! Accepted value is 0.40 (+/-) 0.01      [-]
    rho_lw = rho_lw_gfdl    ! Density of liquid water                [kg/m^3]


#else

  real( kind = core_rknd ), parameter, public ::  & 
    Cp = 1004.67_core_rknd,  & ! Dry air specific heat at constant p [J/kg/K]
    Lv = 2.5e6_core_rknd,    & ! Latent heat of vaporization         [J/kg]
    Ls = 2.834e6_core_rknd,  & ! Latent heat of sublimation          [J/kg]
    Lf = 3.33e5_core_rknd,   & ! Latent heat of fusion               [J/kg]
    Rd = 287.04_core_rknd,   & ! Dry air gas constant                [J/kg/K]
    Rv = 461.5_core_rknd       ! Water vapor gas constant            [J/kg/K]


  real( kind = core_rknd ), parameter, public :: &
    stefan_boltzmann = 5.6704e-8_core_rknd ! Stefan-Boltzmann constant [W/(m^2 K^4)]

  real( kind = core_rknd ), parameter, public :: &
    T_freeze_K = 273.15_core_rknd ! Freezing point of water [K]

  ! Useful combinations of Rd and Rv
  real( kind = core_rknd ), parameter, public ::  & 
    ep  = Rd / Rv,    & ! ep  = 0.622_core_rknd  [-]
    ep1 = (1.0_core_rknd-ep)/ep,& ! ep1 = 0.61_core_rknd   [-]
    ep2 = 1.0_core_rknd/ep        ! ep2 = 1.61_core_rknd   [-]

  real( kind = core_rknd ), parameter, public :: & 
    kappa = Rd / Cp     ! kappa        [-]

  ! Changed g to grav to make it easier to find in the code 5/25/05
  ! real, parameter, public :: grav  = 9.80665_core_rknd ! Gravitational acceleration [m/s^2]
  real( kind = core_rknd ), parameter, public :: & 
    grav = 9.81_core_rknd, & ! Gravitational acceleration     [m/s^2]
    p0   = 1.0e5_core_rknd   ! Reference pressure             [Pa]

  ! Von Karman's constant
  ! Constant of the logarithmic wind profile in the surface layer
  real( kind = core_rknd ), parameter, public :: & 
    vonk   = 0.4_core_rknd,    & ! Accepted value is 0.40 (+/-) 0.01      [-]
    rho_lw = 1000.0_core_rknd    ! Density of liquid water                [kg/m^3]
    
#endif

  real( kind = core_rknd ), parameter, public :: & 
    rho_ice = 917.0_core_rknd    ! Density of ice      [kg/m^3]

  ! Tolerances below which we consider moments to be zero
  real( kind = core_rknd ), parameter, public ::  & 
    w_tol   = 2.e-2_core_rknd, & ! [m/s]
    thl_tol = 1.e-2_core_rknd, & ! [K]
    rt_tol  = max( 1.e-8_core_rknd, epsilon(pi) ), & ! [kg/kg]
    chi_tol = max( 1.e-8_core_rknd, epsilon(pi) ), & ! [kg/kg]
    eta_tol = chi_tol       ! [kg/kg]

  ! Tolerances for use by the monatonic flux limiter.
  ! rt_tol_mfl is larger than rt_tol. rt_tol is extremely small
  ! (1e-8) to prevent spurious cloud formation aloft in LBA.
  ! rt_tol_mfl is larger (1e-4) to prevent the mfl from
  ! depositing moisture at the top of the domain.
  real( kind = core_rknd ), parameter, public :: &
    thl_tol_mfl = 0.2_core_rknd, & ! [K]
    rt_tol_mfl = 1.e-4_core_rknd     ! [kg/kg]

  ! The tolerance for w'^2 is the square of the tolerance for w.
  real( kind = core_rknd ), parameter, public :: &
    w_tol_sqd = w_tol**2 ! [m^2/s^2]

  ! Set tolerances for Khairoutdinov and Kogan rain microphysics to insure
  ! against numerical errors.  The tolerance values for Nc, rr, and Nr insure
  ! against underflow errors in computing the PDF for l_kk_rain.  Basically,
  ! they insure that those values squared won't be less then 10^-38, which is
  ! the lowest number that can be numerically represented.  However, the
  ! tolerance value for rc doubles as the lowest mixing ratio there can be to
  ! still officially have a cloud at that level.  This is figured to be about
  ! 1.0_core_rknd x 10^-7 kg/kg.  Brian; February 10, 2007.
  real( kind = core_rknd ), parameter, public :: & 
    rc_tol  = 1.0E-6_core_rknd, & ! Tolerance value for r_c  [kg/kg]
    Nc_tol  = 1.0E+2_core_rknd, & ! Tolerance value for N_c  [#/kg]
    Ncn_tol = 1.0E+2_core_rknd    ! Tolerance value for N_cn [#/kg]

  real( kind = core_rknd ), parameter, public :: & 
    mvr_cloud_max = 1.6E-5_core_rknd    ! Max. avg. mean vol. rad. cloud    [m]

  real( kind = core_rknd ), parameter, public :: &
    Nc_in_cloud_min = 2.0e+4_core_rknd

  ! Precipitating hydrometeor tolerances for mixing ratios.
  real( kind = core_rknd ), parameter, public :: & 
    rr_tol = 1.0E-10_core_rknd, & ! Tolerance value for r_r [kg/kg]
    ri_tol = 1.0E-10_core_rknd, & ! Tolerance value for r_i [kg/kg]
    rs_tol = 1.0E-10_core_rknd, & ! Tolerance value for r_s [kg/kg]
    rg_tol = 1.0E-10_core_rknd    ! Tolerance value for r_g [kg/kg]

  ! Maximum allowable values for the average mean volume radius of the various
  ! hydrometeor species.
  real( kind = core_rknd ), parameter, public :: & 
    mvr_rain_max    = 5.0E-3_core_rknd, & ! Max. avg. mean vol. rad. rain    [m]
    mvr_ice_max     = 1.3E-4_core_rknd, & ! Max. avg. mean vol. rad. ice     [m]
    mvr_snow_max    = 1.0E-2_core_rknd, & ! Max. avg. mean vol. rad. snow    [m]
    mvr_graupel_max = 2.0E-2_core_rknd    ! Max. avg. mean vol. rad. graupel [m]

  ! Precipitating hydrometeor tolerances for concentrations.
  ! Tolerance value for N_r [#/kg]
  real( kind = core_rknd ), parameter, public :: & 
    Nr_tol = ( one / ( four_thirds * pi * rho_lw * mvr_rain_max**3 ) ) &
             * rr_tol

  ! Tolerance value for N_i [#/kg]
  real( kind = core_rknd ), parameter, public :: & 
    Ni_tol = ( one / ( four_thirds * pi * rho_ice * mvr_ice_max**3 ) ) &
             * ri_tol

  ! Tolerance value for N_s [#/kg]
  real( kind = core_rknd ), parameter, public :: & 
    Ns_tol = ( one / ( four_thirds * pi * rho_ice * mvr_snow_max**3 ) ) &
             * rs_tol

  ! Tolerance value for N_s [#/kg]
  real( kind = core_rknd ), parameter, public :: & 
    Ng_tol = ( one / ( four_thirds * pi * rho_ice * mvr_graupel_max**3 ) ) &
             * rg_tol

  ! Minimum value for em (turbulence kinetic energy)
  ! If anisotropic TKE is enabled, em = (1/2) * ( up2 + vp2 + wp2 );
  ! otherwise, em = (3/2) * wp2.  Since up2, vp2, and wp2 all have
  ! the same minimum threshold value of w_tol_sqd, em cannot be less
  ! than (3/2) * w_tol_sqd.  Thus, em_min = (3/2) * w_tol_sqd.
  real( kind = core_rknd ), parameter, public :: &
    em_min = 1.5_core_rknd * w_tol_sqd  ! [m^2/s^2]

  real( kind = core_rknd ), parameter, public ::  & 
    eps = 1.0e-10_core_rknd ! Small value to prevent a divide by zero

  real( kind = core_rknd ), parameter, public :: &
    max_num_stdevs = 5.0_core_rknd ! Range of standard deviations for statistical significance

  real( kind = core_rknd ), parameter, public ::  &
    zero_threshold = 0.0_core_rknd ! Defining a threshold on a physical quantity to be 0.

  ! The maximum absolute value (or magnitude) that a correlation is allowed to
  ! have.  Statistically, a correlation is not allowed to be less than -1 or
  ! greater than 1, so the maximum magnitude would be 1.
  real( kind = core_rknd ), parameter, public :: &
    max_mag_correlation = 0.99_core_rknd, &      ! Most correlations use this
    max_mag_correlation_flux = 0.99_core_rknd    ! Special for wprtp and wpthlp

  ! "base" smoothing magnitude before scaling for the respective data structure.
  ! See https://github.com/larson-group/clubb/issues/965#issuecomment-1119816722
  ! for a plot on how output behaves with varying min_max_smth_mag
  real( kind = core_rknd ), parameter, public :: &
    min_max_smth_mag = 1.0e-9_core_rknd

  real( kind = core_rknd ), parameter, public :: &
    cloud_frac_min = 0.005_core_rknd ! Threshold for cloud fractions

  real( kind = core_rknd ), parameter, public :: &
    wp2_max = 1000._core_rknd ! Max value for variance clipping

  !-----------------------------------------------------------------------------
  ! Useful conversion factors.
  !-----------------------------------------------------------------------------
  real(kind=core_rknd), parameter, public ::  & 
    sec_per_day = 86400.0_core_rknd, & ! Seconds in a day.
    sec_per_hr  = 3600.0_core_rknd,  & ! Seconds in an hour.
    sec_per_min = 60.0_core_rknd,    & ! Seconds in a minute.
    min_per_hr = 60.0_core_rknd        ! Minutes in an hour.

  real( kind = core_rknd ), parameter, public :: & 
    g_per_kg = 1000.0_core_rknd     ! Grams in a kilogram.

  real( kind = core_rknd ), parameter, public :: &
    pascal_per_mb = 100.0_core_rknd ! Pascals per Millibar

  real( kind = core_rknd ), parameter, public :: & 
    cm3_per_m3   = 1.e6_core_rknd,  & ! Cubic centimeters per cubic meter
    micron_per_m = 1.e6_core_rknd,  & ! Micrometers per meter
    cm_per_m     = 100._core_rknd,  & ! Centimeters per meter
    mm_per_m     = 1000._core_rknd    ! Millimeters per meter

  !-----------------------------------------------------------------------------
  ! Unused variable
  !-----------------------------------------------------------------------------
  real( kind = core_rknd ), parameter, public :: &
    unused_var = -999._core_rknd      ! The standard value for unused variables

!=============================================================================

end module constants_clubb
